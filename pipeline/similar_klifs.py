from importlib import resources
import sys
import inspect
from pathlib import Path
import traceback
import tempfile, os, shutil, time
import signal

from kinoml.core.ligands import Ligand
from kinoml.core.proteins import Protein, KLIFSKinase
from kinoml.core.systems import ProteinLigandComplex
from kinoml.features.complexes import MostSimilarPDBLigandFeaturizer, KLIFSConformationTemplatesFeaturizer, OEDockingFeaturizer

from rdkit import Chem
import pandas as pd
import requests as req
from multiprocessing import Pool, cpu_count
import tqdm
import numpy as np
import dask
import dask.dataframe as dd
import MDAnalysis as mda

# directories
HERE = Path(".").absolute()
CACHE_DIR = HERE / "docking_pipeline"
MOST_SIMILAR = CACHE_DIR / "most_similar.csv"
KLIFS_DIR = CACHE_DIR / "KLIFS"
KLIFS_MAP = CACHE_DIR / "similar_klifs_structures.csv"
TEMP_DIR = CACHE_DIR / "temporary"
DOCKING_DIR = CACHE_DIR / "docking"

def get_system(args):
    ident, uniprot_id, ligand_smiles = args
    protein = KLIFSKinase(uniprot_id=uniprot_id, toolkit="MDAnalysis")
    ligand = Ligand(smiles=ligand_smiles)
    system = ProteinLigandComplex(components=[protein, ligand])

    return system

def output_file(ident):
    return KLIFS_DIR / f'{ident}.csv'

if __name__ == '__main__':
    kinodata = pd.read_csv("activities-chembl31.csv", index_col=0)

    print('Setting up kinoml systems')
    systems = list()
    done = list()
    for _, row in kinodata.iterrows():
        ident = row["activities.activity_id"]
        if output_file(ident).exists():
            done.append(ident)
            continue
        uniprot_id = row["UniprotID"]
        ligand_smiles = row["compound_structures.canonical_smiles"]
        systems.append(get_system((ident, uniprot_id, ligand_smiles)))
    kinodata = kinodata[~kinodata['activities.activity_id'].isin(done)]
    assert len(kinodata) == len(systems), (len(kinodata), len(systems))


    print('Finding most similar PDBs')
    CHUNKSIZE = 8
#     with open(MOST_SIMILAR, "w") as f:
#         f.write("activities.activity_id,similar.ligand_pdb,similar.complex_pdb,similar.chain\n")
    featurized_systems = list()
    for i in tqdm.tqdm(range(0, len(systems), CHUNKSIZE)):
        try:
            print('batch', i / CHUNKSIZE)
            featurizer =KLIFSConformationTemplatesFeaturizer(
                similarity_metric="fingerprint",
                cache_dir=CACHE_DIR,
            )
            featurized_systems = featurizer.featurize(systems[i:i+CHUNKSIZE])
            idents = kinodata['activities.activity_id'].values[i:i+CHUNKSIZE]
        except:
            print('batch', i, 'failed')
            continue
        print(' writing results')
        for ident, system in zip(idents, featurized_systems):
            #try:
            filename = output_file(ident)
            print('result to', filename)
            system.featurizations['last'].to_csv(filename)
            #except:
            #    print('failed writing', KLIFS_DIR / f'{ident}.csv')
