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
    # try:
    protein = Protein(uniprot_id=uniprot_id, toolkit="MDAnalysis")
    ligand = Ligand(smiles=ligand_smiles)
    system = ProteinLigandComplex(components=[protein, ligand])

    return system


def get_most_similar(args):
    featurizer = MostSimilarPDBLigandFeaturizer(
        similarity_metric="fingerprint",
        cache_dir=CACHE_DIR,
    )

    system = featurizer.featurize([system])[0]

    ligand_id = system.protein.expo_id
    pdb_id = system.protein.pdb_id
    chain = system.protein.chain_id

    return ligand_id, pdb_id, chain

if __name__ == '__main__':
    kinodata = pd.read_csv("activity.csv", index_col=0)
    done = pd.read_csv(MOST_SIMILAR)['activities.activity_id'].values
    kinodata = kinodata[~kinodata['activities.activity_id'].isin(done)]

    print('Setting up kinoml systems')
    systems = list()
    for _, row in kinodata.iterrows():
        uniprot_id = row["UniprotID"]
        ligand_smiles = row["compound_structures.canonical_smiles"]
        ident = row["activities.activity_id"]
        systems.append(get_system((ident, uniprot_id, ligand_smiles)))


    print('Finding most similar PDBs')
    CHUNKSIZE = 1
#   with open(MOST_SIMILAR, "w") as f:
#       f.write("activities.activity_id,similar.ligand_pdb,similar.complex_pdb,similar.chain\n")
    featurized_systems = list()
    for i in tqdm.tqdm(range(0, len(systems), CHUNKSIZE)):
        try:
            print('batch', i / CHUNKSIZE)
            featurizer = MostSimilarPDBLigandFeaturizer(
                similarity_metric="fingerprint",
                cache_dir=CACHE_DIR,
                n_processes  = CHUNKSIZE
            )
            featurized_systems = featurizer.featurize(systems[i:i+CHUNKSIZE])
            idents = kinodata['activities.activity_id'].values[i:i+CHUNKSIZE]
            print(' writing results')
            for ident, system in zip(idents, featurized_systems):
                print(ident)
                ligand_id = system.protein.expo_id
                pdb_id = system.protein.pdb_id
                chain = system.protein.chain_id
                with open(MOST_SIMILAR, "a") as f:
                    f.write(",".join(map(str, [ident, ligand_id, pdb_id, chain])) + "\n")
        except:
            print('batch', i, 'failed')
