from importlib import resources
import inspect
from pathlib import Path
import time, sys, os, shutil, tempfile
import socket

from kinoml.core.ligands import Ligand
from kinoml.core.proteins import Protein
from kinoml.core.systems import ProteinLigandComplex
from kinoml.features.complexes import (
    MostSimilarPDBLigandFeaturizer,
    OEDockingFeaturizer,
)

from rdkit import Chem
import pandas as pd
import requests as req
from multiprocessing import Pool, cpu_count
import traceback
import tqdm
import numpy as np
import MDAnalysis as mda
import pathlib


def main():
    print(f"HOST:{socket.gethostname()} - PID:{os.getpid()}")
    HERE = Path(".").absolute()
    activity_id = int(sys.argv[1])
    pdb_filepath = sys.argv[2]
    smiles = sys.argv[3]
    output_dir = pathlib.Path(sys.argv[4])

    protein = Protein.from_file(pdb_filepath)
    ligand = Ligand(smiles=smiles)
    system = ProteinLigandComplex(components=[protein, ligand])

    featurizer = OEDockingFeaturizer(
        output_dir=output_dir, method="Posit", use_multiprocessing=False
    )
    print("start featurization")
    start_time = time.time()
    system = featurizer.featurize([system])[0]
    duration = time.time() - start_time

    print("write result")
    universe = system.featurizations["last"]
    docking_score = universe._topology.docking_score
    posit_probability = universe._topology.posit_probability
    with open(output_dir / "docking.csv", "a") as f:
        f.write(
            ",".join(
                list(
                    map(str, [activity_id, docking_score, posit_probability, duration])
                )
            )
            + "\n"
        )

    print("write", output_dir / f"{activity_id}_complex.pdb")
    with mda.coordinates.PDB.PDBWriter(output_dir / f"{activity_id}_ligand.pdb") as w:
        w.write(universe.select_atoms("resname LIG"))
    print("done")


if __name__ == "__main__":
    main()
