from pathlib import Path
import time, sys, os, glob
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
    ligand = Ligand(smiles=smiles, name='')
    system = ProteinLigandComplex(components=[protein, ligand])

    featurizer = OEDockingFeaturizer(
        output_dir=output_dir, method="Posit", use_multiprocessing=False, all_poses=True,
    )
    print("start featurization")
    start_time = time.time()
    system = featurizer.featurize([system])[0]
    duration = time.time() - start_time

    print("write result")

    all_files = glob.glob(os.path.join(output_dir, '*'))
    for file_path in all_files:
        if not file_path.endswith('_ligand.sdf'):
            os.remove(file_path)

    universe = system.featurizations["last"]
    with open(output_dir / "docking.csv", "a") as f:
        f.write(
            ",".join(
                list(
                    map(str, [activity_id, duration])
                )
            )
            + "\n"
        )

    print("done")


if __name__ == "__main__":
    main()
