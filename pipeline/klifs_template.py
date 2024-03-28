import sys
import numpy as np
import requests as req
from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
import pandas as pd
import tqdm

fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=3, fpSize=2048)
tanimoto = lambda a, b: (a & b).sum() / (a.sum() + b.sum() - (a & b).sum())


@lru_cache(maxsize=100000)
def get_kinase_IDs(uniprot_id):
    resp = req.get(
        "https://klifs.net/api/kinase_ID", params={"kinase_name": uniprot_id}
    )
    resp.raise_for_status()
    result = resp.json()
    return tuple(c["kinase_ID"] for c in result)


@lru_cache(maxsize=100000)
def ligands_list(kinase_IDs):
    ligands = list()
    for kinase_ID in kinase_IDs:
        try:
            ligands.extend(ligands_list_single_kinase(kinase_ID))
        except:
            # print('KLIFS inconsistent for kinase_ID', kinase_ID)
            pass  # inconsistencies in KLIFS
    return ligands


@lru_cache(maxsize=100000)
def ligands_list_single_kinase(kinase_ID):
    resp = req.get(
        "https://klifs.net/api/ligands_list", params={"kinase_ID": kinase_ID}
    )
    resp.raise_for_status()
    return resp.json()


@lru_cache(maxsize=100000)
def get_best_structure(kinase_IDs, ligand_ID):
    resp = req.get(
        "https://klifs.net/api/ligands_list_structures", params={"ligand_ID": ligand_ID}
    )
    resp.raise_for_status()
    return max(
        [struc for struc in resp.json() if struc["kinase_ID"] in kinase_IDs],
        key=(lambda x: x["quality_score"]),
    )["structure_ID"]


@lru_cache(maxsize=100000)
def get_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return None
    try:
        return fpgen.GetFingerprintAsNumPy(mol)
    except:
        return None


def get_template(kinase_ID, ligand_smiles):
    ligs = ligands_list(kinase_ID)
    ref_ligand_fp = get_fp(ligand_smiles)
    if ref_ligand_fp is None:
        raise ValueError("Illegal reference ligand")
    lig_fps = [get_fp(l["SMILES"]) for l in ligs]
    lig_fps = [l for l in lig_fps if l is not None]
    if len(lig_fps) == 0:
        raise ValueError("No valid template smiles")
    similarities = [tanimoto(ref_ligand_fp, fp) for fp in lig_fps]
    best_ligand_idx = np.argmax(similarities)
    best_ligand_ID = ligs[best_ligand_idx]["ligand_ID"]
    return (
        get_best_structure(tuple(kinase_ID), best_ligand_ID),
        similarities[best_ligand_idx],
    )


def main(
    kinodata_file="../data/activities-chembl33_v0.5.csv", template_file="templates.csv"
):
    kinodata = pd.read_csv(kinodata_file, index_col=0)

    templates = pd.read_csv(template_file)

    with open(template_file, "a") as f:
        f.write(
            "activities.activity_id,similar.klifs_structure_id,similar.fp_similarity\n"
        )
        structure_IDs = []
        for activity_id, row in tqdm.tqdm(
            kinodata.set_index("activities.activity_id").iterrows(), total=len(kinodata)
        ):
            if activity_id in templates["activities.activity_id"].values:
                continue
            try:
                kinase_IDs = get_kinase_IDs(row["UniprotID"])
            except:
                structure_ID = -2
                structure_IDs.append(structure_ID)
                f.write(f"{activity_id},{structure_ID},nan\n")
                continue
            lig_smiles = row["compound_structures.canonical_smiles"]
            try:
                structure_ID = get_template(kinase_IDs, lig_smiles)
            except:
                structure_ID = -1
                f.write(f"{activity_id},{structure_ID},nan\n")
                continue
            f.write(f"{activity_id},{structure_ID[0]},{structure_ID[1]}\n")
            structure_IDs.append(structure_ID[0])


if __name__ == "__main__":
    kinodata_file = sys.argv[1]
    template_file = sys.argv[2]
    main(kinodata_file, template_file)
