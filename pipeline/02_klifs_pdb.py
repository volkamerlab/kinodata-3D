"""
Download KLIFS structures for a list of similar PDBs determind by `kinoml`.
"""
from pathlib import Path

import pandas as pd
import requests as req
import numpy as np


def klifs_pdb(
    pdb_id: str,
    chain: str,
    ligand_id: str,
    output_path=Path(".").absolute(),
    lazy=False,
):
    """
    Get the complex PDB from KLIFS.

    Parameters
    ----------
    pdb_id: str
        PDB ID.
    chain: str
        the chain.
    ligand_id: str
        ligand PDB ID
    path: Path, optional
        folder to store the structure in.

    Returns
    -------
    file_path: Path
        path of the PDB file.
    """
    resp = req.get(
        "https://klifs.net/api_v2/structures_pdb_list", {"pdb-codes": pdb_id}
    )
    klifs_info = None
    resp.raise_for_status()
    for info in resp.json():
        if (
            str(info["chain"]).upper() == str(chain).upper()
            and str(info["ligand"]).upper() == str(ligand_id).upper()
        ):
            if klifs_info is None:
                klifs_info = info
            elif klifs_info["quality_score"] < info["quality_score"]:
                klifs_info = info
    if klifs_info is None:
        raise ValueError(f"not found pdb:{pdb_id} chain:{chain} lig_pdb:{ligand_id}")
    structure_ID = klifs_info["structure_ID"]
    filename = output_path / f"{structure_ID}.pdb"
    if lazy and filename.exists():
        return structure_ID
    resp = req.get(
        "https://klifs.net/api_v2/structure_get_pdb_complex",
        {"structure_ID": structure_ID},
    )

    with open(filename, "w") as f:
        f.write(resp.text)

    return structure_ID


if __name__ == "__main__":
    cache_dir = Path(".").absolute() / "data"
    cache_dir.mkdir(exist_ok=True)

    print("read similar.csv")
    df = pd.read_csv(
        cache_dir / "most_similar.csv.gz",
        index_col=0,
        names="ident ligand_pdb pdb chain similarity fp_sim".split(),
    )

    BATCH_SIZE = 1024

    batch_count = int(len(df) / BATCH_SIZE)

    for batch_num in range(0, batch_count):
        batch_start = batch_num * BATCH_SIZE
        batch_end = min((batch_num + 1) * BATCH_SIZE, len(df))

        print(batch_start, batch_end)
        print(f"get pdb batch {batch_num}/{batch_count}")
        with open(f"pdbs.csv", "a") as f:
            for idx, row in df[batch_start:batch_end].iterrows():
                print(f"download {idx}")
                try:
                    filename = klifs_pdb(
                        row["pdb"],
                        row["chain"],
                        row["ligand_pdb"],
                        output_path=cache_dir,
                        lazy=True,
                    )
                    f.write(f"{idx},{filename}\n")
                except Exception as e:
                    with open("download_err.csv", "a") as err:
                        err.write(f"{idx},{str(e)}\n")
