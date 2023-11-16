from typing import List
import pickle
import shutil
import requests as req
import subprocess
import pathlib
import os
import multiprocessing
import time
import psutil
import sys
import signal
import random
from functools import cached_property
from collections import namedtuple

import pandas as pd

HERE = pathlib.Path(".").absolute()

DockingTask = namedtuple("DockingTask", "ident protein smiles")
SimilarKLIFSTask = namedtuple("SimilarKLIFSTask", "ident uniprot_id smiles")


class Job:
    def __init__(self):
        self.start_time = time.time()
        if not self.success:
            self.start()

    def start(self):
        """
        This method is to be implemented by child classes. It needs to set
        `self.process` to a system process doing the task.
        """
        raise NotImplementedError("Subclass must implement")

    @property
    def success(self):
        """
        This method is to be implemented by child classes.
        It returns whether the job has already done (maybe previously).
        The job is not started if `self.success == True` during initialization.
        """
        raise NotImplementedError("Subclass must implement")

    @property
    def pid(self):
        if hasattr(self, "process"):
            return self.process.pid
        else:
            return -1

    @property
    def running(self):
        return self.process.poll()

    @property
    def memory(self):
        """overall process tree's memory [gb]"""
        try:
            return sum(
                child.memory_info().rss / 1024**3
                for child in psutil.Process(self.pid).children(recursive=True)
            )
        except psutil.NoSuchProcess:
            return 0

    @property
    def runtime(self):
        """runtime in minutes"""
        return (time.time() - self.start_time) / 60

    def suicide(
        self, sig=signal.SIGTERM, include_parent=True, timeout=None, on_terminate=None
    ):
        """Kill a process tree (including grandchildren) with signal
        "sig" and return a (gone, still_alive) tuple.
        "on_terminate", if specified, is a callback function which is
        called as soon as a child terminates.
        """
        assert self.pid != os.getpid(), "I won't kill the master"
        if self.pid < 0:
            return
        try:
            parent = psutil.Process(self.pid)
        except psutil.NoSuchProcess:
            return
        children = parent.children(recursive=True)
        if include_parent:
            children.append(parent)
        for p in children:
            try:
                p.send_signal(sig)
            except psutil.NoSuchProcess:
                pass
        gone, alive = psutil.wait_procs(
            children, timeout=timeout, callback=on_terminate
        )
        return (gone, alive)


class SimilarKLIFSJob(Job):
    def __init__(self, task: SimilarKLIFSTask):
        self.ident = task.ident
        self.smiles = task.smiles
        self.uniprot_id = task.uniprot_id
        super().__init__()

    def start(self):
        self.process = subprocess.Popen(
            [
                "conda",
                "run",
                "--no-capture-output",
                "-n",
                "kinoml",
                "python",
                "run_docking.py",
                str(self.ident),
                str(self.uniprot_id),
                str(self.smiles),
                str(self.output_file),
            ],
            close_fds=True,
            shell=False,
        )

    @property
    def output_file(self):
        return self.output_dir / f"{self.ident}.csv"

    @property
    def output_dir(self):
        out_dir = HERE / "docking_pipeline" / "KLIFS"
        out_dir.mkdir(exist_ok=True, parents=True)
        return out_dir

    @property
    def success(self):
        return self.output_file.exists()


class DockingJob(Job):
    def __init__(self, task):
        self.ident = task.ident
        self.protein_filepath = self.output_dir / "protein.pdb"
        shutil.copy2(task.protein, self.protein_filepath)
        self.smiles = task.smiles
        # maybe the docking was done in a previous run
        super().__init__()

    def start(self):
        out = open(self.output_dir / "run.log", "w")
        err = open(self.output_dir / "run.err", "w")
        self.process = subprocess.Popen(
            [
                "conda",
                "run",
                "--no-capture-output",
                "-n",
                "kinoml",
                "python",
                "docking.py",
                str(self.ident),
                str(self.protein_filepath),
                str(self.smiles),
                str(self.output_dir),
            ],
            stdout=out,
            stderr=err,
            close_fds=True,
            shell=False,
        )

    @property
    def output_dir(self):
        out_dir = HERE / "docking_pipeline" / "complexes" / str(self.ident)
        out_dir.mkdir(exist_ok=True, parents=True)
        return out_dir

    @property
    def success(self):
        """Check if docking results are there for `ident`."""
        output_file = self.output_dir / "docking.csv"
        ligand_file = self.output_dir / f"{self.ident}_ligand.pdb"
        success = output_file.exists() and ligand_file.exists()
        return success


class Scheduler:
    def __init__(
        self,
        tasks,
        jobtype,
        capacity=os.cpu_count(),
        proc_mem_limit=10,
        timeout=10,
        total_mem_start_limit=50,
        output_dir=HERE / "docking_pipeline",
    ):
        """
        Parameters
        ----------
        tasks: List[Task]
            docking tasks
        jobtype:
            the job type to process tasks
        capacity: int
            the maximum number of concurrent docking jobs
        proc_mem_limit: int
            memory limit per job in gb
        timeout: int
            job timeout in minutes
        total_mem_start_limit: int
            limit in percentage on memory above which no new jobs are started
        """
        self.capacity = capacity
        self.jobtype = jobtype
        self.proc_mem_limit = proc_mem_limit
        self.timeout = timeout
        self.total_mem_start_limit = total_mem_start_limit
        self.running = list()
        self.waitlist = tasks
        random.shuffle(self.waitlist)
        self.output_dir = output_dir
        self.done = self.try_reading_done()

    def try_reading_done(self):
        done = []
        if (self.success_file).exists():
            done += list(
                pd.read_csv(self.success_file, names=["ident"]).values.flatten()
            )
        if (self.failure_file).exists():
            done += list(
                pd.read_csv(self.failure_file, names="ident reason".split())[
                    "ident"
                ].values
            )
        return done

    @property
    def failure_file(self):
        return self.output_dir / "failures.csv"

    @property
    def success_file(self):
        return self.output_dir / "success.csv"

    def print_status(self):
        print(
            f"|waitlist| = {len(self.waitlist)} |running| = {len(self.running)}",
        )

    def run(self):
        while len(self.waitlist) > 0 or len(self.running) > 0:
            self.print_status()
            time.sleep(1)  # busy wait...
            self.cleanup_running()

            # check overall memory usage
            if psutil.virtual_memory()[2] > self.total_mem_start_limit:
                continue

            self.start_dockings()

    def start_dockings(self):
        for _ in range(self.capacity - len(self.running)):
            if len(self.waitlist) == 0:
                return
            task = self.waitlist.pop()
            if task.ident in self.done:
                continue
            self.running.append(self.jobtype(task))

    def log_fail(self, ident, reason):
        with open(self.failure_file, "a") as f:
            f.write(f"{ident},{reason}\n")

    def log_success(self, ident):
        with open(self.success_file, "a") as f:
            f.write(f"{ident}\n")

    def cleanup_running(self):
        # clean self.running processes
        still_running = list()
        for i, job in enumerate(self.running):
            # check for completion
            if job.success:
                self.log_success(job.ident)
                job.suicide()  # make sure it's dead
                continue

            if not job.running and not job.success and job.runtime > 10:
                self.log_fail(job.ident, "death")
                job.suicide()
                continue

            # check for timeout and out-of-memory
            if job.memory > self.proc_mem_limit:
                self.log_fail(job.ident, "memory")
                job.suicide()
                continue
            if job.runtime > self.timeout:
                self.log_fail(job.ident, "timeout")
                job.suicide()
                continue
            still_running.append(i)
        self.running = [job for i, job in enumerate(self.running) if i in still_running]


class TemplateData:
    def __init__(
        self,
        kinodata_path="activities-chembl31.csv",
        similar_pdb_path="docking_pipeline/most_similar.csv",
    ):
        self.kinodata_path = kinodata_path
        self.similar_pdb_path = similar_pdb_path

    @cached_property
    def kinodata(self):
        # activities.activity_id,assays.chembl_id,target_dictionary.chembl_id,molecule_dictionary.chembl_id,molecule_dictionary.max_phase,activities.standard_type,activities.standard_value,activities.standard_units,compound_structures.canonical_smiles,compound_structures.standard_inchi,component_sequences.sequence,assays.confidence_score,docs.chembl_id,docs.year,docs.authors,UniprotID
        return pd.read_csv(self.kinodata_path, index_col="activities.activity_id")

    @cached_property
    def similar_pdbs(self):
        # activities.activity_id,similar.ligand_pdb,similar.complex_pdb,similar.chain
        return pd.read_csv(self.similar_pdb_path, index_col="activities.activity_id")


def prepare_klifs_tasks(
    data: TemplateData, output_dir=HERE / "docking_pipeline"
) -> List[DockingTask]:
    tasks = list()
    for ident, row in data.kinodata.iterrows():
        similar_file = output_dir / "KLIFS" / f"{ident}.csv"
        cached_tasks = output_dir / "KLIFS" / f"{ident}.tasks"
        if cached_tasks.exists():
            print(f"{ident} cached")
            with open(cached_tasks, "rb") as f:
                ident_tasks = pickle.load(f)
        elif similar_file.exists():
            ident_tasks = []
            similars = pd.read_csv(similar_file)
            for i, similar in similars.iterrows():
                structure_ID = get_klifs_structure_id(
                    similar["pdb_id"],
                    {
                        "chain": similar["chain_id"],
                        "DFG": similar["dfg"],
                        "aC_helix": similar["ac_helix"],
                        "ligand": similar["expo_id"],
                    },
                )
                protein_file = get_klifs_protein(
                    structure_ID,
                    output_path=HERE / "docking_pipeline" / "KLIFS_proteins",
                )
                task = DockingTask(
                    f"{ident}.{i}",
                    protein_file,
                    row["compound_structures.canonical_smiles"],
                )
                print(task)
                ident_tasks.append(tasks.append(task))
                with open(cached_tasks, "wb") as f:
                    pickle.dump(ident_tasks, f)
        tasks.extend(ident_tasks)
    return tasks


def prepare_tasks(
    data: TemplateData, output_dir=HERE / "docking_pipeline"
) -> List[DockingTask]:
    klifs_structures_file = output_dir / "klifs_structures.csv"
    if klifs_structures_file.exists():
        structures = pd.read_csv(
            klifs_structures_file, index_col="activities.activity_id"
        )
    else:
        structures = None
    print("-> populate waitlist")
    tasks, idents, ids = list(), list(), list()
    for ident, row in data.similar_pdbs.iterrows():
        if structures is None or ident not in structures.index:
            continue
        #     try:
        #         structure_ID = get_klifs_structure_id(
        #             row["similar.complex_pdb"],
        #             {
        #                 "chain": row["similar.chain"],
        #                 "ligand": row["similar.ligand_pdb"],
        #             },
        #         )
        #     except ValueError:
        #         with open(output_dir / "klifs_failures.csv", "a") as f:
        #             f.write(f"{ident},no klifs structure\n")
        #         continue
        #     except req.exceptions.HTTPError:
        #         with open(output_dir / "klifs_failures.csv", "a") as f:
        #             f.write(f"{ident},bad request\n")
        #         continue
        #     idents.append(ident)
        #     ids.append(int(structure_ID))
        else:
            structure_ID = int(structures.loc[ident, "similar.klifs_structure_id"])
        protein_file = get_klifs_protein(
            structure_ID, output_path=HERE / "docking_pipeline" / "KLIFS_proteins"
        )
        task = DockingTask(
            ident,
            protein_file,
            data.kinodata.loc[ident, "compound_structures.canonical_smiles"],
        )
        tasks.append(task)
    structure_ids = pd.DataFrame(
        {"activities.activity_id": idents, "similar.klifs_structure_id": ids}
    ).set_index("activities.activity_id")
    if structures is not None:
        structure_ids = pd.concat([structure_ids, structures])

    structure_ids.to_csv(klifs_structures_file)

    return tasks


def get_klifs_structure_id(pdb_id: str, features: dict()):
    """
    Get the complex PDB from KLIFS. Tie-braking via quality score.

    Parameters
    ----------
    pdb_id: str
        PDB ID.
    features: dict
        other fields to match in the klifs query

    Returns
    -------
    structure_ID: int
        KLIFS structure id for the most similar
    """
    resp = req.get(
        "https://klifs.net/api_v2/structures_pdb_list", {"pdb-codes": pdb_id}
    )
    klifs_info = None
    resp.raise_for_status()
    for info in resp.json():
        if all(str(info[k]).upper() == str(v).upper() for k, v in features.items()):
            if klifs_info is None:
                klifs_info = info
            elif klifs_info["quality_score"] < info["quality_score"]:
                klifs_info = info
    if klifs_info is None:
        raise ValueError(f"not found pdb:{pdb_id} with props: {features}")
    return klifs_info["structure_ID"]


def get_klifs_protein(structure_ID: int, output_path=HERE):
    """
    Get the complex mol2 from KLIFS.

    Parameters
    ----------
    structure_ID: int
        KLIFS structure ID
    path: Path, optional
        folder to store the structure in.

    Returns
    -------
    file_path: Path
        path of the mol2 file.
    """
    filename = output_path / f"{structure_ID}.pdb"
    if filename.exists():
        return filename
    pathlib.Path(output_path).mkdir(exist_ok=True, parents=True)
    resp = req.get(
        "https://klifs.net/api_v2/structure_get_pdb_complex",
        {"structure_ID": structure_ID},
    )

    with open(filename, "w") as f:
        f.write(resp.text)

    return filename


def prepare_similar_klifs_tasks(data):
    print("Setting up kinoml systems")
    systems = list()
    done = list()
    for ident, row in data.kinodata.iterrows():
        uniprot_id = row["UniprotID"]
        ligand_smiles = row["compound_structures.canonical_smiles"]
        systems.append(SimilarKLIFSTask(ident, uniprot_id, ligand_smiles))
    return systems


def main_similar_klifs():
    print("-> read data")
    data = TemplateData()
    print("-> prepare docking tasks")
    tasks = prepare_similar_klifs_tasks(data)

    scheduler = Scheduler(
        tasks,
        SimilarKLIFSJob,
        capacity=256,
        proc_mem_limit=10,
        timeout=10,
        total_mem_start_limit=30,
        output_dir=HERE / "docking_pipeline" / "KLIFS",
    )

    print("-> start similar KLIFS search")
    scheduler.run()


def main_docking():
    print("-> read data")
    data = TemplateData()
    print("-> prepare docking tasks")
    tasks = prepare_tasks(data)

    scheduler = Scheduler(
        tasks,
        DockingJob,
        capacity=128,
        proc_mem_limit=10,
        timeout=10,
        total_mem_start_limit=30,
        output_dir=HERE / "docking_pipeline",
    )

    print("-> start docking")
    scheduler.run()


if __name__ == "__main__":
    # main_similar_klifs()
    main_docking()
