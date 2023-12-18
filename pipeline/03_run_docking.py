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
        print('run', self.ident)
        self.process = subprocess.Popen(
            [
                "conda",
                "run",
                "--no-capture-output",
                "-n",
                "kinodata-3D",
                "python",
                "pipeline/docking.py",
                str(self.ident),
                str(self.protein_filepath),
                str(self.smiles),
                str(self.output_dir),
            ],
            # stdout=out,
            # stderr=err,
            close_fds=True,
            shell=False,
        )

    @property
    def output_dir(self):
        out_dir = HERE / "cache" / "complexes" / str(self.ident)
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
        output_dir=HERE / "cache",
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
        return self.output_dir / "docking_failures.csv"

    @property
    def success_file(self):
        return self.output_dir / "docking_successes.csv"

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
            if psutil.virtual_memory().free / 1024 ** 3 <= self.total_mem_start_limit:
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
        kinodata_path="data/activities-chembl31.csv.gz",
        similar_pdb_path="data/most_similar.csv.gz",
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


def prepare_tasks(
    data: TemplateData, output_dir=HERE / "data"
) -> List[DockingTask]:
    klifs_structures_file = output_dir / "klifs_structures.csv.gz"
    structures = pd.read_csv(
        klifs_structures_file, index_col="activities.activity_id"
    )
    print("-> populate waitlist")
    tasks, idents, ids = list(), list(), list()
    for ident, row in data.similar_pdbs.iloc[:100].iterrows():
        if structures is None or ident not in structures.index:
            continue
        else:
            structure_ID = int(structures.loc[ident, "similar.klifs_structure_id"])
        protein_file = HERE / "cache" / f"{structure_ID}.pdb"
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


def main_docking():
    print("-> read data")
    data = TemplateData()
    print("-> prepare docking tasks")
    tasks = prepare_tasks(data)

    scheduler = Scheduler(
        tasks,
        DockingJob,
        capacity=64,
        proc_mem_limit=10,
        timeout=10,
        total_mem_start_limit=64,
        output_dir=HERE / "cache",
    )

    print("-> start docking")
    scheduler.run()


if __name__ == "__main__":
    main_docking()
