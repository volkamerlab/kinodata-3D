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
from pathlib import Path
import tqdm

import pandas as pd

HERE = pathlib.Path(".").absolute()
TEMPLATE_PATH = HERE / "data" / "templates"

DockingTask = namedtuple("DockingTask", "ident protein smiles")


class Job:
    def __init__(self):
        self.start_time = time.time()
        if not self.success and not self.failed:
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
        It returns whether the job has already been done (maybe previously).
        The job is not started if `self.success == True` during initialization.
        """
        raise NotImplementedError("Subclass must implement")

    @property
    def failed(self):
        """
        This method is to be implemented by child classes.
        It returns whether the job has already failed.
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
        if self.pid < 0:
            return False
        return self.process.poll()

    @property
    def memory(self):
        """overall process tree's memory [gb]"""
        if self.pid < 0:
            return 0
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
        self, sig=signal.SIGKILL, include_parent=True, timeout=None, on_terminate=None
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


class DockingJob(Job):
    def __init__(self, task):
        self.ident = task.ident
        self.protein_filepath = self.output_dir / "protein.pdb"
        self.smiles = task.smiles
        self.protein = task.protein
        # maybe the docking was done in a previous run
        super().__init__()

    def start(self):
        # if self.protein_filepath.exists():
        #    return
        shutil.copy2(self.protein, self.protein_filepath)
        out = open(self.output_dir / "run.log", "w")
        err = open(self.output_dir / "run.err", "w")
        self.process = subprocess.Popen(
            [
                "conda",
                "run",
                "--no-capture-output",
                "-n",
                "kinodata-3D",
                "python",
                "docking.py",
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
    def in_progress(self):
        return self.output_dir_name.exists()

    @property
    def output_dir_name(self):
        return HERE / "data" / "cache" / "complexes" / str(self.ident)

    @property
    def output_dir(self):
        self.output_dir_name.mkdir(exist_ok=True, parents=True)
        return self.output_dir_name

    @property
    def success(self):
        """Check if docking results are there for `ident`."""
        output_file = self.output_dir / "docking.csv"
        ligand_file = self.output_dir / f"{self.ident}_ligand.pdb"
        success = output_file.exists() and ligand_file.exists()
        return success

    @property
    def failed(self):
        return (self.output_dir / "fail").exists()

    def fail(self, reason):
        (self.output_dir / "fail").touch()
        (self.output_dir / "fail").write_text(reason)


class Scheduler:
    def __init__(
        self,
        tasks,
        jobtype,
        capacity=os.cpu_count(),
        proc_mem_limit=10,
        timeout=10,
        total_mem_start_limit=50,
        output_dir=HERE / "data" / "cache",
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

    @property
    def done(self):
        return self.try_reading_done()

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
        print("scheduler: start running")
        while len(self.waitlist) > 0 or len(self.running) > 0:
            self.print_status()
            time.sleep(1)  # busy wait...
            self.cleanup_running()

            # check overall memory usage
            if psutil.virtual_memory().free / 1024**3 <= self.total_mem_start_limit:
                continue

            self.start_dockings()

    def start_dockings(self):
        print("start", self.capacity - len(self.running), "jobs")
        for _ in range(self.capacity - len(self.running)):
            if len(self.waitlist) == 0:
                return
            task = self.waitlist.pop()
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
                self.log_fail(job, "death")
                job.suicide()
                continue

            # check for timeout and out-of-memory
            if job.memory > self.proc_mem_limit:
                self.log_fail(job, "memory")
                job.suicide()
                continue
            if job.runtime > self.timeout:
                self.log_fail(job, "timeout")
                job.suicide()
                continue
            still_running.append(i)
        self.running = [job for i, job in enumerate(self.running) if i in still_running]


class TemplateData:
    def __init__(
        self,
        kinodata_path="data/activities-chembl33.csv",
        similar_structures_path="data/templates.csv.gz",
    ):
        self.kinodata_path = kinodata_path
        self.similar_structures_path = similar_structures_path

    @cached_property
    def kinodata(self):
        # activities.activity_id,assays.chembl_id,target_dictionary.chembl_id,molecule_dictionary.chembl_id,molecule_dictionary.max_phase,activities.standard_type,activities.standard_value,activities.standard_units,compound_structures.canonical_smiles,compound_structures.standard_inchi,component_sequences.sequence,assays.confidence_score,docs.chembl_id,docs.year,docs.authors,UniprotID
        return pd.read_csv(self.kinodata_path, index_col="activities.activity_id")

    @cached_property
    def similar_structures(self):
        # activities.activity_id,similar.klifs_structure_id,similar.fp_similarity
        return pd.read_csv(
            self.similar_structures_path, index_col="activities.activity_id"
        )

    @cached_property
    def data(self):
        # activities.activity_id,similar.klifs_structure_id,similar.fp_similarity
        return self.kinodata.join(
            self.similar_structures[
                ~self.similar_structures["similar.klifs_structure_id"] < 0
            ],
            how="inner",
        )


def template_path(structure_id) -> Path:
    TEMPLATE_PATH.mkdir(exist_ok=True)
    return TEMPLATE_PATH / f"{structure_id}.pdb"


def download_template(structure_id) -> Path:
    filename = template_path(structure_id)
    if not filename.exists():
        resp = req.get(
            "https://klifs.net/api_v2/structure_get_pdb_complex",
            {"structure_ID": structure_id},
        )
        with open(filename, "w") as f:
            f.write(resp.text)
    return filename


def download_templates(data: TemplateData):
    print("Download docking templates")
    for structure_id in tqdm.tqdm(data.data["similar.klifs_structure_id"]):
        download_template(structure_id)


def prepare_tasks(data, ident_range=(0, 1e10)) -> List[DockingTask]:
    tasks = list()
    for ident, row in tqdm.tqdm(data.data.iterrows(), total=len(data.data)):
        task = DockingTask(
            ident,
            template_path(row["similar.klifs_structure_id"]),
            data.kinodata.loc[ident, "compound_structures.canonical_smiles"],
        )
        print(ident)
        tasks.append(task)
    print("task preparation done")
    return tasks


if __name__ == "__main__":
    random.seed(int(time.time()))
    data = TemplateData(
        kinodata_path=HERE / "data" / "todo.csv",
        similar_structures_path=HERE / "data" / "templates.csv",
    )

    download_templates(data)

    tasks = prepare_tasks(data)

    print("task prep done")
    output_dir = HERE / "data" / "poses"
    output_dir.mkdir(exist_ok=True, parents=True)
    print("init scheduler")
    scheduler = Scheduler(
        tasks,
        DockingJob,
        capacity=64,
        output_dir=output_dir,
    )

    print("-> start docking")
    scheduler.run()
