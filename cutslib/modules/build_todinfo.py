"""This module is a thin wrapper around the build_todinfo.py in tenki
that creates the todinfo.hdf file required by enki.

Pre-requisite
-------------
 - enki
   module load enki

Parameters
-----------
n_tasks: number of mpi tasks
dataset: name of the dataset
sel: sel passed to build_todinfo

"""
import subprocess


class Module:
    def __init__(self, config):
        self.n_tasks = config.getint("n_tasks", 1)
        self.dataset = config.get("dataset")
        self.sel = config.get("sel")

    def run(self, p):
        import os.path as op, os
        n_tasks = self.n_tasks
        dataset = self.dataset
        sel = self.sel
        ds_root = op.join(p.map_root, dataset)
        cmd = f'enki build_todinfo {ds_root}/todinfo.txt "{sel}"'\
              f' {ds_root}/todinfo.hdf --filedb {ds_root}/filedb.txt'\
              f' --dataset {dataset}'
        if n_tasks > 1:
            cmd = f'mpirun -n {n_tasks} ' + cmd
        print(cmd)
        with open(".tmpfile", "w") as f:
            f.write("#!/bin/bash\n")
            f.write("source ~/.bashrc\n")
            f.write("source activate myenv\n")
            f.write("module load enki\n")
            f.write("module load so_stack\n")
            f.write(cmd)
        os.system("bash .tmpfile")
        # proc = subprocess.Popen(['bash', '.tmpfile'], stdout=subprocess.PIPE)
        # retval = proc.wait()
