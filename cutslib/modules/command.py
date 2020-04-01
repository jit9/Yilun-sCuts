"""This is not a module but a number of utility libraries
for command type module

Options
-------
module_load: list of modules to load
conda_env: if we were to run in a given conda env
use_slurm: whether we want to submit jobs to slurm instead
slurm_nnode: option for sbatch -N
slurm_tpn: option for sbatch --ntasks-per-node
slurm_jobname: option for sbatch -J
slurm_runtime: option for sbatch -t
slurm_qos: option for sbatch --qos
nomp: option for OMP_NUM_THREADS
command: actual command to run

"""

import random, string

def generate_random_filename():
    return '.tmp_' + \
           ''.join(random.choice(string.ascii_letters) for x in range(5))

def create_run_func(config):
    cmd = config.get("command")
    if not cmd:
        raise ValueError("No command given!")
    # temperatory filename
    tfname = generate_random_filename()
    # check for module load
    module_load = config.get("module_load", None)
    if module_load:
        module_loads = module_load.split()
    # check for conda environ
    conda_env = config.get("conda_env", None)
    # check whether to submit job through slurm
    use_slurm = config.getboolean("use_slurm", False)
    if use_slurm:
        nnode = config.getint("slurm_nnode", 1)
        tpn = config.getint("slurm_tpn", 40)
        jobname = config.get("slurm_jobname", "myjob")
        runtime = config.get("slurm_runtime", "1:00:00")
        qos = config.get("slurm_qos", "debug")
        slurm_headers = [
            f'#SBATCH -N {nnode}\n',
            f'#SBATCH --ntasks-per-node={tpn}\n',
            f'#SBATCH -J {jobname}\n',
            f'#SBATCH -t {runtime}\n',
            f'#SBATCH --qos {qos}\n'
        ]
    # check whether to specify openmp
    nomp = config.getint("nomp", 0)
    def run_func(p):
        import os
        with open(tfname, "w") as f:
            f.write("#!/bin/bash\n")
            # write slurm headers if that's what we want
            if use_slurm:
                for h in slurm_headers:
                    f.write(h)
            f.write("source ~/.bashrc\n")
            if conda_env:
                f.write(f"source activate {conda_env}\n")
            if module_load:
                for mod in module_loads:
                    f.write(f"module load {mod}\n")
            # if we want openmp
            if nomp == 0:
                f.write(cmd)
            else:
                f.write(f'OMP_NUM_THREADS={nomp} {cmd}')
        if not use_slurm:
            os.system(f'bash {tfname}')
            os.system(f"rm {tfname}")
        else:
            os.system(f'sbatch {tfname}')
    return run_func
