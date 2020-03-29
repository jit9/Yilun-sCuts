"""This is not a module but a number of utility libraries
for command type module"""

import random, string, os

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
    def run_func(p):
        with open(tfname, "w") as f:
            f.write("#!/bin/bash\n")
            f.write("source ~/.bashrc\n")
            if conda_env:
                f.write(f"source activate {conda_env}\n")
            if module_load:
                for mod in module_loads:
                    f.write(f"module load {mod}\n")
            f.write(cmd)
        os.system(f'bash {tfname}')
        os.system(f"rm {tfname}")
    return run_func
