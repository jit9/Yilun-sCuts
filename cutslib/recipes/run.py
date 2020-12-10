"""Cuts run related recipes"""

import os, glob
from cutslib.environ import CUTS_DIR, CUTS_PYENV
import numpy as np

######################
# interactive action #
######################

def list(all=None):
    """List all the latest cuts parameters and the status of it"""
    # get cuts tags information
    for p in glob.glob(CUTS_DIR+"/pa*"):
        # have an option to ignore a directory that I'm not working on
        # by adding a .cutsignore in the directory
        if os.path.exists(os.path.join(p, '.cutsignore')): continue
        tag = os.path.basename(p)
        # get all cutparam
        cpars = [f.strip() for f in glob.glob(p+"/cutparam*.par")]
        if all != "all":
            # get latest version
            ver_latest = max([int(c.split(".par")[0][-1]) for c in cpars])
            # get latest cutparam
            cpar = "cutparams_v{}".format(ver_latest)
            # get absolute path
            cpar_path = os.path.join(p,cpar+".par")
            # get summary of slurm job
            slm_summary = get_slurm_summary(cpar_path)
            # get slurm progress
            n_tod = int(get_total_tod(cpar_path))
            n_tod_done = int(get_done_tod(cpar_path))
            progress = f"{n_tod_done}/{n_tod} {n_tod_done/n_tod*100:.1f}%"
            # print out result
            print(f"{tag:>10} {cpar:>15} {progress:>10} {slm_summary:>25} {cpar_path}")
        else:
            for cpar_path in cpars:
                cpar = os.path.basename(cpar_path)
                # get summary of slurm job
                slm_summary = get_slurm_summary(cpar_path)
                # get slurm progress
                n_tod = int(get_total_tod(cpar_path))
                n_tod_done = int(get_done_tod(cpar_path))
                progress = f"{n_tod_done}/{n_tod} {n_tod_done/n_tod*100:.1f}%"
                # print out result
                print(f"{tag:>10} {cpar:>15} {progress:>10} {slm_summary:>25} {cpar_path}")

def errors(cpar_path):
    ver = cpar_path.split(".par")[0][-1]  # FIXME
    basedir = os.path.dirname(os.path.abspath(cpar_path))
    run_dir = os.path.join(basedir, f"run_v{ver}")
    os.system(f"cat {run_dir}/slurmjob.log.* | grep '.*Error:'")

def logs(cpar_path):
    ver = cpar_path.split(".par")[0][-1]  # FIXME
    basedir = os.path.dirname(os.path.abspath(cpar_path))
    run_dir = os.path.join(basedir, f"run_v{ver}")
    os.system(f"tail -f {run_dir}/slurmjob.log.*")

def errors_tod(cpar_path, match):
    ver = cpar_path.split(".par")[0][-1]  # FIXME
    basedir = os.path.dirname(os.path.abspath(cpar_path))
    run_dir = os.path.join(basedir, f"run_v{ver}")
    os.system("cat %s/error_list.txt | grep -i %s | awk '{print $2}' | sort | uniq" % (run_dir, match))

def errors_stats(cpar_path):
    ver = cpar_path.split(".par")[0][-1]  # FIXME
    basedir = os.path.dirname(os.path.abspath(cpar_path))
    run_dir = os.path.join(basedir, f"run_v{ver}")
    ntod_total = int(get_total_tod(cpar_path))
    # load error_list.txt
    with open(f"{run_dir}/error_list.txt", "r") as f:
        lines = f.readlines()
    lines = [l.strip() for l in lines]
    # todnames = [l.split(' ')[1] for l in lines]
    # get everything after todname and with whitespace striped
    texts = [' '.join(l.split(' ')[2:]).strip() for l in lines]
    errors = np.unique(texts)
    print(" num| rel%| abs%|text")
    for e in errors:
        ntod = texts.count(e)
        pctg = ntod / len(texts)
        pctg2 = ntod / ntod_total
        print(f"{texts.count(e):4d}|{pctg*100:4.1f}%|{pctg2*100:4.1f}%|{e}")

####################
# utility function #
####################

def get_job_name(cpar_path):
    """Get the slurm job name of a cutparam"""
    return os.popen("cat %s | grep jobname | cut -d'\"' -f2 " % cpar_path).read().strip()

def get_slurm_jobs(cpar_path):
    """Get all slurm jobs info for a given cutparam"""
    from cutslib.environ import CUTS_USER
    # from moby2.util import MobyDict  # this is slower
    # jobname = MobyDict.from_file(cpar_path).get("jobname")
    jobname = get_job_name(cpar_path)
    lines = os.popen("squeue -u %s " % CUTS_USER).read().strip().split('\n')
    return [l.split() for l in lines[1:] if jobname in l]

def get_slurm_summary(cpar_path):
    """Get a text summary of the slurm job status"""
    jobs = get_slurm_jobs(cpar_path)
    tot = len(jobs)
    n_run = 0
    n_pd = 0
    run_nodes = 0
    pd_nodes = 0
    for j in jobs:
        if j[4] == "R":
            n_run += 1
            run_nodes += int(j[6])
        elif j[4] == "PD":
            n_pd += 1
            pd_nodes += int(j[6])
    return "R:j={};n={}|PD:j={};n={}".format(n_run,run_nodes,\
                                             n_pd,pd_nodes)

def get_total_tod(cpar_path):
    """Get the total number of TODs to be done"""
    return os.popen("cat %s | grep '^source_scans.*' | cut -d'\"' -f2 | xargs cat | wc -l" % cpar_path).read().strip()

def get_done_tod(cpar_path):
    """Count the number of TODs completed"""
    # get run_dir
    ver = cpar_path.split('.par')[0][-1]
    run_dir = os.path.join(os.path.dirname(cpar_path),"run_v%s"%ver)
    if os.path.exists(run_dir) and len(glob.glob("%s/*.db*"%run_dir)) > 0:
        # note that this assumes that the data entry starts with 1
        return os.popen("cat %s/*.db* | grep '^1' | sort | uniq | wc -l" % run_dir).read().strip()
    else:
        return 0
