"""Cuts run related recipes"""

import os, glob
from cutslib.environ import CUTS_DIR, CUTS_PYENV

######################
# interactive action #
######################

def submit(cutparam):
    from moby2.util import MobyDict

    # load base parameters
    par = MobyDict.from_file(cutparam)

    # runtime parameters
    tpn = par.get("tpn", 20)  # task per node
    nnode = par.get("nnode", 1)
    jobname = par.get("jobname", "get_cuts")
    nproc = tpn * nnode

    # output parameters
    basedir = os.path.dirname(os.path.abspath(cutparam))
    outdir = os.path.join(basedir, par["outdir"])
    runtime = par.get("runtime")
    qos = par.get("qos")
    partition = par.get("partition", "serial")

    # find list of tods to process
    if not(os.path.isdir(outdir)): os.makedirs(outdir)

    # submit one slurm job on each node, so loop over node here
    for n in range(nnode):
        f = open( '%s/submitjob.sh.%d' % (outdir, n), 'w' )
        f.write( '#!/bin/sh\n' )
        f.write( '#SBATCH -N 1\n')
        f.write( '#SBATCH --ntasks-per-node=40\n')
        f.write( '#SBATCH -J %s%d\n' % (jobname,n))
        f.write( '#SBATCH -t %s\n' % runtime )
        f.write( '#SBATCH --qos %s\n' % qos )
        f.write( '#SBATCH --partition %s\n' % partition)
        f.write( '#SBATCH --output=%s/slurmjob.log.%d\n' % (outdir, n))
        f.write( 'module load gcc\n' )
        f.write( 'module load openmpi\n' )
        f.write( 'source activate %s\n' % CUTS_PYENV )  #FIXME
        f.write( 'cd %s\n' % basedir)
        start, end = n*tpn, (n+1)*tpn
        for i in range(start, end):
            f.write('OMP_NUM_THREADS=20 run_cuts -n %d --index %d -f %s & sleep 1\n' % (nproc, i, cutparam))
        f.write('wait\n')
        f.close()
        os.system("sbatch %s/submitjob.sh.%d\n" % (outdir, n))

def list(all=None):
    """List all the latest cuts parameters and the status of it"""
    # get cuts tags information
    for p in glob.glob(CUTS_DIR+"/*"):
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
