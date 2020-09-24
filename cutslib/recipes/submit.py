"""Collection of slurm submission recipes"""

def update_crit(cutparam):
    submit_command("update_crit", cutparam)

def run_tiger(cutparam):
    submit_command("run_cuts", cutparam, cluster='tiger')

def run_della(cutparam):
    submit_command("run_cuts", cutparam, cluster='della')

def run_tiger_cmd(cutparam, cmd):
    submit_command(cmd, cutparam, cluster='tiger')

def run_della_cmd(cutparam, cmd):
    submit_command(cmd, cutparam, cluster='della')

####################
# utility function #
####################

def submit_command(command, cutparam, jobname=None, cluster='tiger'):
    """submit slurm jobs for given binary script based on todloop"""
    import os
    from moby2.util import MobyDict
    CUTS_PYENV = os.environ.get("CUTS_PYENV","myenv")

    # load base parameters
    par = MobyDict.from_file(cutparam)

    # runtime parameters
    tpn = par.get("tpn", 20)  # task per node
    nnode = par.get("nnode", 1)
    if not jobname:
        jobname = par.get("jobname", command)
    nproc = tpn * nnode

    # output parameters
    basedir = os.path.dirname(os.path.abspath(cutparam))
    outdir = os.path.join(basedir, par["outdir"])
    runtime = par.get("runtime")
    qos = par.get("qos", None)
    nomp = par.get("nomp", 1)
    partition = par.get("partition", None)
    # total task per node including nomp
    # note there is a maximum of 40 logical cores in della
    if cluster == 'tiger':   ttpn = 40
    elif cluster == 'della': ttpn = 32
    else:                    ttpn = 32

    # find list of tods to process
    if not(os.path.isdir(outdir)): os.makedirs(outdir)

    # submit one slurm job on each node, so loop over node here
    for n in range(nnode):
        f = open( '%s/submitjob.sh.%d' % (outdir, n), 'w')
        f.write('#!/bin/sh\n' )
        f.write('#SBATCH -N 1\n')
        f.write('#SBATCH --ntasks-per-node=%d\n' % (ttpn))
        f.write('#SBATCH -J %s%d\n' % (jobname,n))
        f.write('#SBATCH -t %s\n' % runtime)
        if qos:     f.write('#SBATCH --qos %s\n' % qos)
        if partition: f.write('#SBATCH -p %s\n' % partition)
        f.write('#SBATCH --output=%s/slurmjob.log.%d\n' % (outdir, n))
        f.write('module load gcc\n' )
        f.write('module load openmpi\n' )
        f.write('source activate %s\n' % CUTS_PYENV )
        f.write('cd %s\n' % basedir)
        start, end = n*tpn, (n+1)*tpn
        for i in range(start, end):
            f.write('OMP_NUM_THREADS=%d %s -n %d --index %d -f %s & sleep 1\n' % (nomp, command, nproc, i, cutparam))
        f.write('wait\n')
        f.close()
        os.system("sbatch %s/submitjob.sh.%d\n" % (outdir, n))
