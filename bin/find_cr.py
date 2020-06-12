"""Compile cosmic-ray like events from a given list of TODs.

Example configuration:
----------------------
depot = '/projects/ACT/yilung/depot'
source_scans = "/projects/ACT/yilung/depot/TODLists/2019_ar6_scan.txt"
outdir = "/scratch/gpfs/yilung/depot/s19_cr/"
cuts_tag = "pa6_f090_s19_c11_v0"
pcuts_tag = "pa6_f090_s19_c11_v0_partial"
cal_tag = "pa6_f090_s19_c11_v0"
template = "cr_template.npy"

"""
import os.path as op, numpy as np
import moby2, pickle
import argparse

from cutslib.todloop import TODLoop, Routine
from cutslib.routines.tod import LoadTOD
from cutslib.routines.glitch import FindCR
from cutslib import mpi, MobyDict
import cutslib.glitch as gl

# define command line options
parser = argparse.ArgumentParser()
parser.add_argument("param", help="cutparam")
args = parser.parse_args()

# get parameters
cutparam = MobyDict.from_file(args.param)

# initialize loop
loop = TODLoop()

# add list of tod
source_scans = cutparam.get('source_scans')
loop.add_tod_list(source_scans)
outdir = cutparam.get('outdir')

# load tod
loop.add_routine(LoadTOD())
params = {
    'outdir': outdir,
    'cuts_tag': cutparam.get('cuts_tag'),
    'pcuts_tag': cutparam.get('pcuts_tag'),
    'cal_tag': cutparam.get('cal_tag'),
    'template': cutparam.get('template'),
    'force_partial': cutparam.get('force_partial'),
    'glitchp': cutparam.get('glitchp')
}
loop.add_routine(FindCR(**params))
comm = mpi.COMM_WORLD
size = comm.size
rank = comm.rank
if rank == 0:
    import os
    os.makedirs(outdir)
comm.Barrier()
loop.run_fparallel(n_workers=size, rank=rank)
comm.Barrier()
