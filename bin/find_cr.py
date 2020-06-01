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
    'template': cutparam.get('template')
}
loop.add_routine(FindCR(**params))
comm = mpi.COMM_WORLD
size = comm.size
rank = comm.rank
loop.run_fparallel(n_workers=size, rank=rank)

comm.Barrier()
