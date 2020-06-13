"""This script is aim to find excess power in each detector in a
given list of TODs. What i meant by excess power is the residue
after removing the common modes. By default i will look at
the low frequencies from 0.01Hz to 1Hz and remove one common mode
only. The script can generate the excess power in terms of npy
data files as well as plots if specificied

Example config
--------------

[excess_power]
mpi = True
tod_list = tod.txt
outdir = plots

"""

import os, os.path as op, numpy as np
import matplotlib.pyplot as plt
from cutslib import visual as v
from cutslib import analysis as ana
import moby2
from moby2.util.database import TODList
from cutslib.load import quick_transform


class Module:
    def __init__(self, config):
        self.tod_list = config.get('tod_list')
        self.fmin = config.getfloat('fmin', 0.01)
        self.fmax = config.getfloat('fmax', 1)
        self.outdir = config.get('outdir')
        self.n_deproj = config.getint('n_deproj', 1)
        self.glitchp = {'nSig': 10., 'tGlitch' : 0.007, 'minSeparation': 30, \
                        'maxGlitch': 50000, 'highPassFc': 6.0, 'buffer': 50 }
    def run(self, p):
        tod_list = self.tod_list
        fmin = self.fmin
        fmax = self.fmax
        n_deproj = self.n_deproj
        outdir = self.outdir
        if p.rank == 0:
            if not op.exists(outdir):
                os.makedirs(outdir)
        p.comm.Barrier()
        # load tod
        todnames = TODList.from_file(tod_list)
        for i in range(p.rank,len(todnames),p.size):
            tn = todnames[i]
            print(f"{p.rank}: {i+1}/{len(todnames)} {tn}")
            # create directory
            tod_dir = op.join(outdir, op.basename(tn))
            if not op.exists(tod_dir): os.makedirs(tod_dir)
            # load tod and preprocess tod
            tod = moby2.scripting.get_tod({'filename':tn, 'repair_pointing':True})
            quick_transform(tod, steps=['ff_mce','ff_glitch', 'get_iv', 'demean', 'cal', 'detrend'],
                            glitchp=self.glitchp)
            # focus on detectors with valid IV measurements only
            sel = tod.cal.cal != 0
            # get freq-waterfall
            fsw = v.freqSpaceWaterfall(tod, fmin=fmin, fmax=fmax)
            # deproject common modes
            fdata_excess, _ = ana.deproject_modes(fsw.mat, n_modes=n_deproj, preselector=sel)
            # find frequency that contained the largest power excess
            peak_freq = fsw.matfreqs[np.argmax(fdata_excess, axis=1)]
            # save output
            outfile = op.join(tod_dir, 'peaks.npy')
            np.save(outfile, peak_freq)
            outfile = op.join(tod_dir, 'excess.npy')
            np.save(outfile, fdata_excess)
            outfile = op.join(tod_dir, 'freq.npy')
            np.save(outfile, fsw.matfreqs)
        p.comm.Barrier()
