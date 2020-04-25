"""This script is aim to generate various waterfall plots and correlation
matrix for a list of tods.

Example config
--------------

[waterfall]
mpi = True
type = external
file = waterfall.py
tod_list = tod_ar7.txt
fmin = 0.01
cov_fmin = 10
cov_fmax = 20
n_deproj = 3
outdir = plots/ar7/


"""
import os.path as op, numpy as np
import matplotlib.pyplot as plt
from cutslib import visual as v
from cutslib import analysis as ana
import moby2
from moby2.util.database import TODList

def cov_frange(fsw, sel, fmin, fmax, n_deproj=0, plot=True, vmin=-1, vmax=1):
    freq = fsw.matfreqs
    fmask = (freq > fmin) * (freq < fmax)
    fmodes = fsw.mat[np.ix_(sel, fmask)]
    fmodes = ana.deproject_modes(fmodes, n_modes=n_deproj)
    cov = ana.corrmat(fmodes)
    if plot:
        plt.figure(figsize=(10.5,10.5))
        plt.imshow(cov, cmap='jet', origin='lower', vmin=vmin, vmax=vmax)
        plt.colorbar(shrink=0.8)
        plt.xlabel('dets')
        plt.ylabel('dets')
        plt.title(f'correlation [{fmin}Hz, {fmax}Hz]')
    return cov


class Module:
    def __init__(self, config):
        self.tod_list = config.get('tod_list')
        self.fmin = config.getfloat('fmin')
        self.cov_fmin = config.getfloat('cov_fmin', 10)
        self.cov_fmax = config.getfloat('cov_fmax', 10)
        self.outdir = config.get('outdir')
        self.n_deproj = config.getint('n_deproj')

    def run(self, p):
        tod_list = self.tod_list
        fmin = self.fmin
        cov_fmin = self.cov_fmin
        cov_fmax = self.cov_fmax
        n_deproj = self.deproj
        outdir = self.outdir
        # load tod
        todnames = TODList.from_file(tod_list)
        for tn in todnames[p.rank:len(todnames):p.size]:
            print(f"{p.rank}: {tn}")
            tod = moby2.scripting.get_tod({'filename':tn, 'repair_pointing':True})
            # create freq-waterfall object
            fsw = v.freqSpaceWaterfall(tod, fmin=fmin)
            # plot only tes detectors
            sel = tod.info.array_data['det_type'] == 'tes'
            # create waterfall plot and save it
            outfile = op.join(outdir, op.basename(tn)+'_fsw.png')
            fsw.plot(selection=sel, vmin=2, filename=outfile, show=False)
            # create time-waterfall object
            tsw = v.timeSpaceWaterfall(tod)
            outfile = op.join(outdir, op.basename(tn)+'_tsw.png')
            tsw.plot(selection=sel, title=f'Time-domain waterfall:{op.basename(tn)}', filename=outfile)
            # create correlation plot
            cov_frange(fsw, sel, cov_fmin, cov_fmax, n_deproj=deproj);
            outfile = op.join(outdir, op.basename(tn)+'_cov.png')
            plt.savefig(outfile)
            plt.close()
        p.comm.Barrier()
