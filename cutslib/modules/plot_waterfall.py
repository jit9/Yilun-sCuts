"""This script is aim to generate various waterfall plots and correlation
matrix for a list of tods.

Example config
--------------
[waterfall]
mpi=True
tod_list = tod.txt
fmin = 0.01
outdir = plots/ar7

plot_fsw = True
plot_tsw = True
plot_tod = True
plot_fft = True

"""
import os, os.path as op, numpy as np
import matplotlib.pyplot as plt
from cutslib import visual as v
from cutslib import analysis as ana
import moby2
from moby2.util.database import TODList
from cutslib.load import quick_transform

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
        self.fmin = config.getfloat('fmin', 0.01)
        self.fmax = config.getfloat('fmax', 100)
        self.cov_fmin = config.getfloat('cov_fmin', 10)
        self.cov_fmax = config.getfloat('cov_fmax', 10)
        self.outdir = config.get('outdir')
        self.n_deproj = config.getint('n_deproj', 0)
        self.glitchp = {'nSig': 10., 'tGlitch' : 0.007, 'minSeparation': 30, \
                        'maxGlitch': 50000, 'highPassFc': 6.0, 'buffer': 50 }
        # decide on what to plot
        self.plot_fsw = config.getboolean('plot_fsw', False)
        self.plot_tsw = config.getboolean('plot_tsw', False)
        self.plot_cov = config.getboolean('plot_cov', False)
        self.plot_tod = config.getboolean('plot_tod', False)
        self.plot_fft = config.getboolean('plot_fft', False)

    def run(self, p):
        tod_list = self.tod_list
        fmin = self.fmin
        fmax = self.fmax
        cov_fmin = self.cov_fmin
        cov_fmax = self.cov_fmax
        n_deproj = self.n_deproj
        outdir = self.outdir
        if p.rank == 0:
            if not op.exists(outdir):
                os.makedirs(outdir)
        p.comm.Barrier()
        # load tod
        todnames = TODList.from_file(tod_list)
        for tn in todnames[p.rank:len(todnames):p.size]:
            print(f"{p.rank}: {tn}")
            # create directory
            tod_dir = op.join(outdir, op.basename(tn))
            if not op.exists(tod_dir): os.makedirs(tod_dir)
            tod = moby2.scripting.get_tod({'filename':tn, 'repair_pointing':True})
            # fill partial cuts
            mce_cuts = moby2.tod.get_mce_cuts(tod)
            moby2.tod.fill_cuts(tod, mce_cuts, no_noise=True)
            pcuts = moby2.tod.get_glitch_cuts(tod=tod, params=self.glitchp)
            moby2.tod.fill_cuts(tod, pcuts, no_noise=True)
            # get detectors to plot
            sel = tod.info.array_data['det_type'] == 'tes'
            # quick transform
            quick_transform(tod, steps=['demean'])
            # remove detectors with zero reading
            nonzero_sel = tod.data[:,::100].any(axis=1)
            sel *= nonzero_sel
            if self.plot_fsw:
                # create freq-waterfall object
                fsw = v.freqSpaceWaterfall(tod, fmin=fmin, fmax=fmax)
                # create waterfall plot and save it
                outfile = op.join(tod_dir, 'fsw_col.png')
                fsw.plot(selection=sel, vmin=2, filename=outfile, show=False)
                outfile = op.join(tod_dir, 'fsw_row.png')
                fsw.plot(selection=sel, vmin=2, filename=outfile, show=False, rowDominance=True)
            if self.plot_tsw:
                # create time-waterfall object
                tsw = v.timeSpaceWaterfall(tod)
                outfile = op.join(tod_dir, 'tsw_col.png')
                tsw.plot(selection=sel, title=f'Time-domain waterfall:{op.basename(tn)} CD', filename=outfile)
                outfile = op.join(tod_dir, 'tsw_row.png')
                tsw.plot(selection=sel, title=f'Time-domain waterfall:{op.basename(tn)} RD', filename=outfile, rowDominance=True)
            # create correlation plot
            if self.plot_cov:
                cov_frange(fsw, sel, cov_fmin, cov_fmax, n_deproj=n_deproj);
                outfile = op.join(tod_dir, 'cov.png')
                plt.savefig(outfile)
                plt.close()
            if self.plot_tod:
                plt.figure(figsize=(8,6))
                v.plot_tod(tod, det=np.where(sel)[0])
                outfile = op.join(tod_dir, 'tod.png')
                plt.savefig(outfile)
                plt.close()
            if self.plot_fft:
                assert self.plot_fsw, "Need to compute freq-space waterfall first!"
                # make fft plots
                plt.figure(figsize=(8,6))
                plt.plot(fsw.matfreqs, fsw.mat[sel,:].T, alpha=0.1, c='k');
                plt.axvline(0.05,ls='--',c='r', alpha=0.1)
                plt.ylim([1e2,1e16])
                plt.ylabel('fft power (DAQ)')
                plt.yscale('log')
                plt.xlabel('freqs (Hz)')
                plt.xscale('log')
                outfile = op.join(tod_dir, 'fft.png')
                plt.savefig(outfile)
                plt.close()
        p.comm.Barrier()
