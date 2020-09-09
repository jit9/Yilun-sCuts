"""This script aims to implement a simpler version of the
routines in the standard cuts pipeline

"""
import numpy as np
import moby2
from moby2.util import Depot
from moby2.scripting import products
from moby2.tod.cuts import CutsVector
from cutslib import MobyDict, preselect, analysis as ana, visual as v
from cutslib.todloop import Routine
from cutslib.analysis import CutsManager, PathologyManager
from cutslib.pathologies import get_detector_params
from cutslib.tools import nextregular
from cutslib.util import dets2sel


class PathologySimple(Routine):
    def __init__(self, depot, det_param, cal_param, jmp_param, lf_param, hf_param, gl_param, out_param):
        self.depot = depot
        self.out_param = out_param
        self.det_param = det_param
        self.cal_param = cal_param
        self.jmp_param = jmp_param
        self.lf_param  = lf_param
        self.hf_param  = hf_param
        self.gl_param  = gl_param
    def initialize(self):
        # load detectors: only needs to run once
        live = MobyDict.from_file(self.det_param['live'])
        excl = MobyDict.from_file(self.det_param['exclude'])
        dark = MobyDict.from_file(self.det_param['dark'])
        self.live = live['det_uid']
        self.excl = excl['det_uid']
        self.dark = dark['det_uid']
        self.depot = Depot(self.depot)
    def execute(self, store):
        tod = store.get('tod')
        # simple preprocessing and downsampling
        moby2.tod.remove_mean(tod)
        moby2.tod.remove_filter_gain(tod)
        tod = tod.copy(resample=2, resample_offset=1)
        # initialize manager for cuts and pathologies
        cman = CutsManager.for_tod(tod)
        pman = PathologyManager.for_tod(tod)
        # find mce cuts and glitch cuts
        cuts_partial = moby2.tod.get_glitch_cuts(tod=tod, params=self.gl_param)
        mce_cuts = moby2.tod.get_mce_cuts(tod)
        # store in manager
        cman.add('mce', mce_cuts)
        cman.add('glitch', cuts_partial.copy())
        # fill glitches before proceeding
        cuts_partial.merge_tod_cuts(mce_cuts)
        moby2.tod.fill_cuts(tod, cuts_partial, extrapolate=False, no_noise=True)
        # find detectors that are all zeros
        zero_sel = tod.data[:,::100].any(axis=1)
        cman.add('zero_sel', zero_sel)
        # remove detector that fluctuates over the full
        # dynamic range
        range_sel = np.std(tod.data,axis=1) < 1e8
        cman.add('range_sel', range_sel)
        # analyze scan
        scan_params = ana.analyze_scan(tod)
        for f_ in ['scan_flags', 'turn_flags']:
            cman.add(f_, CutsVector.from_mask(scan_params[f_]))
        # get ff
        ff_ = moby2.detectors.RelCal.from_dict(self.cal_param['flatfield'])
        ff_sel, ff = ff_.get_property('cal', tod.det_uid, default=0)
        _, stable = ff_.get_property('stable', tod.det_uid, default=False)
        # get biasstep
        resp = products.get_calibration(self.cal_param['config'], tod.info).cal
        resp_sel = resp != 0
        # save various cuts
        cman.add('resp_sel', resp_sel)
        cman.add('ff_sel', ff_sel)
        cman.add('stable', stable)
        # also save the pathologies
        pman.add('resp', resp)
        pman.add('ff', ff)
        # calibrate tod
        cal = resp * ff
        moby2.tod.apply_calibration(tod.data, self.live, cal[self.live])
        # find jump
        jmp = moby2.libactpol.find_jumps(tod.data, self.jmp_param['dsStep'], self.jmp_param['window'])
        pman.add('jump', jmp)
        # fft transform
        trend = moby2.tod.detrend_tod(tod)
        fsw = v.freqSpaceWaterfall(tod, fmin=0.01, fmax=100)
        fdata = fsw.mat**0.5  # fsw given in amplitude^2
        freqs = fsw.matfreqs
        # low-freq analysis
        fmask = (freqs > self.lf_param['fmin']) * (freqs < self.lf_param['fmax'])
        dets4presel = dets2sel(self.live, len(tod.det_uid)) * (cal != 0)
        preselector = preselect.by_median(min_corr=self.lf_param['min_corr'], dets=dets4presel)
        cm = ana.analyze_common_mode(fsw.mat[:,fmask], preselector=preselector, pman=pman)
        # high-freq analysis
        fmask = (freqs > self.hf_param['fmin']) * (freqs < self.hf_param['fmax'])
        nm = ana.analyze_detector_noise(fsw.mat[:,fmask], pman=pman,
                                        n_deproject=self.hf_param['n_deproject'])
        # narrow pathologies to the live candidates only
        pman.restrict_dets(self.live)
        # save pathologies
        # bind useful info such as cuts manager, det lists
        pman.cman = cman
        pman.live = self.live
        pman.excl = self.excl
        pman.dark = self.dark
        # drop the correlation matrices to save disk space
        pman.drop('cc').drop('c')
        self.depot.write_object(pman, tod=tod, **self.out_param)
