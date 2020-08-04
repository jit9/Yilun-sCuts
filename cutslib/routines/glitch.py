
import moby2, pickle
import os.path as op, numpy as np
from cutslib import glitch as gl
from cutslib import Depot
from cutslib.load import quick_transform
from cutslib.todloop import Routine


class FindCR(Routine):
    def __init__(self, **params):
        """A routine to perform a template search in the glitches. It
        extract the snippets of each glitch that match the given
        template above certain snr. The results are written as a pickle.

        Parameters
        ----------
        cuts_tag (str)     : tag for det cuts
        pcuts_tag (str)    : tag for partial cuts
        cal_tag (str)      : tag for calibration
        mce_ndet_lim (int) : filter for mce glitch with ndet above this limit
        mce_len_lim (int)  : filter for mce glitch with len below this limit
        snr_lim (int)      : signal-to-noise threshold for template matching
        template (str)     : path to a template to match (npy file)
        outdir (str)       : output dir, needs to exist
        force (bool)       : when True existing file will be overridden
        glitchp (dict)     : glitch parameters
        steps (list)       : list of pre-processing steps
        uncut_only (bool)  : whether to use the uncut detectors only

        """
        Routine.__init__(self)
        self.cuts_tag = params.get('cuts_tag')
        self.pcuts_tag = params.get('pcuts_tag')
        self.cal_tag = params.get('cal_tag')
        self.filter_mce = params.get('filter_mce', False)
        self.mce_ndet_lim = params.get('mce_ndet_lim', 200)
        self.mce_len_lim = params.get('mce_len_lim', 100)
        self.snr_lim = params.get('snr_lim', 20)
        self.template = params.get('template')
        self.outdir = params.get('outdir', 'out')
        self.force = params.get('force', False)
        self.force_partial = params.get('force_partial', False)
        self.glitchp = params.get('glitchp', None)
        self.steps = params.get('steps',['demean','cal'])
        self.uncut_only = params.get('uncut_only', True)
    def initialize(self):
        self.depot = moby2.util.Depot(Depot().root)
        self.template = np.load(self.template)
    def execute(self, store):
        outfile = op.join(self.outdir, self.get_name()+'.pkl')
        # skip if already done
        if op.exists(outfile) and not self.force:
            self.logger.info(f"{outfile} exists, skipping...")
            return
        # check whether we have cuts before running anything
        todname = self.get_name()
        first5 = todname[:5]
        cuts_file = op.join(Depot().root, 'TODCuts', self.cuts_tag, first5, todname+'.cuts')
        cuts_exist = op.exists(cuts_file)
        if not cuts_exist: return
        # load tod
        self.logger.info("Loading TOD")
        tod = moby2.scripting.get_tod({'filename': self.get_name(), 'repair_pointing': True})
        # load cuts
        cuts = self.depot.read_object(moby2.TODCuts, tag=self.cuts_tag, tod=tod)
        tod.cuts = cuts
        # get partial cuts
        # 1. fill mce cuts
        self.logger.info("Filling MCE cuts")
        mce_cuts = moby2.tod.get_mce_cuts(tod)
        moby2.tod.fill_cuts(tod, mce_cuts, no_noise=True)
        # 2. generate partial cuts
        if not self.force_partial:
            # load partial cuts
            self.logger.info("Loading partial cuts")
            pcuts = self.depot.read_object(moby2.TODCuts, tag=self.pcuts_tag, tod=tod)
        else:  # force partial
            self.logger.info("Generating partial cuts")
            pcuts = moby2.tod.get_glitch_cuts(tod=tod, params=self.glitchp)
        tod.partial = pcuts
        if 'cal' in self.steps:
            # load calibration
            self.logger.info("Loading calibration")
            cal = self.depot.read_object(moby2.Calibration, tag=self.cal_tag, tod=tod)
            tod.cal = cal
        # demean and calibrate if needed
        self.logger.info("Transforming tod")
        quick_transform(tod, steps=self.steps)
        self.logger.info("Extracting snippets")
        # count number of dets glitching at each sample
        dets = tod.cuts.get_uncut() if self.uncut_only else None
        count = gl.glitch_det_count(tod.partial, None)
        # get rid of mce cuts in the partial cuts
        # this can be done by filtering out time ranges that more than 200 dets
        # are affected by cuts and also check for a cuts shorter than bufferred
        # duration
        if self.filter_mce:
            excess = count > self.mce_ndet_lim  # 200 by default
            if np.sum(excess) > 0:
                cv = gl.CutsVector.from_mask(excess)
                len_mask = (cv[:,1] - cv[:,0]) < self.mce_len_lim  # 100 by default
                cv = cv[len_mask]
                gl.fill_cv(count, cv)
        # get "events" that correspond to a localized range of time that many dets
        # are affected
        events = gl.CutsVector.from_mask(count>0)
        dets = tod.cuts.get_mask() if self.uncut_only else np.ones_like(tod.det_uid, dtype=bool)
        snippets = gl.affected_snippets_from_cv(tod, tod.partial, events, dets)
        # filter cr by template
        cr_snippets = list(filter(lambda s: s.max_snr_template(self.template)>self.snr_lim, snippets))
        # write file
        with open(outfile, "wb") as f:
            pickle.dump(cr_snippets, f)
        self.logger.info(f"Written: {outfile}")
