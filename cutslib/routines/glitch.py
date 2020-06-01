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

        """
        Routine.__init__(self)
        self.cuts_tag = params.get('cuts_tag')
        self.pcuts_tag = params.get('pcuts_tag')
        self.cal_tag = params.get('cal_tag')
        self.mce_ndet_lim = params.get('mce_ndet_lim', 200)
        self.mce_len_lim = params.get('mce_len_lim', 100)
        self.snr_lim = params.get('snr_lim', 20)
        self.template = params.get('template')
        self.outdir = params.get('outdir', 'out')
    def initialize(self):
        self.depot = moby2.util.Depot(Depot().root)
        self.template = np.load(self.template)
    def execute(self, store):
        tod = store.get(self.inputs['tod'])
        # load cuts
        cuts = self.depot.read_object(moby2.TODCuts, tag=self.cuts_tag, tod=tod)
        tod.cuts = cuts
        # load partial cuts
        pcuts = self.depot.read_object(moby2.TODCuts, tag=self.pcuts_tag, tod=tod)
        tod.partial = pcuts
        # load calibration
        cal = self.depot.read_object(moby2.Calibration, tag=self.cal_tag, tod=tod)
        tod.cal = cal
        # demean and calibrate
        quick_transform(tod, steps=['demean','cal'])
        # count number of dets glitching at each sample
        count = gl.glitch_det_count(tod.partial, tod.cuts.get_uncut())
        # get rid of mce cuts in the partial cuts
        # this can be done by filtering out time ranges that more than 200 dets
        # are affected by cuts and also check for a cuts shorter than bufferred
        # duration
        excess = count > self.mce_ndet_lim  # 200 by default
        if np.sum(excess) > 0:
            cv = gl.CutsVector.from_mask(excess)
            len_mask = (cv[:,1] - cv[:,0]) < self.mce_len_lim  # 100 by default
            cv = cv[len_mask]
            gl.fill_cv(count, cv)
        # get "events" that correspond to a localized range of time that many dets
        # are affected
        events = gl.CutsVector.from_mask(count>0)
        snippets = gl.affected_snippets_from_cv(tod, tod.partial, events, tod.cuts.get_mask())
        # filter cr by template
        cr_snippets = list(filter(lambda s: s.max_snr_template(self.template)>self.snr_lim, snippets))
        # write file
        outfile = op.join(self.outdir, self.get_name()+'.pkl')
        with open(outfile, "wb") as f:
            pickle.dump(cr_snippets, f)
        self.logger.info(f"Written: {outfile}")
