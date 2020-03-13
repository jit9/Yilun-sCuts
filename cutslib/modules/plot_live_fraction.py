"""This script aims to produce an array plot showing the fraction
of time each detector is alive in the given list. """

import moby2
import pickle, os.path as op, numpy as np
from cutslib.visual import array_plots
from matplotlib import cm

class Module:
    def __init__(self, config):
        pass

    def run(self, p):
        # get afs
        freq = p.i.freq
        array = p.i.ar
        season = p.i.season
        # load pickle file
        print("Loading pickle file: %s" % p.o.pickle_file)
        with open(p.o.pickle_file, "rb") as f:
            data = pickle.load(f)
        # get detector lists
        cutParam = moby2.util.MobyDict.from_file(p.i.cutParam)
        ld_dict = cutParam.get_deep(('pathologyParams','detectorLists','live'))
        exclude_dict = cutParam.get_deep(('pathologyParams','detectorLists','exclude'))
        ld_data = moby2.util.MobyDict.from_file(ld_dict)
        exclude_data = moby2.util.MobyDict.from_file(exclude_dict)
        dets = [d for d in ld_data['det_uid'] if d not in exclude_data['det_uid']]
        # make array plots
        outfile = op.join(p.o.root, "live_frac.png")
        live_frac = np.sum(data['sel'], axis=1)*1.0/data['sel'].shape[-1]

        print("Saving plot: %s" % outfile)
        array_plots(live_frac[dets], dets, array=array, season=season,
                    display='save', save_name=outfile, title="Live fraction: %s" % p.tag,
                    cmap=cm.copper)
