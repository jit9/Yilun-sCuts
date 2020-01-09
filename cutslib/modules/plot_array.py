"""In this script I try to produce an array plot of certain pathological
quantity. For example, it would be interesting to see the gain on an
array plot, or the correlation on the array plot"""

import moby2
import pickle, os.path as op
import numpy as np
import numpy.ma as ma
from moby2.analysis.tod_ana.visual import array_plots
from matplotlib import pyplot as plt


def init(config):
    global calibrate, targets, shared_depot
    calibrate = config.getboolean("calibrate", False)
    targets = config.get("targets", None)
    shared_depot = config.get("shared_depot", None)

def run(p):
    global calibrate, targets, shared_depot
    freq = p.i.freq
    array = p.i.ar
    season = p.i.season
    pickle_file = p.i.pickle_file
    ad = moby2.tod.ArrayData.from_fits_table(op.join(shared_depot, 'ArrayData/{}/{}/default.fits'.format(season, array)))
    dets = ad['det_uid'][ad['nom_freq']==freq]

    with open(pickle_file, "rb") as f:
        data = pickle.load(f)

    # target to plot
    if calibrate:
        data['rmsLive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['normLive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['MFELive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['DELive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['jumpLive'] *= data['resp'] * data['ff'][:,np.newaxis]

    # if targets are not specified, all are calculated
    if not targets:
        targets = ['corrLive', 'rmsLive', 'kurtLive', 'skewLive',
                   'normLive', 'darkRatioLive', 'MFELive',
                   'gainLive', 'DELive', 'jumpLive']

    # Run estimators on each target pathological param
    for target in targets:
        target_values = ma.array(data[target], mask=np.logical_not(data['sel']))

        # estimators
        target_mean = np.mean(target_values, axis=1)
        target_median = np.median(target_values, axis=1)
        target_std = np.std(target_values, axis=1, ddof=1)
        # more can be added
        # their output paths are managed by the proj script

        # mean
        outfile = p.o.patho.array.mean + "/" + target + "_mean.png"
        print("Saving plot: %s" % outfile)
        array_plots(target_mean[dets], dets, array=array,
                    season=season, display='save', save_name=outfile,
                    title=target+"_mean", pmin=1e-10)

        # median
        outfile = p.o.patho.array.median + "/" + target + "_median.png"
        print("Saving plot: %s" % outfile)
        array_plots(target_median[dets], dets, array=array,
                    season=season, display='save', save_name=outfile,
                    title=target+"_median", pmin=1e-10)

        # std
        outfile = p.o.patho.array.std + "/" + target + "_std.png"
        print("Saving plot: %s" % outfile)
        array_plots(target_std[dets], dets, array=array,
                    season=season, display='save', save_name=outfile,
                    title=target+"_std", pmin=1e-10)
