"""In this script I try to produce an array plot of certain pathological
quantity. For example, it would be interesting to see the gain on an
array plot, or the correlation on the array plot"""

import moby2
import pickle, os.path as op
import numpy as np
import numpy.ma as ma
from cutslib.visual import array_plots
from cutslib import SharedDepot
from matplotlib import pyplot as plt


class Module:
    def __init__(self, config):
        self.calibrate = config.getboolean("calibrate", False)
        self.targets = config.get("targets", None)
        self.estimator_name = config.get("estimator", "mean")
        self.gain_l = config.getfloat("gain_l", None)
        self.gain_h = config.getfloat("gain_h", None)
        self.corr_l = config.getfloat("corr_l", None)
        self.corr_h = config.getfloat("corr_h", None)
        self.rms_l = config.getfloat("rms_l", None)
        self.rms_h = config.getfloat("rms_h", None)
        self.kurt_l = config.getfloat("kurt_l", None)
        self.kurt_h = config.getfloat("kurt_h", None)
        self.skew_l = config.getfloat("skew_l", None)
        self.skew_h = config.getfloat("skew_h", None)
        self.norm_l = config.getfloat("norm_l", None)
        self.norm_h = config.getfloat("norm_h", None)
        self.de_l = config.getfloat("de_l", None)
        self.de_h = config.getfloat("de_h", None)

    def run(self, p):
        calibrate = self.calibrate
        targets = self.targets
        estimator_name = self.estimator_name
        gain_l = self.gain_l
        gain_h = self.gain_h
        corr_l = self.corr_l
        corr_h = self.corr_h
        rms_l = self.rms_l
        rms_h = self.rms_h
        kurt_l = self.kurt_l
        kurt_h = self.kurt_h
        skew_l = self.skew_l
        skew_h = self.skew_h
        norm_l = self.norm_l
        norm_h = self.norm_h
        de_l = self.de_l
        de_h = self.de_h

        freq = p.i.freq
        array = p.i.ar
        season = p.i.season
        pickle_file = p.i.pickle_file

        # load shared depot
        shared_depot = SharedDepot()
        ad = moby2.tod.ArrayData.from_fits_table(
            op.join(shared_depot.root,
                    'ArrayData/{}/{}/default.fits'.format(season, array)))
        dets = ad['det_uid'][ad['nom_freq']==freq]

        with open(pickle_file, "rb") as f:
            data = pickle.load(f)

        # target to plot
        # if targets are not specified, all are calculated
        if not targets:
            targets = ['corrLive', 'rmsLive', 'kurtLive', 'skewLive',
                       'normLive', 'darkRatioLive', 'MFELive',
                       'gainLive', 'DELive', 'jumpLive']
        else:
            targets = targets.split(',')

        if calibrate:
            data['rmsLive'] *= data['resp'] * data['ff'][:,np.newaxis]
            data['normLive'] *= data['resp'] * data['ff'][:,np.newaxis]
            data['MFELive'] *= data['resp'] * data['ff'][:,np.newaxis]
            data['DELive'] *= data['resp'] * data['ff'][:,np.newaxis]
            data['jumpLive'] *= data['resp'] * data['ff'][:,np.newaxis]

        # Get estimator
        if estimator_name == "mean":
            estimator = lambda x: np.mean(x, axis=1)
        elif estimator_name == "median":
            estimator = lambda x: np.median(x, axis=1)
        elif estimator_name == "std":
            estimator = lambda x: np.std(x, axis=1, ddof=1)
        else:
            raise NotImplementedError("Estimator %s is not implemented!" %
                                      estimator_name)

        # Run estimators on each target pathological param
        target = 'gainLive'
        if target in targets:
            target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
            target_est = estimator(target_values)
            outfile = p.o.patho.array.root + "/" + target + "_%s.png" % estimator_name
            print("Saving plot: %s" % outfile)
            array_plots(target_est[dets], dets, array=array,
                        season=season, display='save', save_name=outfile,
                        title=target+"_%s" % estimator_name, pmin=gain_l, pmax=gain_h)

        target = 'corrLive'
        if target in targets:
            target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
            target_est = estimator(target_values)
            outfile = p.o.patho.array.root + "/" + target + "_%s.png" % estimator_name
            print("Saving plot: %s" % outfile)
            array_plots(target_est[dets], dets, array=array,
                        season=season, display='save', save_name=outfile,
                        title=target+"_%s" % estimator_name, pmin=corr_l, pmax=corr_h)

        target = 'rmsLive'
        if target in targets:
            target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
            target_est = estimator(target_values)
            outfile = p.o.patho.array.root + "/" + target + "_%s.png" % estimator_name
            print("Saving plot: %s" % outfile)
            array_plots(target_est[dets], dets, array=array,
                        season=season, display='save', save_name=outfile,
                        title=target+"_%s" % estimator_name, pmin=rms_l, pmax=rms_h)

        target = 'kurtLive'
        if target in targets:
            target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
            target_est = estimator(target_values)
            outfile = p.o.patho.array.root + "/" + target + "_%s.png" % estimator_name
            print("Saving plot: %s" % outfile)
            array_plots(target_est[dets], dets, array=array,
                        season=season, display='save', save_name=outfile,
                        title=target+"_%s" % estimator_name, pmin=kurt_l, pmax=kurt_h)

        target = 'skewLive'
        if target in targets:
            target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
            target_est = estimator(target_values)
            outfile = p.o.patho.array.root + "/" + target + "_%s.png" % estimator_name
            print("Saving plot: %s" % outfile)
            array_plots(target_est[dets], dets, array=array,
                        season=season, display='save', save_name=outfile,
                        title=target+"_%s" % estimator_name, pmin=skew_l, pmax=skew_h)

        target = 'normLive'
        if target in targets:
            target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
            target_est = estimator(target_values)
            outfile = p.o.patho.array.root + "/" + target + "_%s.png" % estimator_name
            print("Saving plot: %s" % outfile)
            array_plots(target_est[dets], dets, array=array,
                        season=season, display='save', save_name=outfile,
                        title=target+"_%s" % estimator_name, pmin=norm_l, pmax=norm_h)

        target = 'MFELive'
        if target in targets:
            target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
            target_est = estimator(target_values)
            outfile = p.o.patho.array.root + "/" + target + "_%s.png" % estimator_name
            print("Saving plot: %s" % outfile)
            array_plots(target_est[dets], dets, array=array,
                        season=season, display='save', save_name=outfile,
                        title=target+"_%s" % estimator_name, pmin=mfe_l, pmax=mfe_h)

        target = 'DELive'
        if target in targets:
            target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
            target_est = estimator(target_values)
            outfile = p.o.patho.array.root + "/" + target + "_%s.png" % estimator_name
            print("Saving plot: %s" % outfile)
            array_plots(target_est[dets], dets, array=array,
                        season=season, display='save', save_name=outfile,
                        title=target+"_%s" % estimator_name, pmin=de_l, pmax=df_h)
