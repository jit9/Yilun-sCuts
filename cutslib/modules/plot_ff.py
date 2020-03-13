"""This module plots the flatfield on an array

Options:
    pmin, pmax: plot ranges
    flatfield: if flatfield is other than the output one from param
    freq: if freq is other than the given
    tag: if tag is other than the default tag from param
"""
import moby2.util
import numpy as np
import os.path as op
from cutslib import visual as v

class Module:
    def __init__(self, config):
        self.pmin = config.getfloat("pmin", None)
        self.pmax = config.getfloat("pmax", None)
        self.flatfield = config.get("flatfield", None)
        self.tag = config.get("tag", None)
        self.freq = config.getfloat("freq", None)

    def run(self, p):
        pmin = self.pmin
        pmax = self.pmax
        flatfield = self.flatfield
        tag = self.tag
        freq = self.freq

        ####################
        # output flatfield #
        ####################

        freq = p.i.freq
        if not flatfield:
            ff_name = p.i.ff
        else:
            ff_name = flatfield
        print("Plotting: %s" % ff_name)
        ff = moby2.util.MobyDict.from_file(op.join(p.o.root, ff_name))

        det_uid = np.asarray(ff['det_uid'], dtype=int)
        cal = np.asarray(ff['cal'], dtype=float)
        stable = np.asarray(ff['stable'], dtype=int)

        if pmin is None:
            pmin = cal.min()
        if pmax is None:
            pmax = cal.max()

        # save flatfield plot
        if not tag:
            tag = p.tag
        outfile = p.o.ff + "/ff_%s_cal_output.png" % tag
        print("Saving plot: %s" % outfile)
        v.array_plots(cal, det_uid, season=p.i.season, array=p.i.ar, fr=freq,
                      pmin=pmin, pmax=pmax, title='Flatfield (output) %s' % tag,
                      display='save', save_name=outfile)

        # save stable detector plot
        outfile = p.o.ff + "/ff_%s_stable_output.png" % tag
        print("Saving plot: %s" % outfile)
        v.array_plots(stable, det_uid, title = 'Stable %s' % tag,
                      season=p.i.season, array=p.i.ar, display='save', fr=freq,
                      save_name=outfile)

        ###################
        # input flatfield #
        ###################

        cutParam = moby2.util.MobyDict.from_file(p.i.cutParam)
        ff_name = cutParam.get_deep(('pathologyParams','calibration','flatfield'))
        ff = moby2.util.MobyDict.from_file(ff_name)

        det_uid = np.asarray(ff['det_uid'], dtype=int)
        cal = np.asarray(ff['cal'], dtype=float)
        stable = np.asarray(ff['stable'], dtype=int)

        if pmin is None:
            pmin = cal.min()
        if pmax is None:
            pmax = cal.max()

        # save flatfield plot
        if not tag:
            tag = p.tag
        outfile = p.o.ff + "/ff_%s_cal_input.png" % tag
        print("Saving plot: %s" % outfile)
        v.array_plots(cal, det_uid, season=p.i.season, array=p.i.ar, fr=freq,
                      pmin=pmin, pmax=pmax, title='Flatfield (input) %s' % tag,
                      display='save', save_name=outfile)

        # save stable detector plot
        outfile = p.o.ff + "/ff_%s_stable_input.png" % tag
        print("Saving plot: %s" % outfile)
        v.array_plots(stable, det_uid, title = 'Stable %s' % tag,
                      season=p.i.season, array=p.i.ar, display='save', fr=freq,
                      save_name=outfile)
