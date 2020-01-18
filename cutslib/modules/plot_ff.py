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

def init(config):
    global pmin, pmax, flatfield, tag, freq
    pmin = config.getfloat("pmin", None)
    pmax = config.getfloat("pmax", None)
    flatfield = config.get("flatfield", None)
    tag = config.get("tag", None)
    freq = config.getfloat("freq", None)

def run(p):
    global pmin, pmax, flatfield, tag, freq
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
    outfile = p.o.ff + "/ff_%s_cal.png" % tag
    print("Saving plot: %s" % outfile)
    v.array_plots(cal, det_uid, season=p.i.season, array=p.i.ar, fr=freq,
                  pmin=pmin, pmax=pmax, title='Flatfield %s' % tag,
                  display='save', save_name=outfile)

    # save stable detector plot
    outfile = p.o.ff + "/ff_%s_stable.png" % tag
    print("Saving plot: %s" % outfile)
    v.array_plots(stable, det_uid, title = 'Stable %s' % tag,
                  season=p.i.season, array=p.i.ar, display='save', fr=freq,
                  save_name=outfile)
