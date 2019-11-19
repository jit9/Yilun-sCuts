"""This module plots the flatfield on an array"""

import os
import moby2.util
import numpy as np
from moby2.analysis.tod_ana import visual as v

def init(config):
    global pmin, pmax
    pmin = config.getfloat("pmin", None)
    pmax = config.getfloat("pmax", None)

def run(p):
    global pmin, pmax
    ff_name = p.i.ff
    ff = moby2.util.MobyDict.from_file(ff_name)

    det_uid = np.asarray(ff['det_uid'], dtype=int)
    cal = np.asarray(ff['cal'], dtype=float)
    stable = np.asarray(ff['stable'], dtype=int)

    if pmin is None:
        pmin = cal.min()
    if pmax is None:
        pmax = cal.max()

    # save flatfield plot
    outfile = p.o.ff + "/ff_%s_cal.png" % p.tag
    print("Saving plot: %s" % outfile)
    v.array_plots(cal, det_uid, season=p.i.season, array=p.i.ar,
                  pmin=pmin, pmax=pmax, title='Flatfield %s' % p.tag,
                  display='save', save_name=outfile)

    # save stable detector plot
    outfile = p.o.ff + "/ff_%s_stable.png" % p.tag
    print("Saving plot: %s" % outfile)
    v.array_plots(stable, det_uid, title = 'Stable %s' % p.tag,
                  season=p.i.season, array=p.i.ar, display='save',
                  save_name=outfile)


