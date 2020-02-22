"""This script plots all planet tods calibration components"""
import pickle as pickle
import numpy as np
import moby2
from cutslib import visual as v

def init(config):
    pass


def run(p):
    ff_dict = moby2.util.MobyDict.from_file(p.i.ff)
    dets = ff_dict['det_uid']

    with open(p.i.pickle_file, "rb") as f:
        data = pickle.load(f)

    resp = np.ma.array(data['resp'], mask=~data['sel'].astype(bool))
    gain = np.ma.array(data['gainLive'], mask=~data['sel'].astype(bool))
    for i in range(len(data['name'])):
        outfile = p.o.cal.resp+"/{}_resp.png".format(data['name'][i])
        print("Saving plot: %s" % outfile)
        v.array_plots(resp[dets,i], dets, season=p.i.season,
                      array=p.i.ar, title=data['name'][i] + " resp",
                      display='save', save_name=outfile, pmin=1e-16, pmax=2e-16)

        outfile = p.o.cal.gain_inv+"/{}_gain_inv.png".format(data['name'][i])
        print("Saving plot: %s" % outfile)
        v.array_plots((1./gain)[dets,i], dets, season=p.i.season,
                     array=p.i.ar, title=data['name'][i] + " gainLive^-1",
                     pmin=1e-10, display='save', save_name=outfile)
