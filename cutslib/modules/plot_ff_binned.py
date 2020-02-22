"""This module wraps around generateFlatfield script in moby2 to
create a flatfield dictionary that can be used as the input to the
next iteration of cuts
"""

import fitsio
import os.path as op
import pickle, numpy as np, os, sys
import moby2
import pandas as pd
from cutslib.visual import array_plots


def get_pwv(tod_list, obs_catalog):
    # load catalog to get loadings
    npcat = fitsio.read(obs_catalog)
    npcat = npcat.byteswap().newbyteorder()
    catalog = pd.DataFrame.from_records(npcat)
    catalog.index = pd.to_datetime(catalog.date)
    loadings = []
    for d in tod_list:
        l = catalog[catalog.tod_name==d].loading.values
        assert len(l) == 1  # sanity check
        loadings.append(l)
    loadings = np.hstack(loadings)
    return loadings


def getArrayStats(gains, sels, selTODs, use_sel=True, gain_limit=10,
                  min_samples=50, sigmas=5):
    sels = np.array(sels, dtype=bool)
    Ndets = gains.shape[0]
    means = np.zeros(Ndets)
    stds = np.zeros(Ndets)
    for d in range(Ndets):
        calib = 1./np.array(gains[d])
        calib[np.isnan(calib)] = 0.0
        sel = (np.abs(calib) < gain_limit)*(np.abs(calib) > 1./gain_limit)*selTODs
        if use_sel and len(sels[d]) > 0: sel *= np.array(sels[d])
        if np.array(sel, dtype = int).sum() < min_samples:
            # print("Failed minSample test %d" % d)
            continue
        tmp = np.sort(calib[sel])
        n = len(calib[sel])
        m = tmp[int(n/2)]
        q25 = tmp[int(n/4)]
        q75 = tmp[int((3*n)/4)]
        s = 0.741*(q75-q25)
        sel = (calib < m+sigmas*s)*(calib > m-sigmas*s)*(calib > 0)
        means[d] = np.median(calib[sel])
        stds[d] = calib[sel].std()
    return means, stds


def init(config):
    global obs_catalog, gain_limit, sigmas, min_samples, use_sel, pmin, pmax
    obs_catalog = config.get("obs_catalog", None)
    gain_limit = config.getfloat("gain_limit", 10.)
    sigmas = config.getfloat("sigmas", 5.)
    min_samples = config.getint("min_samples", 50)
    use_sel = config.getboolean("use_sel", True)
    pmin = config.getfloat("pmin", 0.6)
    pmax = config.getfloat("pmax", 1.3)

def run(p):
    global obs_catalog, gain_limit, sigmas, min_samples, use_sel, pmin, pmax
    # load parameters
    params = moby2.util.MobyDict.from_file(p.i.cutparam)
    cutParam = moby2.util.MobyDict.from_file(p.i.cutParam)
    ffpar = params["ff_params"]
    ff_name = cutParam.get_deep(('pathologyParams','calibration','flatfield'))
    ff = moby2.util.MobyDict.from_file(ff_name)
    # get det_uid from input flatfield
    det_uid = np.asarray(ff['det_uid'], dtype=int)
    # load pickle file
    print('Loading data')
    with open(p.i.pickle_file, "rb") as f:
        pf = pickle.Unpickler(f)
        data = pf.load()
    # get pwv for each of the tods
    loadings = get_pwv(data['name'], obs_catalog)
    gains = data["gainLive"].copy()
    sel = np.asarray(data['sel'],dtype=bool)*np.asarray(data['respSel'],dtype=bool)
    # create bins of loadings
    # [0,1), [1,2), [2,3), [3,4)
    nbins = 4
    lbins = np.arange(nbins+1)
    lbins_l = lbins[:-1]
    lbins_h = lbins[1:]
    # produce flatfield for each bin
    for i in range(nbins):
        # get bin edges
        bin_l, bin_h = lbins_l[i], lbins_h[i]
        # get list of tods inside this bin
        tod_sel = (loadings >= bin_l) * (loadings < bin_h)
        m, s = getArrayStats(gains, sel, tod_sel, use_sel, gain_limit,
                             min_samples, sigmas)
        # set color range
        if not pmin:
            pmin = m.min()
        if not pmax:
            pmax = m.max()
        # plot flatfield on the array and save figure
        outfile = op.join(p.o.ff, "ff_binned_%d.png" % i)
        print("Saving plot: %s" % outfile)
        array_plots(m[det_uid],det_uid,season=p.i.season,array=p.i.ar,fr=p.i.freq,
                    pmin=pmin,pmax=pmax,display='save',save_name=outfile,
                    title="Flatfield: loading bin [%.1f,%.1f)" % (bin_l,bin_h))
