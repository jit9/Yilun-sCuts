"""This script aims to produce a plot of live detectors vs. loading.
"""

import fitsio
import moby2, sys, numpy as np, ephem, os
from moby2.util.database import TODList
from moby2.instruments import actpol
from moby2.util import ctime as ct
from moby2.scripting.products import get_filebase
import pandas as pd
import matplotlib.pyplot as plt
import os.path as op
import seaborn as sns
from cutslib.pathologies_tools import pathoList, get_pwv


class Module:
    def __init__(self, config):
        # allow optionally specify an alternative cuts_db
        cuts_db = config.get('cuts_db', None)
        obs_catalog = config.get('obs_catalog', None)

    def run(self, p):
        cuts_db = self.cuts_db
        obs_catalog = self.obs_catalog

        if cuts_db is not None:
            pl = pathoList( cuts_db )
        else:
            cpar = moby2.util.MobyDict.from_file(p.i.cutparam)
            pl = pathoList( str(p.i.db) )
        #pl.addPWV2results()
        pl.removeDuplicates()
        pl.addEphem()
        Ndata = pl.ndata

        keys = ['todName', 'liveDets', 'hour', 'hourAfterSunset', 'hourAfterSunrise']
        PL = pd.DataFrame.from_dict( {k:pl.data[k] for k in keys} )

        filename = obs_catalog
        npcat = fitsio.read(filename)
        npcat = npcat.byteswap().newbyteorder()
        catalog = pd.DataFrame.from_records(npcat)
        catalog.index = pd.to_datetime(catalog.date)
        sel = np.logical_and( catalog.obs_type != 'stare', catalog.season == p.i.season)
        sel = np.logical_and( sel, catalog.array == p.i.ar)
        output = pd.merge(catalog[sel], PL, left_on='tod_name', right_on='todName', how='left')
        output.index = pd.to_datetime(output.ctime, unit='s')
        output.sort_index(inplace=True)
        output['PWV'] = output.pwv
        output.PWV[~np.isfinite(output.PWV)] = 0
        output['flag'] = np.zeros(len(output))
        output['loading_c'] = np.zeros(len(output))
        output.flag[~np.isnan(output.liveDets)] += 1
        ndets = output.liveDets.max()
        output = output[output.flag==1]
        # start to make the binned fraction plot
        # generate bins
        binsize = 0.25
        nbins = 16
        edges = binsize * np.arange(nbins+1)
        centers = (edges[1:] + edges[:-1])*1.0/2
        # calculate mean and std for each bin
        ms = np.zeros(nbins)
        ls = np.zeros(nbins)
        hs = np.zeros(nbins)
        means = np.zeros(nbins)
        stds = np.zeros(nbins)
        for i in range(nbins):
            bin_l, bin_h = edges[i], edges[i+1]
            # select tods inside the loading bin
            sel = np.logical_and(output.loading >= bin_l, output.loading < bin_h)
            l,m,h = output[sel].liveDets.quantile([0.25,0.5,0.75])
            print(l,m,h)
            ls[i] = l
            ms[i] = m
            hs[i] = h
            means[i] = output[sel].liveDets.mean()
            stds[i] = output[sel].liveDets.std(ddof=1)
        yerr = np.stack([ms-ls,hs-ms],axis=0)*100./ndets
        plt.figure(figsize=(15,8))
        plt.subplot(121)
        plt.errorbar(centers, ms*100./ndets, yerr=yerr, fmt='kx')
        plt.xlabel("loading / (pwv / sin(alt))")
        plt.ylabel("uncut percentage / \%")
        plt.title("25/50/75 quantiles")
        plt.subplot(122)
        plt.errorbar(centers, means*100./ndets, yerr=stds*100./ndets, fmt='kx')
        plt.xlabel("loading / (pwv / sin(alt))")
        plt.ylabel("uncut percentage / \%")
        plt.title("mean/std")
        outfile = op.join(p.o.root, 'ld_vs_loading.png')
        print("Writing figure: %s" % outfile)
        plt.savefig(outfile)
