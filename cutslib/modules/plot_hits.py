"""This script plots the hits map of the source list"""

import numpy as np
import matplotlib.pyplot as plt
import moby2
from moby2.util.database import TODList
import os.path as op
import os, pickle

def init(config):
    global exclude, ra_l, ra_h, dec_l, dec_h, ra_bins, dec_bins, force
    exclude = config.get("exclude_list")
    ra_l = config.getfloat("ra_l",-2)
    ra_h = config.getfloat("ra_h",2)
    dec_l = config.getfloat("dec_l",-0.4)
    dec_h = config.getfloat("dec_h",0.4)
    ra_bins = config.getint("ra_bins",100)
    dec_bins = config.getint("dec_bins",100)
    force = config.getboolean("force", False)

def run(p):
    global exclude, ra_l, ra_h, dec_l, dec_h, ra_bins, dec_bins, force
    # check if pickle file is present
    pkl_file = op.join(p.o.root, 'hits.pkl')
    if not os.path.exists(pkl_file) or force:
        # get list of tods
        par = moby2.util.MobyDict.from_file(p.i.cutparam)
        tod_list = TODList.from_file(par.get('source_scans'))
        print("Found %d TODs" % len(tod_list))
        exclude_list = TODList.from_file(exclude)
        print("Excluding %d TODs" % len(exclude_list))
        tod_list -= exclude_list
        # load each tod and get the ra and dec for each tod
        # setup histogram
        x = np.linspace(ra_l,ra_h,ra_bins)
        y = np.linspace(dec_l,dec_h,dec_bins)
        X, Y = np.meshgrid(x, y)
        first = True
        for i,n in enumerate(tod_list):
            print("Processing [%d/%d]: %s" % (i,len(tod_list),n))
            tod = moby2.scripting.get_tod({'filename': n, 'repair_pointing': True, 'read_data': False})
            # check scan error
            error = np.sum(np.logical_or(tod.alt > np.pi/2, tod.alt < 0)) > 0
            if error:
                print("Alt falls outside range of 0 and pi/2, skipping")
                continue

            ra, dec = moby2.pointing.get_coords(tod.ctime,tod.az,tod.alt)
            if first:  # use the first as a basis and add others on it
                H, _, __ = np.histogram2d(ra, dec, bins=(x, y))
                H = H.T
                first = False
            else:  # add other tods onto this histogram
                H_new, _, __ = np.histogram2d(ra, dec, bins=(x, y))
                H_new = H_new.T
                H += H_new
        # convert to degrees
        X = X / np.pi * 180.
        Y = Y / np.pi * 180.
        # dump the data into pickle file
        hits_dict = {
            'X': X,
            'Y': Y,
            'H': H,
        }
        outfile = op.join(p.o.root, 'hits.pkl')
        print("Saving data: %s" % outfile)
        with open(outfile, "w") as f:
            pickle.dump(hits_dict, f)
    else:  # if pickle file exists
        print("Found existing: %s" % pkl_file)
        with open(pkl_file, "r") as f:
            hits_dict = pickle.load(f)
        X,Y,H = hits_dict['X'], hits_dict['Y'], hits_dict['H']
    # plot hist
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(X,Y,H, cmap="hot")
    plt.colorbar(mesh)
    ax.set_xlabel("Ra / degree")
    ax.set_ylabel("Dec / degree")
    plt.gca().invert_xaxis()
    outfile = op.join(p.o.root, 'hits.png')
    print("Saving plot: %s" % outfile)
    plt.savefig(outfile)
