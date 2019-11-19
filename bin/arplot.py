"""In this script I try to produce an array plot of certain pathological 
quantity. For example, it would be interesting to see the gain on an 
array plot, or the correlation on the array plot"""

import moby2
import pickle
import numpy as np
import numpy.ma as ma
from moby2.analysis.tod_ana.visual import array_plots
from matplotlib import pyplot as plt

# input file
calibrate = True
target = 'gainLive'
freq = 150
version = "v1"
pickle_file = "/home/yguan/cuts/s17/c11/pa4_f150/pa4_f150_s17_c11_{}_results.pickle".format(version)

ad = moby2.tod.ArrayData.from_fits_table('/home/lmaurin/actpol_data_shared/ArrayData/2017/ar4/default.fits')
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

    # calculate the average of the target patho param
    target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
    target_values *= data['ff'][:,None]
    target_m = np.std(target_values, axis=1)
    # target_m = data[target]
    # target_m = np.std(target_values, axis=1, ddof=1)

    # array_plots(target_m[dets], dets, array="ar4", season="2017", display='save', 
    #            save_name=target+"_std.png", title=target+"_std", pmin=1e-10)

    array_plots(target_m[dets], dets, array="ar4", season="2017", display='save', 
                save_name=target+"_{}_std.png".format(version), title=target, pmin=1e-10, pmax=0.2)
