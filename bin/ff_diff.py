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
calibrate = False
freq = 150
v1 = "v1"
v2 = "v2"
pickle_file_1 = "/home/yguan/cuts/s17/c11/pa4_f150/pa4_f150_s17_c11_{}_results.pickle".format(v1)
pickle_file_2 = "/home/yguan/cuts/s17/c11/pa4_f150/pa4_f150_s17_c11_{}_results.pickle".format(v2)

ad = moby2.tod.ArrayData.from_fits_table('/home/lmaurin/actpol_data_shared/ArrayData/2017/ar4/default.fits')
dets = ad['det_uid'][ad['nom_freq']==freq]

with open(pickle_file_1, "rb") as f:
    data = pickle.load(f)

    # target to plot
    if calibrate:
        data['rmsLive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['normLive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['MFELive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['DELive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['jumpLive'] *= data['resp'] * data['ff'][:,np.newaxis]

    # calculate the average of the target patho param
    target_1 = data['ff']
    # target_m = np.std(target_values, axis=1, ddof=1)

    # array_plots(target_m[dets], dets, array="ar4", season="2017", display='save', 
    #            save_name=target+"_std.png", title=target+"_std", pmin=1e-10)

with open(pickle_file_2, "rb") as f:
    data = pickle.load(f)

    # target to plot
    if calibrate:
        data['rmsLive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['normLive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['MFELive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['DELive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['jumpLive'] *= data['resp'] * data['ff'][:,np.newaxis]

    # calculate the average of the target patho param
    target_2 = data['ff']

diff = np.abs((target_2 - target_1)) * 1.0 / target_1
name = "ff_diff_{}_{}".format(v1, v2)
array_plots(diff[dets], dets, array="ar4", season="2017", display='save', 
            save_name=name+".png", title=name, pmin=1e-10, pmax=0.1)
