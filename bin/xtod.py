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
pickle_file = "/home/yguan/cuts/s17/c11/pa4_f150/pa4_f150_s17_c11_v0_results.pickle"

# target to plot
target = 'DELive'
freq = 150

ad = moby2.tod.ArrayData.from_fits_table('/home/lmaurin/actpol_data_shared/ArrayData/2017/ar4/default.fits')
dets = ad['det_uid'][ad['nom_freq']==freq]

with open(pickle_file, "rb") as f:
    data = pickle.load(f)

    # calculate the average of the target patho param
    target_values = ma.array(data[target], mask=np.logical_not(data['sel']))
    target_m = np.mean(target_values, axis=0)
    times = [int(tod_name.split('.')[0]) for tod_name in data['name']]
    
    plt.plot(times, target_m, 'kx')
    plt.xlabel("ctime")
    plt.ylabel(target)
    plt.title("{} vs time".format(target))
    plt.savefig("{}_ts.png".format(target))
    plt.show()
    plt.close()

