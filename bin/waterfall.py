"""This script aims to generate a fft waterfall plot for a given
TOD"""

import moby2
import numpy as np
from matplotlib import pyplot as plt
import pickle
from moby2.analysis.tod_ana import visual as v

##############
# parameters #
##############

tod_name = "1500033549.1500091842.ar4"
pickle_file = "/home/yguan/cuts/s17/c11/pa4_f150/pa4_f150_s17_c11_v7_results.pickle"
depot = '/data/actpol/depot_yguan'
tag = 'pa4_f150_s17_c11_v7'
n_downsample = 10
season = '2017'
array = 'ar4'
freq = 150.

########
# main #
########
tag_partial = tag + "_partial"
fb = moby2.scripting.get_filebase()
filename = fb.filename_from_name(tod_name)[0]

# load data
tod = moby2.scripting.get_tod({
    'filename': filename,
    'repair_pointing': True,
})
moby2.tod.detrend_tod(tod)
moby2.tod.remove_mean(tod)

depot = moby2.util.Depot(depot_path=depot)
cuts = depot.read_object(moby2.TODCuts, tag=tag, tod=tod)
partial_cuts = depot.read_object(moby2.TODCuts, tag=tag_partial, tod=tod)
cal = depot.read_object(moby2.Calibration, tag=tag, tod=tod)

# calibrate
tod.data[cal.det_uid] *= cal.cal[:,None]
# fill cuts
moby2.tod.fill_cuts(tod, cuts, no_noise=True)
moby2.tod.fill_cuts(tod, partial_cuts, no_noise=True)

# manual fft 
# tdata = tod.data[:,::n_downsample]
# fdata = np.fft.fft(tdata, axis=1)
# freqs = np.fft.fftfreq(tdata.shape[1], 1./400*n_downsample)

# get the cal and cuts from pickle data
# idx = np.where(np.array(pickle_data['name']) == tod_name)[0][0]
# cuts_mask = pickle_data['sel'][:,idx]
# cuts_mask = cuts_mask.astype(bool)
# cal = pickle_data['cal'][:,idx]

# get the detector mask
# ardata = moby2.scripting.get_array_data({
#     'season': season,
#     'array_name': array
# })
# freq_mask = (ardata['nom_freq'] == freq)

# plot
# mask = np.logical_and(freq_mask, cuts_mask)
# plt.plot(freqs, np.abs(fdata)[mask, :].T, 'k', alpha=0.1)
# plt.xlabel('Freq (Hz)')
# plt.ylabel('FFT')
# plt.title(tod_name)
# plt.xscale('log')
# plt.yscale('log')
# plt.show()

ld = cuts.get_uncut()
obj = v.freqSpaceWaterfall(tod=tod, fmax=40)
plt.plot(obj.matfreqs, obj.mat[ld].T)
plt.show()

