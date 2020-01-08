"""This script aims to provide a more advanced tod selection
algorithm that selects a particular region and include the uranus
tods, it allows trimming the tods based on ra and dec.

Example:
>> python select_tod_adv.py tod.par

An example parameter file is given by:
# =========================================
output_dir = "/scratch/gpfs/yilung/depot/TODLists"
tag = '{season}_{array}_wide01hn'
season = '2017'
array = 'ar6'
obs_detail = ['wide_01h_n']
hwp_epoch = False
night = 1
pwv_upper = 2
ra_lower = -1
ra_upper = 1
calibration = 'uranus'
nlimit = 1000
# ==========================================
"""

import moby2
import sys, os
import fitsio
import pandas as pd
import numpy as np
from moby2.ephem import ACTEphem
import matplotlib.pyplot as plt

fb = moby2.scripting.get_filebase()

# utility function
def write_list_to_file(output_dir, filename, lst):
    """A utility function to save tod list into disk"""
    print('Output %i tods to %s' %(len(lst), filename))
    with open(os.path.join(output_dir, filename), "w") as f:
        for n in lst:
            f.write("%s\n" % n)

# load parameter file
params = moby2.util.MobyDict.from_file( sys.argv[1] )
output_dir = params.get("output_dir", ".")

# load catalog of tods and convert to pandas dataframe
moby2_conf = moby2.util.MobyConfig()
filename = moby2_conf.get_deep(('obs_catalog','filename'))

npcat = fitsio.read(filename)
npcat = npcat.byteswap().newbyteorder()
catalog = pd.DataFrame.from_records(npcat)
catalog.index = pd.to_datetime(catalog.date)

# select tods based on the specified criteria
# 1st iteration
sel = np.ones(catalog.shape[0], dtype=bool)
for k in params.keys():
    if k in catalog.columns:
        p = params[k]
        if type(p) == list:
            sel *= catalog[k].isin(p)
        elif type(p) == bool:
            if p:
                sel *= catalog[k] != 'none'
            else:
                sel *= catalog[k] == 'none'
        else:
            sel *= catalog[k] == p
        print('%s: left with %i TODs' % (k, sel.sum()))
    # else:
    #     print('%s: not a valid parameter, skip it' %k)

# pre-selection
catalog = catalog[sel]

# get ra and dec for each tod
def compute_radec(row):
    ctime = row['ctime']
    alt = row['alt']
    az = row['az']
    # load ephem
    eph = ACTEphem()
    eph.set_ctime(ctime)
    ra, dec = eph.altaz_to_radec(alt, az)
    if ra>np.pi: ra -= 2*np.pi  # make sure range is -pi to pi
    return ra, dec

# check whether each entry is at night
def is_night(row):
    """night is defined to be CLT 10pm to 10am. This function returns
    1 for TOD observed at night and 0 for day. Note that UTC time has
    4 hour difference with CLT time"""
    hour_utc = row['hour_utc']
    if hour_utc >= 2 and hour_utc <= 14:
        return 1
    else:
        return 0

ra = catalog.apply(lambda row: compute_radec(row)[0], axis=1)
dec = catalog.apply(lambda row: compute_radec(row)[1], axis=1)
night = catalog.apply(lambda row: is_night(row), axis=1)
catalog['ra'] = ra
catalog['dec'] = dec
catalog['night'] = night

# second iteration selection
# select tods based on the specified criteria
sel = np.ones(catalog.shape[0], dtype=bool)

for k in params.keys():
    if k in ['ra','dec','night']:
        p = params[k]
        if type(p) == list:
            sel *= catalog[k].isin(p)
        elif type(p) == bool:
            if p:
                sel *= catalog[k] != 'none'
            else:
                sel *= catalog[k] == 'none'
        else:
            sel *= catalog[k] == p
        print('%s: left with %i TODs' % (k, sel.sum()))

catalog = catalog[sel]

# third iteration selection
# upper and lower are both closed interval
for k in ['pwv', 'ra', 'dec']:
    k_upper = "%s_upper" % k
    if k_upper in params.keys():
        v_upper = params.get(k_upper)
        sel = catalog[k] <= v_upper
        catalog = catalog[sel]
        print("%s: left with %d TODs" % (k_upper, sel.sum()))
    k_lower = "%s_lower" % k
    if k_lower in params.keys():
        v_lower = params.get(k_lower)
        sel = catalog[k] >= v_lower
        catalog = catalog[sel]
        print("%s: left with %d TODs" % (k_lower, sel.sum()))

# compile list of tods
tod_list = list(catalog['tod_name'])

if 0:
    plt.plot(catalog['ra'], catalog['dec'], 'r.')
    plt.savefig("test.png")

# load each tod and get the ra and dec for each tod
if 0:
    # setup histogram
    x = np.linspace(-2,2.5,100)
    y = np.linspace(-0.4,0.4,100)
    X, Y = np.meshgrid(x, y)
    first = True
    for n in tod_list:
        print("Processing %s" % n)
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

    fig, ax = plt.subplots()
    ax.pcolormesh(X,Y,H)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    plt.savefig("hits.png")

# 4th iteration: look at uranus
if "calibration" in params.keys():
    planet = params.get("calibration")
    print("include %s calibration tods..." % planet)
    catalog = pd.DataFrame.from_records(npcat)
    catalog.index = pd.to_datetime(catalog.date)
    # select tods based on the specified criteria
    # 1st iteration
    sel = np.ones(catalog.shape[0], dtype=bool)
    for k in params.keys():
        if k in ['season', 'array', 'hwp_epoch']:
            p = params[k]
            if type(p) == list:
                sel *= catalog[k].isin(p)
            elif type(p) == bool:
                if p:
                    sel *= catalog[k] != 'none'
                else:
                    sel *= catalog[k] == 'none'
            else:
                sel *= catalog[k] == p
            print('%s: left with %i TODs' % (k, sel.sum()))
    # choose the given planet
    sel *= catalog['obs_detail'] == planet
    print('%s: left with %i TODs' % (planet, sel.sum()))
    catalog = catalog[sel]
    new_tod_list = list(catalog['tod_name'])
    # save the uranus tods separately
    output_filename = params['tag'].format(**params)+"_%s.txt" % planet
    write_list_to_file(output_dir, output_filename, new_tod_list)
    print("Add %s tods: %d" % (planet, len(new_tod_list)))
    tod_list.extend(new_tod_list)

print("Total number of TODs selected: %d" % len(tod_list))

if "nlimit" in params.keys():
    nlimit = params.get("nlimit")
    print("nlimit set: %d" % nlimit)
    tod_list = tod_list[:nlimit]  # choose the last nlimit ensures
                                  # that calibration tods are included
    print("Total number of TODs selected: %d" % len(tod_list))

# save the list
output_filename = params['tag'].format(**params)+".txt"
write_list_to_file(output_dir, output_filename, tod_list)


print("Done")
######################################################################
# end of select_tod_adv.py
