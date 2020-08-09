#!/bin/env python

"""This script aims to patch bias-step failures in the form of an entire card
failed. It will find the previous timepoint that it didn't fail and use the
measurements from that instead. This will recover roughly 5% of data.

Example:
> python patch_bs.py --season 2017 --array pa4 --depot /projects/ACT/yilung/depot/biasstep/ --tag calibration_20200808 --ignore-bias-orig -V

"""
import numpy as np
import glob
from collections import Counter
from tqdm import tqdm
from cutslib import MobyDict
import moby2
import os.path as op
import pandas as pd
import pickle
import argparse
import os

#############
# interface #
#############

# define command-line arguments
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-D','--depot', help='depot for biassteps', default='/projects/ACT/feynman-users/spho/')
parser.add_argument('-S','--season', help='season to look at', required=True)
parser.add_argument('-A','--array', help='array to look at', required=True)
parser.add_argument('-T','--tag', help='tag to look at', default='calibration_20190802')
parser.add_argument('--ignore-bias-orig', help='whether to ignore DET_BIAS_ORG line', action="store_true")
parser.add_argument('--plot', help='whether to produce debug plot', action="store_true")
parser.add_argument('-F', '--force', help='whether to force overwriting cache', action="store_true")
parser.add_argument('--cache-name', help='filename for caching', default="cache_bsinfo.pkl")
parser.add_argument('-V', '--verbose', help='more verbose', action="store_true")
parser.add_argument('-W', '--write', help='write to file', action="store_true")
args = parser.parse_args()

#####################
# utility functions #
#####################

def find_fix_map(bcs):
    to_ = np.where(bcs==0)[0]
    from_ = to_ - 1
    dups = set(to_).intersection(from_)
    while len(dups) > 0:
        for dup in dups:
            from_[from_==dup] -= 1
        dups = set(to_).intersection(from_)
    return {'from':from_, 'to':to_}

def patch_bs(from_path, to_path, bc, season, array, write):
    """patch bias-step files

    Parameters
    ----------
    from_path: where to take bias-step from
    to_path: where to apply bias-step to
    bc: bias card to be patched (bc1, bc2, bc3)
    season/array: to get the array data
    write: flag to turn on actual writing

    """
    from_bs = MobyDict.from_file(from_path)
    to_bs = MobyDict.from_file(to_path)
    ad = moby2.scripting.get_array_data({'season':season,'array_name':array})
    bc_dets = np.where(ad['bias_card'] == bc)[0]
    # find alive dets in the given bias card by intersection
    to_add = list(set(from_bs['det_uid']).intersection(bc_dets))
    # find the index of each alive det and corresponding cal
    idx = [from_bs['det_uid'].index(d) for d in to_add]
    cal = [from_bs['cal'][i] for i in idx]
    calRMS = [from_bs['calRMS'][i] for i in idx]
    # add them to the target bs file
    to_bs['det_uid'] += to_add
    to_bs['cal'] += cal
    to_bs['calRMS'] += calRMS
    # sort list
    idx = np.argsort(to_bs['det_uid'])
    to_bs['det_uid'] = list(np.array(to_bs['det_uid'])[idx])
    to_bs['cal'] = list(np.array(to_bs['cal'])[idx])
    to_bs['calRMS'] = list(np.array(to_bs['calRMS'])[idx])
    if write: to_bs.write_to_file(to_path)

########
# main #
########

# search path
search_path = op.join(args.depot, args.season, args.tag, args.array, "*/*.cal")

# get array data for bias_card info
ad = moby2.scripting.get_array_data({'season':'2017','array_name':'ar4'})

# remove det_bias_org column from bias card 2 if necessary
if args.ignore_bias_orig:
    ad['bias_card'][ad['bias_name'] == 'DET_BIAS_ORG'] = 'none'

if not args.force and op.exists(args.cache_name):
    print("Use cached bias-step metadata")
    with open(args.cache_name, "rb") as f:
        data = pickle.load(f)
    bc1s = data['bc1s']
    bc2s = data['bc2s']
    bc3s = data['bc3s']
    ctimes = data['ctimes']
    ndets = data['ndets']
    fnames = data['fnames']
else:
    # collect metadata
    print("Collecting metadata")
    bc1s = []
    bc2s = []
    bc3s = []
    ctimes = []
    ndets = []
    fnames = []
    for fn in tqdm(glob.glob(search_path)):
        data = MobyDict.from_file(fn)
        ctr = Counter(ad['bias_card'][data['det_uid']])
        bc1 = bc2 = bc3 = 0
        if 'bc1' in ctr: bc1 = ctr['bc1']
        if 'bc2' in ctr: bc2 = ctr['bc2']
        if 'bc3' in ctr: bc3 = ctr['bc3']
        ctime = int(op.basename(fn).split('.')[0])
        bc1s.append(bc1)
        bc2s.append(bc2)
        bc3s.append(bc3)
        ctimes.append(ctime)
        ndets.append(len(data['det_uid']))
        fnames.append(fn)
    # sort by ctime
    idx = np.argsort(ctimes)
    bc1s = np.asarray(bc1s)[idx]
    bc2s = np.asarray(bc2s)[idx]
    bc3s = np.asarray(bc3s)[idx]
    ctimes = np.asarray(ctimes)[idx]
    ndets = np.asarray(ndets)[idx]
    fnames = np.asarray(fnames)[idx]
    # cache
    with open(args.cache_name, "wb") as f:
        to_save = {
            'bc1s': bc1s, 'bc2s': bc2s, 'bc3s': bc3s,
            'ctimes': ctimes, 'ndets': ndets, 'fnames': fnames
        }
        pickle.dump(to_save, f)
        print("Cache saved")

# produce plots
if args.plot:
    print("Generating debug plots")
    from matplotlib import pyplot as plt
    # det count plot
    plt.figure(figsize=(8,6))
    plt.scatter(ctimes, ndets, alpha=0.5, s=2, c='k')
    plt.scatter(ctimes, bc1s, alpha=0.5, s=2, c='r', label='bc1')
    plt.scatter(ctimes, bc2s, alpha=0.5, s=2, c='g', label='bc2')
    plt.scatter(ctimes, bc3s, alpha=0.5, s=2, c='b', label='bc3')
    plt.ylabel('# of dets')
    plt.xlabel('ctime')
    plt.legend(loc='upper left')
    plt.savefig("debug_det_count.png")
    plt.close()
    # bias card failure plot
    plt.figure(figsize=(8,6))
    plt.plot(ctimes, np.cumsum(bc1s == 0), alpha=0.5, c='r', label='bc1')
    plt.plot(ctimes, np.cumsum(bc2s == 0), alpha=0.5, c='g', label='bc2')
    plt.plot(ctimes, np.cumsum(bc3s == 0), alpha=0.5, c='b', label='bc3')
    plt.ylabel('Accumulation number of Biascard failure')
    plt.xlabel('ctime')
    plt.legend()
    plt.savefig("debug_acc_bc_fail.png")
    plt.close()

# come up with fix table
fixes = {
    'bc1':find_fix_map(bc1s),
    'bc2':find_fix_map(bc2s),
    'bc3':find_fix_map(bc3s)
}
fixes_df = []
for bc in fixes:
    for from_, to_ in zip(fixes[bc]['from'], fixes[bc]['to']):
        fixes_df.append({
            'from_ctime': ctimes[from_],
            'to_ctime': ctimes[to_],
            'from_path': fnames[from_],
            'to_path': fnames[to_],
            'bc_fix': bc
        })
fixes_df = pd.DataFrame(fixes_df)
outfile = f"fix_{args.season}_{args.array}.csv"
fixes_df.to_csv(outfile, index=False)

# show the fix table
print("Showing fix table")
os.system(f"cat {outfile} | less")
if input("proceed to fix? y/n ") == 'y':
    print("Fixing...")
    if not args.write: print("write flag off: performing a dry-run")
    for i, row in fixes_df.iterrows():
        if args.verbose: print(f"{row['bc_fix']}: {row['from_ctime']} -> {row['to_ctime']}")
        patch_bs(row['from_path'],row['to_path'],row['bc_fix'],args.season,args.array, args.write)
else:
    print("Fixing...")
    import time
    time.sleep(1)
    print("Just kidding :)")
