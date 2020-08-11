#!/bin/env python

"""This script aims to patch bias-step failures in the form of an entire card
failed. It will find the previous timepoint that it didn't fail and use the
measurements from that instead. This will recover roughly 5% of data.

Example:
> python patch_bs.py --season 2017 --array pa4 --depot /projects/ACT/yilung/depot/biasstep/ --tag calibration_20200808 -V --plot

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
parser.add_argument('--include-bias-orig', help='whether to include DET_BIAS_ORG line', action="store_true")
parser.add_argument('--include-pxl-htr', help='whether to include PXL_HTR line', action="store_true")
parser.add_argument('--plot', help='whether to produce debug plot', action="store_true")
parser.add_argument('-F', '--force', help='whether to force overwriting cache', action="store_true")
parser.add_argument('--cache-name', help='filename for caching, no extension', default="cache_bsinfo")
parser.add_argument('-V', '--verbose', help='more verbose', action="store_true")
parser.add_argument('-W', '--write', help='write to file', action="store_true")
parser.add_argument('-Y', '--yes', help='skip asking y/n and proceed with yes', action="store_true")
args = parser.parse_args()

#####################
# utility functions #
#####################

def find_fix_map(ctimes, ivtimes, bcs):
    to_ = np.where(bcs==0)[0]
    from_ = []
    for i in tqdm(to_):
        # get the iv time for the biasstep to patch
        ivt = ivtimes[i]
        # sort biasstep time from nearest to furthest
        idx = np.argsort(np.abs(ctimes - ctimes[i]))
        # narrow down to biassteps taken under the same iv
        idx = [i_ for i_ in idx if i_ in np.where(ivtimes == ivt)[0]]
        # if we can't find any, indicate with a None
        to_append = -1
        if len(idx) > 1:
            for j in range(1,len(idx)):
                if idx[j] not in to_: to_append=idx[j]; break
        from_.append(to_append)
    return {'from':np.asarray(from_), 'to':np.asarray(to_)}

def find_iv(ctimes, ivtimes):
    """find iv to use for each time in given in the ctimes list

    Parameters
    ----------
    ctimes: list / array of ctimes when bias-step are taken
    ivtimes: list / array of times when iv are taken

    Returns
    -------
    array of iv times to use for each time in ctimes list

    """
    ivs = np.zeros_like(ctimes)
    for i in range(len(ctimes)):
        # find nearest preceeding iv to use
        try: ivs[i] = ivtimes[np.where(ivtimes < ctimes[i])[0][-1]]
        except: ivs[i] = -1
    return ivs

def patch_bs(from_path, to_path, bc, ad, write):
    """patch bias-step files

    Parameters
    ----------
    from_path: where to take bias-step from
    to_path: where to apply bias-step to
    bc: bias card to be patched (bc1, bc2, bc3)
    ad: array_data dictionary from moby2
    write: flag to turn on actual writing

    """
    from_bs = MobyDict.from_file(from_path)
    to_bs = MobyDict.from_file(to_path)
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

# ivtimes lookup table
iv_lookup = {
    '2017_pa4': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S5_mce1',
    '2017_pa5': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S5_mce2',
    '2017_pa6': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S5_mce3',
    '2018_pa4': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S6_mce1',
    '2018_pa5': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S6_mce2',
    '2018_pa6': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S6_mce3',
    '2019_pa4': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S7_mce1',
    '2019_pa5': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S7_mce2',
    '2019_pa6': '/projects/ACT/yilung/depot/biasstep/IV_ctime_S7_mce3',
}
ivtimes = np.sort(np.loadtxt(iv_lookup[f"{args.season}_{args.array}"]))

# get array data for bias_card info
ad = moby2.scripting.get_array_data({'season':args.season,'array_name':args.array})

# remove det_bias_org column from counting if necessary
if not args.include_bias_orig:
    ad['bias_card'][ad['bias_name'] == 'DET_BIAS_ORG'] = 'none'
# remove PXL_HTR column from counting if necessary
if not args.include_pxl_htr:
    ad['bias_card'][ad['bias_name'] == 'PXL_HTR'] = 'none'

cache_name =f"{args.cache_name}_{args.season}_{args.array}.pkl"
if not args.force and op.exists(cache_name):
    print("Use cached bias-step metadata")
    with open(cache_name, "rb") as f:
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
    with open(cache_name, "wb") as f:
        to_save = {
            'bc1s': bc1s, 'bc2s': bc2s, 'bc3s': bc3s,
            'ctimes': ctimes, 'ndets': ndets, 'fnames': fnames
        }
        pickle.dump(to_save, f)
        print("Cache saved")

# match ivtimes to biasstep times
ivtimes = find_iv(ctimes, ivtimes)
# skip bias-steps without iv to use
mask = ivtimes != -1
bc1s = bc1s[mask]
bc2s = bc2s[mask]
bc3s = bc3s[mask]
ctimes = ctimes[mask]
ivtimes = ivtimes[mask]
ndets = ndets[mask]
fnames = fnames[mask]

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
    plt.savefig(f"{args.season}_{args.array}_det_count.png")
    plt.close()
    # bias card failure plot
    plt.figure(figsize=(8,6))
    plt.plot(ctimes, np.cumsum(bc1s == 0), alpha=0.5, c='r', label='bc1')
    plt.plot(ctimes, np.cumsum(bc2s == 0), alpha=0.5, c='g', label='bc2')
    plt.plot(ctimes, np.cumsum(bc3s == 0), alpha=0.5, c='b', label='bc3')
    plt.ylabel('Accumulation number of Biascard failure')
    plt.xlabel('ctime')
    plt.legend()
    plt.savefig(f"{args.season}_{args.array}_acc_bc_fail.png")
    plt.close()

# come up with fix table
fixes = {
    'bc1':find_fix_map(ctimes, ivtimes, bc1s),
    'bc2':find_fix_map(ctimes, ivtimes, bc2s),
    'bc3':find_fix_map(ctimes, ivtimes, bc3s)
}

fixes_df = []
for bc in fixes:
    for from_, to_ in zip(fixes[bc]['from'], fixes[bc]['to']):
        if from_ == -1: continue
        fixes_df.append({
            'from_ctime': ctimes[from_],
            'to_ctime': ctimes[to_],
            'from_path': fnames[from_],
            'to_path': fnames[to_],
            'bc_fix': bc
        })
fixes_df = pd.DataFrame(fixes_df)
outfile = f"fix_{args.season}_{args.array}.csv"
fixes_df.sort_values(by=['to_ctime']).to_csv(outfile, index=False, columns=['to_ctime','from_ctime','bc_fix'], sep=' ', na_rep='-')

# show the fix table
print("Showing fix table")
if not args.yes: os.system(f"cat {outfile} | less")
if args.yes or input("proceed to fix? y/n ") == 'y':
    print("Start fixing...")
    if not args.write: print("write flag off: performing a dry-run")
    for i, row in fixes_df.iterrows():
        if args.verbose: print(f"{row['bc_fix']}: {row['from_ctime']} -> {row['to_ctime']}")
        patch_bs(row['from_path'], row['to_path'], row['bc_fix'], ad, args.write)
else:
    print("Fixing...")
    import time
    time.sleep(1)
    print("Just kidding :)")
