"""This script generates detector list for a given array and season"""

import argparse, numpy as np
import os, os.path as op
import moby2
parser = argparse.ArgumentParser()
parser.add_argument("--season")
parser.add_argument("--array")
parser.add_argument("-o", "--odir", default="out")
parser.add_argument("--depot", help="actpol shared depots", default="/home/yguan/actpol_data_shared")
parser.add_argument("--suffix", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.suffix: suffix = f"_{args.suffix}"
else: suffix = ""

def get_data_path(season, array, shared_depot=None):
    """the idea is to load the fits table from a given season and array as input"""
    season = season.lower().replace('s','20')
    array = array.lower().replace('pa','ar')
    if not shared_depot: shared_depot = moby2.user_cfg.get_deep(('depots','actpol_shared','path'))
    return op.join(shared_depot, 'ArrayData', season, array, 'default.fits')

# get fits file
fpath = get_data_path(args.season, args.array, args.depot)
ad = moby2.util.MobyFitsTable.from_fits(fpath)

# assume that we are working with dichoric arrays, we will find the
# two frequency codes
fcodes = np.unique(ad.fcode)
# remove the dummy fcode 'f000'
fcodes = [fc for fc in fcodes if fc != 'f000']
print("Found fcodes:", fcodes)

for fc in fcodes:
    # first, make a list of live detectors: this means all tes detectors
    live = ad.det_uid[(ad.det_type == 'tes') * (ad.fcode == fc)]
    # next, make a list of dark detectors: this means all dark_tes detectors
    dark = ad.det_uid[(ad.det_type == 'dark_tes') * (ad.fcode == fc)]

    # make sure none of the list has length of 0, otherwise call for help
    for name, dets in zip(['live', 'dark'], [live, dark]):
        print(f"Found {name}:{fc}: {len(dets):5d}")
        assert len(dets) > 0, f"No detectors found for {name}, something is wrong!"

    # now it looks like everything is okay, generating the dictionary of detector lists
    for name, dets in zip(['live', 'dark'], [live, dark]):
        ofile = op.join(args.odir, f'{name}_{fc}{suffix}.dict')
        print("Writing:", ofile)
        moby2.util.MobyDict({'det_uid': list(dets)}).write(ofile)

# last write exclude list which are NC detectors by default
exclude = ad.det_uid[ad.det_type == 'NC']
ofile = op.join(args.odir, f'exclude{suffix}.dict')
print("Writing:", ofile)
moby2.util.MobyDict({'det_uid': list(exclude)}).write(ofile)
