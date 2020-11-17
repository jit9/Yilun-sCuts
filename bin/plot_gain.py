"""This script compiles gain as a function of time for each detector and
save them as data file for furthur analysis

Example:
python plot_gain.py --tag pa4_f150_s19_c11_v1 --nrand 20 --crange 1555000000:1556000000 --hour -v
"""

import argparse, os, os.path as op
import numpy as np
import matplotlib.pyplot as plt

from cutslib import SeasonStats
from cutslib.glitch import PixelReader

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--odir', default='out')
parser.add_argument('-v', action='store_true', default=False)
parser.add_argument('--depot', default=None)
parser.add_argument('--tag')
parser.add_argument('--dets', help='det_uids in comma seperated format, default to all', default=None)
parser.add_argument('--nrand', help='choose n random dets', type=int, default=None)
parser.add_argument('--crange', help='range of ctimes seperated by colon', default=None)
parser.add_argument('--hour', help='use hour since beginning', action='store_true', default=False)
parser.add_argument('--rmax', help="maximum A/B ratio to plot", type=int, default=3)
args = parser.parse_args()
# args preprocessing
if not op.exists(args.odir): os.makedirs(args.odir)

# first load the season stats from the given tag and depot
if args.v: print("Loading season stats...")
ss = SeasonStats(tag=args.tag, depot=args.depot, verbose=args.v, planet=False, rundb=False)

# default list of detector
if args.dets: dets = np.array(args.dets.split(',')).astype(int)
else:
    dets = np.where(ss.ff_sel * ss.tes_sel)[0]
    if args.nrand: dets = np.random.choice(dets, args.nrand, replace=False)
# plot calibration with uncut dets
cal = np.ma.array(np.abs(ss.cal))
cal[~ss.sel] = np.ma.masked
# apply absolute calibration per tod to account for pwv effects
cal *= ss.abscal[None,:]
# xaxis
if args.hour:
    tfunc = lambda x: (x - ss.ctime.min())/3600
    ctime = tfunc(ss.ctime)
    xlabel = f"hours since {ss.ctime.min()} [h]"
else:
    ctime = ss.ctime
    xlabel = "ctime [s]"
lines = plt.plot(ctime, cal[dets].T, '.', markersize=1)
plt.xlabel(xlabel)
plt.ylabel('calibration [uK/DAC]')
# legend
if len(dets)<=20:
    plt.legend(iter(lines), dets, bbox_to_anchor=(1.1,1),
               loc="upper left", ncol=int(np.ceil(len(dets)/10)))
# also plot pwv
pwv_ax = plt.gca().twinx()
idx = np.argsort(ctime)
pwv_ax.plot(ctime[idx], ss.pwv[idx], 'k-', alpha=0.2)
pwv_ax.set_ylabel('pwv / sin(alt) [mm]')
# ctime range
if args.crange:
    cstart, cend = args.crange.split(':')
    cstart = max(ss.ctime.min(), int(cstart))
    cend   = min(ss.ctime.max(), int(cend))
    if args.hour: cstart, cend = tfunc(cstart), tfunc(cend)
    plt.xlim([cstart, cend])
plt.title('calibration=ff*biasstep*abscal')
ofile = op.join(args.odir, f'gain_{args.tag}.png')
if args.v: print("Saving:", ofile)
plt.savefig(ofile, bbox_inches='tight')
plt.close()

# next make a plot of A/B relative calibrations
# for all pairs (for now)
dets = np.where(ss.ff_sel * ss.tes_sel)[0]
pr = PixelReader(season=ss.season, array=ss.array)
# find all pixels that have both detectors listed in dets
pf1 = [pr.get_f1(p) for p in pr.get_pixels()
       if np.sum(np.isin(pr.get_f1(p), dets))==2]
pf2 = [pr.get_f2(p) for p in pr.get_pixels()
       if np.sum(np.isin(pr.get_f2(p), dets))==2]
# one of f1 and f2 will have nonempty list corresponding
# to the freq of interests
if args.v: print(f"Found pairs: f1 [{len(pf1)}] f2 [{len(pf2)}]")
pairs = np.array(pf1) if len(pf1) > 0 else np.array(pf2)
if args.nrand:
    ridx = np.random.choice(np.arange(len(pairs)), args.nrand, replace=False)
    pairs = pairs[ridx]
for pair in pairs:
    AB_ratio = cal[pair[0]] / cal[pair[1]]
    if np.all(AB_ratio < 1): AB_ratio = 1/AB_ratio
    plt.plot(ctime, AB_ratio, '.', markersize=1, label=f"{pair[0]}-{pair[1]}")
plt.ylim([0,args.rmax])
plt.xlabel(xlabel)
plt.ylabel('A/B ratio')
plt.legend(bbox_to_anchor=(1.1,1), loc="upper left", ncol=int(np.ceil(len(pairs)/10)))
# ctime range
if args.crange: plt.xlim([cstart, cend])
plt.title('calibration A/B ratio (>1)')
# also plot pwv
pwv_ax = plt.gca().twinx()
idx = np.argsort(ctime)
pwv_line = pwv_ax.plot(ctime[idx], ss.pwv[idx], 'k-', alpha=0.2, label='pwv/sin(alt)')
pwv_ax.set_ylabel('pwv / sin(alt) [mm]')
plt.legend(loc='best')
ofile = op.join(args.odir, f'gain_AB_{args.tag}.png')
if args.v: print("Saving:", ofile)
plt.savefig(ofile, bbox_inches='tight')
plt.close()
