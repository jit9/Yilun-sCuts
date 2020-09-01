#!/bin/env python

"""This script aims to take a input mask based on a fits file and project
the mask into TODCuts objects and store them under a given tag"""

import numpy as np, sys, os
from enlib import enmap, config, log, pmat, mpi, utils, errors, scan as enscan
from enact import actscan, filedb, files, actdata
from moby2.tod.cuts import CutsVector, TODCuts
from moby2.util import Depot

config.default("downsample", 1, "Factor with which to downsample the TOD")
config.default("verbosity",  1, "Verbosity for output. Higher means more verbose. 0 outputs only errors etc. 1 outputs INFO-level and 2 outputs DEBUG-level messages.")

parser = config.ArgumentParser(os.environ["HOME"] + "/.enkirc")
parser.add_argument("sel")
parser.add_argument("mask")
parser.add_argument("--buffer", type=int, default=10)
parser.add_argument("--widen", type=int, default=0, help="widen the mask by the specified arcmin")
parser.add_argument("--depot", help="depot to store the output file", required=True)
parser.add_argument("--tag", help="tag under which the cuts will be stored", required=True)
parser.add_argument("prefix",nargs="?")
args = parser.parse_args()

filedb.init()
ids = filedb.scans[args.sel]

comm  = mpi.COMM_WORLD
dtype = np.float64

# Load input mask
imask = enmap.read_map(args.mask)
# Widen if necessary
if args.widen > 0:
    radius = args.widen*utils.arcmin
    imask = (~imask).distance_transform(rmax=radius)<radius
# Expand to 3 components, as the pointing code expects that
mask = enmap.zeros((3,)+imask.shape[-2:], imask.wcs, dtype)
mask[0] = imask.reshape((-1,)+imask.shape[-2:])[0]
del imask

# Setup depot
depot = moby2.util.Depot(args.depot)

# Set up logging
utils.mkdir(root + "log")
logfile   = root + "log/log%03d.txt" % comm.rank
log_level = log.verbosity2level(config.get("verbosity"))
L = log.init(level=log_level, file=logfile, rank=comm.rank)
L.info("Initialized")

# Loop through each scan
myinds = np.arange(comm.rank, len(ids), comm.size)
for ind in myinds:
    id = ids[ind]
    entry = filedb.data[id]
    try:
        d = actdata.read(entry, ["point_offsets","boresight","site","array_info"])
        d = actdata.calibrate(d, exclude=["autocut", "fftlen"])
    except (errors.DataMissing, AttributeError) as e:
        print("Skipping %s (%s)" % (id, str(e)))
        continue
    # Build a projector between samples and mask. This
    # requires us to massage d into scan form. It's getting
    # annoying that scan and data objects aren't compatible.
    bore = d.boresight.T.copy()
    bore[:,0] -= bore[0,0]
    scan = enscan.Scan(
        boresight = bore,
        offsets = np.concatenate([np.zeros(d.ndet)[:,None],d.point_offset],1),
        comps = np.concatenate([np.ones(d.ndet)[:,None],np.zeros((d.ndet,3))],1),
        mjd0 = utils.ctime2mjd(d.boresight[0,0]),
        sys = "hor", site = d.site)
    scan.hwp_phase = np.zeros([len(bore),2])
    bore_box = np.array([np.min(d.boresight,1),np.max(d.boresight,1)])
    scan.entry = d.entry
    # Is the source above the horizon? If not, it doesn't matter how close
    # it is.
    L.debug("Processing %s" % str(id))
    pmap = pmat.PmatMap(scan, mask)
    tod  = np.zeros([scan.ndet, scan.nsamp], dtype=dtype)
    pmap.forward(tod, mask)
    # cut nonzero
    tod = np.rint(tod)
    # write cuts
    dets = np.array([int(det.split('_')[-1]) for det in d.dets])
    cuts = TODCuts(det_uid=dets, nsamps=scan.nsamp, sample_offset=0)
    for i in range(len(scan.dets)):
        cuts.cuts[i] = CutsVector.from_mask(tod[i])
    # add buffer if that's what we want
    cuts.buffer(args.buffer)
    # write cuts to file
    depot.write_object(cuts, tag=args.tag, tod=entry.id)

L.info("Done")
