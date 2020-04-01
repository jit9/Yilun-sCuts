"""This script creates a simple binned map of a given TOD.
It depends on enki module
"""
import matplotlib
matplotlib.use('agg')

import os, numpy as np, sys, pipes

from enlib import enmap, pmat, fft, bench
from enlib import config, mpi, utils, errors
from enact import filedb, actdata, actscan
from enlib import enplot
from enact import nmat_measure
from cutslib.util import mkdir

opj = os.path.join
comm = mpi.COMM_WORLD
dtype = np.float32

def get_afs(scan_id):
    freq = scan_id[-4:]
    season = str(12 + int(filedb.extractors['season'](scan_id[:-5])))
    pa = filedb.extractors['ar'](scan_id[:-5])
    return pa, freq, season

def write_settings(outdir):
    '''Store parameters for this run. '''

    if comm.rank == 0:
        config.save(opj(outdir, 'config.txt'))

        with open(opj(outdir, 'env.txt'), 'w') as f:
            for k,v in os.environ.items():
                f.write('{}: {}\n'.format(k, v))

        with open(opj(outdir, 'args.txt'), 'w') as f:
            argstring = " ".join([pipes.quote(a) for a in sys.argv[1:]])
            f.write(argstring + "\n")

def write_scans(outdir, filelist):
    '''
    Store scan indices.

    Parameters
    ----------
    outdir : str
        Absolute path to output directory.
    filelist : array-like of str
        Sequence of scan indices (XXXXXXXXXX.XXXXXXXXXX.arX:fXXX).
    '''

    if comm.rank == 0:
        with open(opj(outdir, 'ids.txt'), 'w') as f:
            for idx in filelist:
                f.write("{}\n".format(idx))

def allocate_output(area):
    '''
    Parameters
    ----------
    area : str
        Absolute path to footprint map.

    Returns
    -------
    cutmap : (1, ny, nx) enmap
    '''

    shape, wcs = enmap.read_map_geometry(area)
    omap = enmap.zeros(shape, wcs, dtype)
    div  = enmap.zeros((3,)+shape, wcs, dtype)
    return omap, div

def read_metadata(entry):
    '''
    Parameters
    ----------
    entry : filedb.data object

    Returns
    -------
    data : enlib.dataset.DataSet instance
    '''
    with bench.show("read"):
        data = actdata.read(entry)
    with bench.show("calibrate"):
        data = actdata.calibrate(data, exclude=['autocut'])
    if data.ndet == 0 or data.nsamp == 0:
        raise errors.DataMissing("No data in tod")
    return data

if __name__ == '__main__':

    parser = config.ArgumentParser(os.environ["HOME"]+"/.enkirc")
    parser.add_argument("area")
    parser.add_argument("sel")
    parser.add_argument("odir")
    parser.add_argument("--nrandom", type=int, default=100)
    args = parser.parse_args()

    filedb.init()
    ids = filedb.scans[args.sel]
    # create outdir if it doesn't exist
    outdir = args.odir
    mkdir(outdir, comm)
    write_settings(outdir)
    # choose nrandom ids to work on
    if comm.rank == 0:
        ids = np.random.choice(ids, args.nrandom, replace=False)
    else:
        ids = None
    ids = comm.bcast(ids, root=0)
    for scan_id in range(comm.rank, len(ids), comm.size):
        id = ids[scan_id]
        entry = filedb.data[id]
        print(f'{comm.rank}:{scan_id}:{entry.id}')
        # read tod data
        try:
            data = read_metadata(entry)
        except errors.DataMissing as e:
            print('rank {:3d}: {} skipped ({})'.format(
                comm.rank, scan_id, e.message))
            continue
        # prepare output
        omap, div = allocate_output(args.area)
        # get pointing matrix
        scan = actscan.ACTScan(entry)
        P = pmat.PmatMap(scan, omap)
        # estimate noise model
        with bench.show("fft"):
            ft = fft.rfft(data.tod) * data.nsamp**0.5
        with bench.show("nmat_measure"):
            noise = nmat_measure.detvecs_jon(ft, scan.srate)
        # apply noise model
        # copy tod so not to overwrite
        with bench.show("remove noise model"):
            tod = data.tod.copy()
            noise.apply(tod)
        # make simple binned map
        tod = tod.astype(dtype)
        with bench.show("map"):
            P.backward(tod, omap)
        # fix unit
        # represent white noise properties of the map
        todw = tod*0.+1
        noise.white(todw)
        with bench.show("map"):
            P.backward(todw, div)
        # intensity is good, polarization can be not accurate
        div = div[0]
        map_bin = omap / div  # only T map is good
        del tod, div, omap, todw
        # write fits files
        with bench.show("write"):
            # prepare outdir for each scan
            pa, freq, season = get_afs(id)
            outdir_scan = opj(outdir, f'pa{pa}_{freq}_s{season}', entry.id)
            mkdir(outdir_scan)
            enmap.write_map(opj(outdir_scan, f'tod_{entry.tag}.fits'),
                            map_bin)
        del map_bin
