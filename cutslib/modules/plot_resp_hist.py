"""This script aims to produce a histogram of the number of detectors
that has bias step calibration. I will estimate the percentage of the
number with resp against the total tes detector in the given array"""

import json, os.path as op, numpy as np
import matplotlib.pyplot as plt

import moby2
from moby2.util.database import TODList
from cutslib.pathologies_tools import get_pwv
from cutslib.pathologies import Pathologies, get_pathologies


def init(config):
    global todname, tod_list, limit, debug
    todname = config.get("tod",None)
    tod_list = config.get("tod_list",None)
    limit = config.getint("limit", None)
    debug = config.getboolean("debug", True)

def run(p):
    global todname, tod_list, limit, debug
    # load cut parameters
    params = moby2.util.MobyDict.from_file(p.i.cutparam)
    cutParams = moby2.util.MobyDict.from_file(p.i.cutParam)

    obsnames = TODList()
    if todname:
        obsnames.append(todname)
    elif tod_list:
        obsnames = TODList.from_file(tod_list)
    else:
        obsnames = TODList.from_file(params.get("source_scans"))

    # remove unprepared tods
    depot_file = p.i.db
    if op.isfile(depot_file):
        done = TODList.from_file(depot_file)
        undone = obsnames - done
        obsnames -= undone

    if limit and (limit<len(obsnames)):
        obsnames = obsnames[:limit]

    # get the list of tes detectors. first, load array_data
    array_data = moby2.scripting.get_array_data({
        'instrument': 'actpol',
        'array_name': p.i.ar,
        'season': p.i.season
    })
    tesSel = (array_data['nom_freq']== p.i.freq) * (array_data['det_type'] == 'tes')
    depot = moby2.util.Depot(p.depot)
    fracs = []
    for i,obs in enumerate(obsnames):
        print("[%d/%d] %s" % (i,len(obsnames),obs))
        tod = moby2.scripting.get_tod({'filename':obs,
                                       'read_data': False})
        # check whether relevant files exist:
        if op.isfile(depot.get_full_path(Pathologies, tod=tod, tag=p.tag)):
            # load all relevant patholog results
            patho = get_pathologies({'depot': p.depot,
                                     'tag': p.tag}, tod=tod)
            # for the first data let's load the calibration data
            if i == 0:
                # number to compare to
                # freq + tes
                sel = tesSel
                # freq + tes + ff
                sel *= patho.calData['ffSel']
                # freq + tes + ff + stable
                sel *= patho.calData['stable']
            respSel = patho.calData['respSel']
            fracs.append(np.sum(respSel * sel)*1./np.sum(sel))
    # make plot
    plt.figure(figsize=(8,6))
    pdf = fracs / np.sum(fracs)*len(fracs)
    plt.hist(fracs, normed=True, cumulative=True, bins=100, histtype='step')
    plt.xlabel("Fraction of dets", fontsize=14)
    plt.title("cdf of percentage of the detectors with\nvalid bias step measurements",
              fontsize=14, x=0.5, y=0.85)
    plt.ylabel("Fraction of TODs", fontsize=14)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    outfile = op.join(p.o.root, "resp_hist.png")
    print("Saving: %s" % outfile)
    plt.savefig(outfile)
    # save data for furthur processing
    outfile = op.join(p.o.root, "resp_fracs.npy")
    print("Saving: %s" % outfile)
    np.savetxt(outfile, np.array(fracs))
