"""This script aims to produce a histogram of the number of detectors
that has bias step calibration. I will estimate the percentage of the
number with resp against the total tes detector in the given array"""

import json, os.path as op, numpy as np
import matplotlib.pyplot as plt

import moby2
from moby2.util.database import TODList
from cutslib.pathologies_tools import get_pwv
from cutslib.pathologies import Pathologies, get_pathologies

class Module:
    def __init__(self, config):
        self.todname = config.get("tod",None)
        self.tod_list = config.get("tod_list",None)
        self.limit = config.getint("limit", None)
        self.debug = config.getboolean("debug", True)
        self.force = config.getboolean("force", False)

    def run(self, p):
        todname = self.todname
        tod_list = self.tod_list
        limit = self.limit
        debug = self.debug
        force = self.force

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
        # if file exists, one may not want to redo the whole processing again
        # unless forced to do
        outfile = op.join(p.o.cal.resp, "resp_fracs.npy")
        if force or (not op.exists(outfile)):
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
        else:
             fracs = np.loadtxt(outfile)
        # make plot
        plt.figure(figsize=(8,6))
        plt.plot(np.linspace(0,1,len(fracs)), sorted(fracs), 'k', lw=2)
        plt.xlabel("Fraction of TOD", fontsize=16)
        plt.ylabel("Fraction of Dets", fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.title("Fraction of data with valid\nbias-step measurements",
                  fontsize=14, x=0.70, y=0.1)
        outfile = op.join(p.o.cal.resp, "resp_hist.png")
        print("Saving: %s" % outfile)
        plt.savefig(outfile)
        # save data for furthur processing
        outfile = op.join(p.o.cal.resp, "resp_fracs.npy")
        print("Saving: %s" % outfile)
        np.savetxt(outfile, np.array(fracs))
