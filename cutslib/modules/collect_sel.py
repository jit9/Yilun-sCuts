"""This module aims to anwser the question: how much data is cut at
different pwv. To answer this, this module collects the relevant sels
for each tod in the source_scans and stores them for plotting in the
future. This will be the pre-requisite module for all modules that
look at number of dets being cut.

"""

import moby2
import os.path as op

from cutslib import Catalog
from cutslib.pathologies import get_pathologies, Pathologies
from cutslib import util

import pickle

class Module:
    def __init__(self, config):
        pass

    def run(self, p):
        all_keys = ["gainLive", "corrLive", "normLive", "rmsLive",
                    "kurtLive", "skewLive", "MFELive", "DELive"]
        # first load the catalog and narrow down to the list
        # of tods specified in the source_scans variable in the
        # cutparams file
        cpar = moby2.util.MobyDict.from_file(p.i.cutparam)
        source_scans = cpar.get("source_scans")
        # load catalog file
        catalog = Catalog()
        # narrow down to the source list
        catalog.narrow_down(source_scans)

        # next i can load in the collected results pickle file
        # but actually it doesn't contain much information except
        # the final cut, i can get a good deal of information
        # from the pathoReport object, but it doesn't contain
        # information from the bias_step information, though i
        # would be nice to have bias_step valid detectors as
        # a reference number to evaluate how well cuts are doing.
        # Therefore the best approach is to load directly the
        # pathology object associated with each TODs and parse
        # the information from there.
        base_arr = catalog.data[['tod_name','loading']].values
        n_tot = base_arr.shape[0]
        depot = moby2.util.Depot(p.depot)
        res = {
            'sel': [],
            'resp_sel': [],
            'presel': [],
            'tod_sel': [],
            'pwv': [],
            'tod': [],
        }
        for key in all_keys:
            res[key] = []
        initialized = False
        for i in range(p.rank, n_tot, p.size):
            obs = base_arr[i,0]
            pwv = base_arr[i,1]
            res['pwv'].append(pwv)
            res['tod'].append(obs)
            print(f"{p.rank:3d} {obs} {i:>5d}/{n_tot:>5d} {pwv:.2f}")
            try:
                tod = moby2.scripting.get_tod({'filename':obs,
                                               'read_data': False})
            except Exception as e:
                # data missing
                res['tod_sel'].append(False)
                continue
            if op.isfile(depot.get_full_path(Pathologies, tod=tod, tag=p.tag)):
                # load all relevant patholog results
                patho = get_pathologies({'depot': p.depot,
                                         'tag': p.tag}, tod=tod)
                # for the first data let's load the calibration data
                if not initialized:
                    ff_sel = patho.calData['ffSel']
                    ff_stable = patho.calData['stable']
                    live_cand = patho.liveCandidates
                    initialized = True
                resp_sel = patho.calData['respSel']
                # get final cuts
                patho.makeNewSelections()
                sel = patho.liveSel
                presel = patho.preLiveSel
                # store all cuts
                res['sel'].append(sel)
                res['resp_sel'].append(resp_sel)
                res['presel'].append(presel)
                for k in all_keys:
                    res[k].append(patho.crit[k]['sel'])
                res['tod_sel'].append(True)
            else:
                res['tod_sel'].append(False)

        # gather results from individual processes
        res['sel'] = util.allgatherv(res['sel'],p.comm)
        res['resp_sel'] = util.allgatherv(res['resp_sel'],p.comm)
        res['presel'] = util.allgatherv(res['presel'],p.comm)
        res['tod'] = util.allgatherv(res['tod'],p.comm)
        res['pwv'] = util.allgatherv(res['pwv'],p.comm)
        res['tod_sel'] = util.allgatherv(res['tod_sel'],p.comm)
        for k in all_keys:
            res[k] = util.allgatherv(res[k], p.comm)
        if p.rank == 0:
            # store common results
            res['ff_sel'] = ff_sel
            res['ff_stable'] = ff_stable
            res['live_cand'] = live_cand
            array_data = moby2.scripting.get_array_data({
                'instrument': 'actpol',
                'array_name': p.i.ar,
                'season': p.i.season
            })
            tes_sel = (array_data['nom_freq']== p.i.freq) * \
                (array_data['det_type'] == 'tes')
            res['tes_sel'] = tes_sel
            # write out
            outfile = op.join(p.o.patho.root, 'sel.pickle')
            with open(outfile, "wb") as f:
                pickle.dump(res, f)
            print(f"Written to: {outfile}")
