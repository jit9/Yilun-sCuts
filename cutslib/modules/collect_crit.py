"""This is a wrapper to the collectSeasonCrit script in moby2. It is
wrapped in this framework for convenience"""
import os, os.path as op
import moby2
import pickle, numpy as np, sys, os
from moby2.util.database import TODList
from cutslib.pathologies_tools import get_pathologies
from cutslib import Catalog
from cutslib.pathologies import Pathologies
from cutslib import util

class Module:
    def __init__(self, config):
        self.limit = config.getint("limit", None)
        self.split = config.getint("split", 1)

    def run(self, p):
        limit = self.limit
        all_keys = ["gainLive", "corrLive", "normLive", "rmsLive",
                    "kurtLive", "skewLive", "MFELive", "DELive", "jumpLive"]
        # read parameters
        cpar = moby2.util.MobyDict.from_file(p.i.cutparam)
        depot_file = p.i.db
        source_scans = cpar.get("source_scans")
        catalog = Catalog(filename=p.catalog)
        # narrow down to the source list
        catalog.narrow_down(source_scans)
        base_arr = catalog.data[['tod_name','loading']].values
        depot = moby2.util.Depot(p.depot)
        if limit:
            base_arr = base_arr[:limit,:]
        # the list is usually huge so it's worthwhile to split into
        # sublists and process each one
        base_arr_splits = np.array_split(base_arr, self.split)
        for i_split in range(self.split):
            base_arr = base_arr_splits[i_split]
            res = {
                'tods': [],
                'sel': [],
                'psel': [],
                'scanf': [],
                'resp': [],
                'resp_sel': [],
                'cal': [],
                'ctimes': [],
                'pwv': [],
                'alt': [],
                'tod_sel': []
            }
            for key in all_keys:
                res[key] = []
                res[f"{key}_sel"] = []
            n_tot = base_arr.shape[0]
            print(f"Collecting Criteria for {n_tot} files in [{i_split+1}/{self.split}] split")
            for i in range(p.rank, n_tot, p.size):
                obs = base_arr[i,0]
                pwv = base_arr[i,1]
                print(f"{p.rank:3d} {obs} {i:>5d}/{n_tot:>5d} {pwv:.2f}")
                # try loading, if failed add it to tod_sel
                try:
                    tod = moby2.scripting.get_tod({"filename": obs,
                                                   "read_data": False})
                except:
                    # failed to process
                    res['tod_sel'].append(False)
                    continue
                if os.path.isfile(depot.get_full_path(Pathologies, tod=tod, tag=p.tag)) and \
                   os.path.isfile(depot.get_full_path(moby2.TODCuts, tod=tod, tag=p.tag)):
                    pa = get_pathologies(tod, cpar)
                    # get all crit cuts
                    # for final cuts we will get from TODCuts objects instead of here
                    # because the externally included cuts will not be stored here
                    pa.makeNewSelections()
                    # store pathologies crits
                    for k in all_keys:
                        if "values" in pa.crit[k]:
                            res[k].append(pa.crit[k]["values"])
                        if "sel" in pa.crit[k]:
                            res[f"{k}_sel"].append(pa.crit[k]["sel"])
                    # fix_tod_length(tod, pa.offsets)
                    # get preselection
                    res['psel'].append(pa.preLiveSel)
                    resp, _, _, re_sel, _, _ = pa.getpWCalibration()
                    res['tods'].append(tod.info.name)
                    res['scanf'].append(pa.scan_freq)
                    res['resp'].append(resp)
                    res['resp_sel'].append(re_sel)
                    res['cal'].append(resp*pa.calData['ff'])
                    res['ctimes'].append(tod.info.ctime)
                    res['pwv'].append(pwv)
                    res['alt'].append(np.mean(tod.alt))
                    res['tod_sel'].append(True)
                    # get final cuts
                    cuts = depot.read_object(moby2.TODCuts, tod=tod, tag=p.tag)
                    res['sel'].append(cuts.get_mask())
                else:
                    res['tod_sel'].append(False)

            p.comm.Barrier()
            data = {}
            data['name'] = util.allgatherv(res['tods'], p.comm)
            data["scan_freq"] = util.allgatherv(res['scanf'], p.comm)
            data['sel'] = util.allgatherv(res['sel'], p.comm).T
            data['psel'] = util.allgatherv(res['psel'], p.comm).T
            data["resp"] = util.allgatherv(res['resp'], p.comm).T
            data["resp_sel"] = util.allgatherv(res['resp_sel'], p.comm).T
            data["cal"] = util.allgatherv(res['cal'], p.comm).T
            data["ctime"] = util.allgatherv(res['ctimes'], p.comm)
            data["alt"] = util.allgatherv(res['alt'], p.comm)
            data["pwv"] = util.allgatherv(res['pwv'], p.comm)
            data['tod_sel'] = util.allgatherv(res['tod_sel'],p.comm)
            for k in all_keys:
                data[k] = util.allgatherv(np.array(res[k]), p.comm).T
                ksel = f'{k}_sel'
                data[ksel] = util.allgatherv(np.array(res[ksel]), p.comm).T
            if p.rank == 0:
                data["live"] = pa.liveCandidates
                data["dark"] = pa.origDark
                data["ff"] = pa.calData['ff']
                data["ff_sel"] = pa.calData['ffSel']
                data["stable"] = pa.calData['stable']
                array_data = moby2.scripting.get_array_data({
                    'instrument': 'actpol',
                    'array_name': p.i.ar,
                    'season': p.i.season
                })
                tes_sel = (array_data['nom_freq']== p.i.freq) * \
                    (array_data['det_type'] == 'tes')
                data['tes_sel'] = tes_sel
                outfile = p.o.pickle_file+f'.{i_split}'
                print("Saving data: %s" % outfile)
                with open(outfile, 'wb') as f:
                    pkl = pickle.Pickler(f,2)
                    pkl.dump(data)
            p.comm.Barrier()
        # merge all splits
        if p.rank == 0:
            data_splits = []
            print("Collecting and removing split data...")
            for i in range(self.split):
                infile = p.o.pickle_file+f".{i}"
                with open(infile, "rb") as f:
                    data = pickle.load(f)
                data_splits.append(data)
                os.remove(infile)
            data = {}
            for k in ['name','scan_freq','sel','psel','resp','resp_sel',\
                      'cal','ctime','alt','pwv','tod_sel']:
                data[k] = np.hstack([d[k] for d in data_splits])
            for k in all_keys:
                data[k] = np.hstack([d[k] for d in data_splits])
                ksel = f'{k}_sel'
                data[ksel] = np.hstack([d[ksel] for d in data_splits])
            # those that only needs one time
            for k in ['live','dark','ff','ff_sel','stable','tes_sel']:
                data[k] = data_splits[0][k]
            outfile = p.o.pickle_file
            print("Saving data: %s" % outfile)
            with open(outfile, 'wb') as f:
                pkl = pickle.Pickler(f,2)
                pkl.dump(data)
