"""This module merges the cuts for cmb, it aims to replace the script
mergeCuts4cmb in moby2 with added mpi support

"""

import moby2
from cutslib import MobyDict, TODList
import os.path as op

class Module:
    def __init__(self, config):
        pass

    def run(self, p):
        # get basic tags
        cpar = MobyDict.from_file(p.i.cutparam)
        tag_cuts = cpar.get('tag_out')
        tag_planet = cpar.get('tag_planet')
        tag_cmb = cpar.get('tag_cmb')
        # get source list
        source_scans = cpar.get('source_scans')
        obsnames = TODList.from_file(source_scans)
        # get depot
        depot_path = cpar.get('depot')
        depot = moby2.util.Depot(depot_path)

        n_tot = len(obsnames)
        for i in range(p.rank,n_tot,p.size):
            obs = obsnames[i]
            print(f"{p.rank:2d} {i:5d}/{n_tot} {obs}")
            params_loadtod = {
                'filename': obs,
                'read_data': False,
            }
            try:
                tod = moby2.scripting.get_tod(params_loadtod)
            except: continue
            cmb_path = depot.get_full_path(moby2.TODCuts, tag=tag_cmb, tod=tod)
            if not op.exists(cmb_path):
                cuts_path = depot.get_full_path(moby2.TODCuts, tag=tag_cuts, tod=tod)
                planet_path = depot.get_full_path(moby2.TODCuts, tag=tag_planet, tod=tod)
                if op.exists(cuts_path):
                    cuts = depot.read_object(moby2.TODCuts, tag = tag_cuts, tod = tod)
                    if op.exists(planet_path):
                        planet_cuts = depot.read_object(moby2.TODCuts, tag=tag_planet, tod=tod)
                        cuts.merge_tod_cuts(planet_cuts)
                    depot.write_object(cuts, tag=tag_cmb, force=True, tod=tod, make_dirs=True)
        p.comm.Barrier()
