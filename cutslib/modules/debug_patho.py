"""This script aims to act as a debugger or boilerplate for all
pathological studies of a list of TODs
website to manipulate and visualize"""

import moby2
import json, os.path as op, numpy as np
from moby2.util.database import TODList
from moby2.analysis.tod_ana import pathologies
from cutslib.pathologies_tools import get_pwv
from cutslib.pathologies import Pathologies
Pathologies = pathologies.Pathologies


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

    for obs in obsnames:
        depot = moby2.util.Depot(p.depot)
        tod = moby2.scripting.get_tod({'filename':obs,
                                       'read_data': False})
        # check whether relevant files exist:
        if op.isfile(depot.get_full_path(Pathologies, tod=tod, tag=p.tag)) and \
           op.isfile(depot.get_full_path(moby2.TODCuts, tod=tod, tag=p.tag)) and \
           op.isfile(depot.get_full_path(moby2.Calibration, tod=tod, tag=params["tag_cal"])):

            # load all relevant patholog results
            patho = moby2.scripting.get_pathologies({'depot': p.depot,
                                                     'tag': p.tag}, tod=tod)
            cuts = depot.read_object(moby2.TODCuts, tod=tod, tag=p.tag)
            calo = depot.read_object(moby2.Calibration, tod=tod, tag=params["tag_cal"])

            if debug:
                import ipdb;ipdb.set_trace()

            ###############
            # Debug below #
            ###############
