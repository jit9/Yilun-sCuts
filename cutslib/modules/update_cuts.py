"""Update cuts, just like the todloop based update_crit but
is more about manually running in an interactive node for a quick
run"""

import os.path as op, os
import moby2
from cutslib import TODList
from cutslib.pathologies_tools import reportPathologies
import traceback

class Module:
    def __init__(self, config):
        # whether we want to overwrite existing db file
        self.overwrite = config.getboolean('overwrite', False)

    def run(self, p):
        overwrite = self.overwrite
        cpar = moby2.util.MobyDict.from_file(p.i.cutparam)
        # check whether db file exists
        # if op.exists(p.i.db):
        #     if not overwrite:
        #         raise RuntimeError("db file exists!")
        #     else:  # remove it
        #         os.remove(p.i.db)
        # load tod source list
        source_scans = cpar.get("source_scans")
        obsnames = TODList.from_file(source_scans)
        # remove existing obs in the db files
        if op.exists(p.i.db):
            completed = TODList.from_file(p.i.db)
            print(f"Remove {len(completed)} completed TODs from the list")
            obsnames -= completed
        # initialize reportPathologies object to collect results
        report = reportPathologies(p.i.cutparam)
        report.depot_file += f".{p.rank}"
        for i in range(p.rank,len(obsnames),p.size):
            obs = obsnames[i]
            print(f"{p.rank:2d} {i:5d}/{len(obsnames)}: {obs}")
            # load data without actually loading it
            try: report.appendResult(obs)
            except Exception as e:
                print(type(e))
                traceback.print_exc()
                continue
        p.comm.Barrier()
