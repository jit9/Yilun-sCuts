import os.path as op, sys, numpy as np
from cutslib.pathologies_tools import pathoList
from moby2.util import MobyDict

class Module:
    def __init__(self, config):
        self.include_time = config.get("include_time", None)

    def run(proj):
        include_time = self.include_time

        params = MobyDict.from_file(proj.i.cutparam)
        datafile = str(proj.i.db)
        pl = pathoList(datafile)
        pl.addPWV2results()
        outfile = op.join(proj.o.root, "selTODs_%s.txt" % params.get("tag_out"))
        print("Saving file: %s" % outfile)
        sp = params.get("selParams",{"liveDets":{"gt":100},"PWV": {"lt":3}})
        tod_list = pl.data['todName']
        sel = pl.selectGoodTODs(sp)
        if include_time:
            tod_sel = select_tod_in_time(tod_list, include_time)
            pl.outputTODList(outfile,np.logical_and(sel, tod_sel))
        else:
            pl.outputTODList(outfile,sel)

####################
# utility function #
####################

def select_tod_in_time(tod_list, include_time):
    # load include time file
    with open(include_time, "r") as f:
        lines = f.readlines()
    t_start = np.array([int(t.split('.')[0]) for t in tod_list])
    mask = np.zeros(t_start.shape).astype(bool)

    for l in lines:
        l_start = int(l.split()[0])
        l_end = int(l.split()[1])
        mask[np.logical_and(t_start > l_start, t_start < l_end)] = True

    return mask
