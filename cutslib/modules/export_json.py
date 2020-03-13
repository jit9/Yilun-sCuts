"""This script exports the pathology res into json format for the
website to manipulate and visualize"""

import moby2
import json, os.path as op, numpy as np, os
from moby2.util.database import TODList
from cutslib.pathologies_tools import get_pwv
from cutslib.pathologies import Pathologies, get_pathologies

class Module:
    def init(self, config):
        self.todname = config.get("tod",None)
        self.tod_list = config.get("tod_list",None)
        self.limit = config.getint("limit", None)
        self.debug = config.getboolean("debug", False)

    def run(self, p):
        todname = self.todname
        tod_list = self.tod_list
        limit = self.limit
        debug = self.debug

        # load cut parameters
        params = moby2.util.MobyDict.from_file(p.i.cutparam)

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

        # store metadata
        metadata = {
            'tod_list': [],
            'dimensions': ['det_uid','array_x','array_y','row','col','MFELive',\
                           'skewLive','corrLive','rmsLive','gainLive','DELive',\
                           'normLive','kurtLive','ff','resp','presel','pol_family',\
                           'bias_line', 'optical_sign', 'sel']}
        for obs in obsnames:
            res = parse_stats(obs, p)
            if res:
                metadata['tod_list'].append(res)

        # dump metadata as well
        outfile = op.join(p.o.patho.viz, "metadata.json")
        oldfile = outfile
        existing = False
        while op.isfile(oldfile):
            oldfile += ".old"
            existing = True
        if existing:
            print(f"Found existing: renaming existing file to {oldfile}")
            os.system(f"mv {outfile} {oldfile}")
        print("Writing: %s" % outfile)
        with open(outfile,"w") as f:
            f.write(json.dumps(metadata))

    def parse_stats(todname, p):
        """This function parses the useful stats from a given TOD
        Args:
            todname: name of the tod
            p: proj parameter as in run(proj)
        """
        # load tod
        depot = moby2.util.Depot(p.depot)
        tod = moby2.scripting.get_tod({'filename':todname,
                                       'read_data': False})
        # check whether relevant files exist:
        if op.isfile(depot.get_full_path(Pathologies, tod=tod, tag=p.tag)) and \
           op.isfile(depot.get_full_path(moby2.TODCuts, tod=tod, tag=p.tag)):
           # and op.isfile(depot.get_full_path(moby2.Calibration, tod=tod, tag=params["tag_cal"])):

            # load pathology and cuts
            patho = get_pathologies({'depot': p.depot,
                                     'tag': p.tag}, tod=tod)
            cuts = depot.read_object(moby2.TODCuts, tod=tod, tag=p.tag)
            lsel = np.zeros(patho.ndet)
            lsel[cuts.get_uncut()] = 1

            # empty dictionary to store relevant res
            res = {}
            # store array data
            res.update(patho.tod.info.array_data)
            # store pathology stats
            res['corrLive'] = patho.crit['corrLive']['values']
            res['rmsLive'] = patho.crit['rmsLive']['values']
            res['kurtLive'] = patho.crit['kurtLive']['values']
            res['skewLive'] = patho.crit['skewLive']['values']
            res['MFELive'] = patho.crit['MFELive']['values']
            res['normLive'] = patho.crit['normLive']['values']
            res['gainLive'] = patho.crit['gainLive']['values']
            res['jumpLive'] = patho.crit['jumpLive']['values']
            res['DELive'] = patho.crit['DELive']['values']
            res['ff'] = patho.calData['ff']
            res['resp'] = patho.calData['resp']
            res['presel'] = patho.preLiveSel
            res['sel'] = lsel
            # pwv = get_pwv([tod.info.ctime])
            if debug:
                import ipdb;ipdb.set_trace()
            # generate export json file
            export = {}
            export['tag'] = p.tag
            # export['pwv'] = float(pwv[0])
            export['source'] = []
            for i in range(len(res['det_uid'])):
                if res['nom_freq'][i] == p.i.freq:
                    export['source'].append([
                        int(res['det_uid'][i]),
                        float(res['array_x'][i]),
                        float(res['array_y'][i]),
                        int(res['row'][i]),
                        int(res['col'][i]),
                        float(res['MFELive'][i]),
                        float(res['skewLive'][i]),
                        float(res['corrLive'][i]),
                        float(res['rmsLive'][i]),
                        float(res['gainLive'][i]),
                        float(res['DELive'][i]),
                        float(res['normLive'][i]),
                        float(res['kurtLive'][i]),
                        float(res['ff'][i]),
                        float(res['resp'][i])*1e16,
                        int(res['presel'][i]),
                        res['pol_family'][i],
                        int(res['bias_line'][i]),
                        int(res['optical_sign'][i]),
                        int(res['sel'][i]),
                    ])
            outfile = op.join(p.o.patho.viz, "%s.json" % todname)
            print("Writing: %s" % outfile)
            with open(outfile,"w") as f:
                f.write(json.dumps(export))
            return todname
        return None
