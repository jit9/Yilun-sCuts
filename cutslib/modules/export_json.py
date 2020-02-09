"""This script exports the pathology res into json format for the
website to manipulate and visualize"""

import moby2
import json, os.path as op
from moby2.util.database import TODList

def init(config):
    global todname, tod_list, limit
    todname = config.get("tod",None)
    tod_list = config.get("tod_list",None)
    limit = config.getint("limit", None)


def run(p):
    global todname, tod_list, limit
    # load cut parameters
    params = moby2.util.MobyDict.from_file(p.i.cutparam)

    obsnames = []
    if todname:
        obsnames.append(todname)
    elif tod_list:
        obsnames = TODList.from_file(tod_list)
    else:
        obsnames = TODList.from_file(params.get("source_scans"))

    if limit:
        obsnames = obsnames[:limit]

    for obs in obsnames:
        parse_stats(obs, p)


def parse_stats(todname, p):
    """This function parses the useful stats from a given TOD
    Args:
        todname: name of the tod
        p: proj parameter as in run(proj)
    """
    # load tod
    tod = moby2.scripting.get_tod({'filename':todname,
                                   'read_data': False})
    # load pathology
    patho = moby2.scripting.get_pathologies({'depot': p.depot,
                                             'tag': p.tag}, tod=tod)
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

    export = {}
    export['tag'] = p.tag
    export['dimensions'] = ['det_uid','array_x','array_y','row','col','MFELive',\
                            'skewLive','corrLive','rmsLive','gainLive','DELive',\
                            'normLive','kurtLive','ff','resp','presel','pol_family',
                            'bias_line', 'optical_sign']
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
                res['pol_family'][i].decode('utf-8'),
                int(res['bias_line'][i]),
                int(res['optical_sign'][i])
            ])
    outfile = op.join(p.o.patho.viz, "%s.json" % todname)
    print("Writing: %s" % outfile)
    with open(outfile,"w") as f:
        f.write(json.dumps(export))
