"""This is a wrapper to the collectSeasonCrit script in moby2. It is
wrapped in this framework for convenience"""
import os
import moby2
from moby2.analysis.tod_ana import pathologies
import pickle, numpy,sys, os
from moby2.util.database import TODList
from moby2.scripting.pathologies_tools import fix_tod_length, get_pwv
Pathologies = pathologies.Pathologies

def init(config):
    pass

def run(p):
    # read parameters
    params = moby2.util.MobyDict.from_file(p.i.cutparam)
    depot = moby2.util.Depot(params.get('depot'))
    cutParams = moby2.util.MobyDict.from_file(p.i.cutParam)
    pathop = cutParams['pathologyParams']

    # get tod list
    obsnames = TODList.from_file(params.get("source_scans"))
    depot_file = p.i.db
    if os.path.isfile(depot_file):
        done = TODList.from_file(depot_file)
        undone = obsnames - done
        obsnames -= undone

    print("Collecting Criteria for %d files"%len(obsnames))
    tods = []
    sel = []
    scanf = []
    resp = []
    respSel = []
    cal = []
    ctimes = []
    criteria = {}
    alt = []
    if params.get("keys") is None:
        keys = ["corrLive", "DELive", "MFELive","darkRatioLive",# "ampLive",
                "rmsLive", "skewLive", "kurtLive",
                "normLive", "gainLive", "jumpLive",
                "corrDark", "DEDark", "normDark",
                "gainDark", "rmsDark", "jumpDark"]
    else:
        keys = params.get("keys")

    n = 0
    N = len(obsnames)
    for obs in obsnames:
        print("Collecting %s: %d out of %d"%(obs.split("/")[-1], n, N))
        n += 1
        try:
            tod = moby2.scripting.get_tod({"filename":obs, "read_data":False})
        except:
            print("Failed")
            continue
        if os.path.isfile(depot.get_full_path(Pathologies, tod=tod, tag=params['tag_patho'])) and \
           os.path.isfile(depot.get_full_path(moby2.TODCuts, tod=tod, tag=params["tag_out"])) and \
           os.path.isfile(depot.get_full_path(moby2.Calibration, tod=tod, tag=params["tag_cal"])):
            calo = depot.read_object(moby2.Calibration, tod=tod, tag=params["tag_cal"])
            if len(calo.cal)==0:
                print("No calibration available")
                continue
            pa = depot.read_object(Pathologies, tod=tod, tag=params['tag_patho'],
                                   options = {"tod":tod, "params":pathop},
                                   structure = params.get("structure"))
            for crit in keys:
                if "values" in pa.crit[crit]:
                    criteria.setdefault(crit,[]).append( pa.crit[crit]["values"] )
            fix_tod_length(tod, pa.offsets)
            cuts = depot.read_object(moby2.TODCuts, tod=tod, tag=params["tag_out"])
            re, ff, _, re_sel, _, stable = pa.getpWCalibration()
            lsel = numpy.zeros(pa.ndet)
            lsel[cuts.get_uncut()] = True
            sel.append(lsel)
            tods.append(tod.info.name)
            scanf.append(pa.scan_freq)
            resp.append(re)
            respSel.append(re_sel)
            cal.append(calo.get_property("cal", det_uid=tod.det_uid, default=1)[1])
            ctimes.append(tod.info.ctime)
            alt.append(numpy.mean(tod.alt))

    data = {}
    data['name'] = tods
    data['sel'] = numpy.array(sel).T
    data["live"] = pa.liveCandidates
    data["dark"] = pa.origDark
    data["scan_freq"] = scanf
    data["resp"] = numpy.array(resp).T
    data["respSel"] = numpy.array(respSel).T
    data["ff"] = ff
    data["cal"] = numpy.array(cal).T
    data["stable"] = stable
    data["ctime"] = numpy.array(ctimes)
    data["pwv"] = get_pwv(data["ctime"])
    data["alt"] = alt
    for k in list(criteria.keys()):
        data[k] = numpy.array(criteria[k]).T
    outfile = p.o.pickle_file
    print("Saving data: %s" % outfile)
    f = open(outfile, 'wb')
    p = pickle.Pickler(f,2)
    p.dump(data)
    f.close()
