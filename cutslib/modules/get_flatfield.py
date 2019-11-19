"""This module wraps around generateFlatfield script in moby2 to
create a flatfield dictionary that can be used as the input to the
next iteration of cuts
"""

import os
import moby2.util
import cPickle, numpy, os, sys
np=numpy

def init(config):
    pass

def getArrayStats(gains, sels, ffpar, selTODs):
    sels = numpy.array(sels, dtype=bool)
    Ndets = gains.shape[0]
    means = np.zeros(Ndets)
    stds = np.zeros(Ndets)
    for d in range(Ndets):
        # calib = numpy.array(gains[d])
        calib = 1./numpy.array(gains[d])
        calib[numpy.isnan(calib)] = 0.0
        sel = (numpy.abs(calib) < ffpar['gainLimit'])*(numpy.abs(calib) > 1./ffpar['gainLimit'])*selTODs
        if ffpar['useSel'] and len(sels[d]) > 0: sel *= numpy.array(sels[d])
        if numpy.array(sel, dtype = int).sum() < ffpar['minSamples']:
            print "Failed minSample test %d"%d
            continue
        tmp = numpy.sort(calib[sel])
        n = len(calib[sel])
        m = tmp[int(n/2)]
        q25 = tmp[int(n/4)]
        q75 = tmp[int((3*n)/4)]
        s = 0.741*(q75-q25)
        sel = (calib < m+ffpar['sigmas']*s)*(calib > m-ffpar['sigmas']*s)*(calib > 0)
        means[d] = np.median(calib[sel])
        stds[d] = calib[sel].std()
    return means, stds


def normalizeGains(gains,sels,stable):
    for i in xrange(gains.shape[1]):
        sel = np.array(sels[:,i],bool)*stable
        if sel.sum() > 0: gains[:,i] /= np.median(gains[:,i][sel])

def run(proj):
    params = moby2.util.MobyDict.from_file(proj.i.cutparam)
    ffpar = params["ff_params"]

    print 'Loading data'
    f = open(proj.i.pickle_file)
    p = cPickle.Unpickler(f)
    data = p.load()
    f.close()

    if ffpar.get('selectedTODs') is not None:
        names = []
        f = open(proj.i.root+"/"+ffpar['selectedTODs'])
        for l in f:
            if l[0] != '#':
                names.append(l.split('\n')[0].split()[0].split('/')[-1].split(".zip")[0])
        f.close()
        selTODs = []
        names = numpy.array(names)
        for name in data['name']:
            if numpy.any(names == name): selTODs.append(True)
            else: selTODs.append(False)
        selTODs = numpy.array(selTODs)
    else:
        selTODs = numpy.ones(len(data['name']), dtype = bool)

    print "Using %d TODs for flatfield"%(selTODs.sum())

    cutParams = moby2.util.MobyDict.from_file(proj.i.cutParam)
    old_ff_file = cutParams["pathologyParams"]["calibration"]["flatfield"]

    Ndets = data['gainLive'].shape[0]
    ff_old = moby2.detectors.FlatField.from_dict(old_ff_file)
    cal0 = ff_old.get_property("cal",det_uid=range(Ndets),default=1.)[1]
    calRMS0 = ff_old.get_property("cal",det_uid=range(Ndets),default=0.)[1]
    stable0 = ff_old.get_property("stable",det_uid=range(Ndets),default=False)[1]

    gains = data["gainLive"].copy()
    sel = np.asarray(data['sel'],dtype=bool)*np.asarray(data['respSel'],dtype=bool)
    if ffpar.get("normalize",True):
        normalizeGains(gains, sel, stable0)

    for it in xrange(10):
        m,s = getArrayStats(gains, sel, ffpar, selTODs)
        mm = m.copy()
        mm[m==0]=1
        st = (s/mm < 0.2)*(m > 0)
        normalizeGains(gains, sel, st)

    # Create new flatfield object
    ff = moby2.util.MobyDict.fromkeys(["det_uid","stable","cal","calRMS"])
    ff["det_uid"] = []
    ff["stable"] = []
    ff["cal"] = []
    ff["calRMS"] = []


    # Compute new flatfield values
    ld = np.where(m != 0)[0]
    ff["det_uid"] = list(ld)
    if ffpar.get("updateStable",True): ff["stable"] = list(st[ld])
    else: ff["stable"] = list(stable0[ld])
    if ffpar.has_key("maxRMS"):
        good = s[ld]/m[ld] < ffpar["maxRMS"]
        if ffpar.has_key("maxCal"):
            good *= cal0[ld]*m[ld] < ffpar["maxCal"]
        cal = cal0[ld]
        calRMS = calRMS0[ld]
        cal[good] = m[ld][good]#*cal0[ld][good]
        calRMS[good] = s[ld][good]#*cal0[ld][good]
        ff["cal"] = list(cal)
        ff["calRMS"] = list(calRMS)
    else:
        ff["cal"] = list(m[ld])#*cal0[ld])
        ff["calRMS"] = list(s[ld])#*cal0[ld])
    ff_out_file = proj.i.root + "/ff_%s.dict"%params["tag_patho"]
    ff.write_to_file(ff_out_file)

