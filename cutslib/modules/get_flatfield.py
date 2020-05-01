"""This module wraps around generateFlatfield script in moby2 to
create a flatfield dictionary that can be used as the input to the
next iteration of cuts.

Options:
--------
legacy: True if we are using loic's gain in pickle file. The main difference
  is that i no longer remove the flatfield in the gain in the pickle file,
  so in legacy mode the gain will come directly from the gainLive field in
  the pickle file but for not legacy mode the gain will be the gain in the
  pickle file multiplying with the flatfield.

"""

import os.path as op
import moby2.util
import pickle, numpy as np, os, sys


class Module:
    def __init__(self, config):
        self.legacy = config.get("legacy",False)

    def run(self, proj):
        legacy = self.legacy

        params = moby2.util.MobyDict.from_file(proj.i.cutparam)
        ffpar = params["ff_params"]

        print('Loading data')
        f = open(proj.i.pickle_file, "rb")
        p = pickle.Unpickler(f)
        data = p.load()
        f.close()

        if ffpar.get('selectedTODs') is not None:
            names = []
            f = open(op.join(proj.o.root,ffpar['selectedTODs']), "r")
            for l in f:
                if l[0] != '#':
                    names.append(l.split('\n')[0].split()[0].split('/')[-1].split(".zip")[0])
            f.close()
            selTODs = []
            names = np.array(names)
            for name in data['name']:
                if np.any(names == name): selTODs.append(True)
                else: selTODs.append(False)
            selTODs = np.array(selTODs)
        else:
            selTODs = np.ones(len(data['name']), dtype = bool)

        print("Using %d TODs for flatfield"%(selTODs.sum()))

        cutParams = moby2.util.MobyDict.from_file(proj.i.cutParam)
        old_ff_file = cutParams["pathologyParams"]["calibration"]["flatfield"]

        Ndets = data['gainLive'].shape[0]
        ff_old = moby2.detectors.FlatField.from_dict(old_ff_file)
        cal0 = ff_old.get_property("cal",det_uid=list(range(Ndets)),default=1.)[1]
        calRMS0 = ff_old.get_property("cal",det_uid=list(range(Ndets)),default=0.)[1]
        stable0 = ff_old.get_property("stable",det_uid=list(range(Ndets)),default=False)[1]

        if legacy:
            gains = data["gainLive"].copy()
        else:
            gains = data["gainLive"]*data["ff"][:,None]
        sel = np.asarray(data['sel'],dtype=bool)*np.asarray(data['resp_sel'],dtype=bool)
        if ffpar.get("normalize",True):
            normalizeGains(gains, sel, stable0)

        for it in range(10):
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
        if "maxRMS" in ffpar:
            good = s[ld]/m[ld] < ffpar["maxRMS"]
            if "maxCal" in ffpar:
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
        ff_out_file = op.join(proj.o.root, "ff_%s.dict"%params["tag_patho"])
        ff.write_to_file(ff_out_file)

#####################
# utility libraries #
#####################

def getArrayStats(gains, sels, ffpar, selTODs):
    sels = np.array(sels, dtype=bool)
    Ndets = gains.shape[0]
    means = np.zeros(Ndets)
    stds = np.zeros(Ndets)
    for d in range(Ndets):
        # notice the division here
        calib = 1./np.array(gains[d])
        calib[np.isnan(calib)] = 0.0
        sel = (np.abs(calib) < ffpar['gainLimit'])*(np.abs(calib) > 1./ffpar['gainLimit'])*selTODs
        if ffpar['useSel'] and len(sels[d]) > 0: sel *= np.array(sels[d])
        if np.array(sel, dtype = int).sum() < ffpar['minSamples']:
            # print("Failed minSample test %d"%d)
            continue
        tmp = np.sort(calib[sel])
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
    for i in range(gains.shape[1]):
        sel = np.array(sels[:,i],bool)*stable
        if sel.sum() > 0: gains[:,i] /= np.median(gains[:,i][sel])
