"""Grab some of the codes written for SO"""
import pickle, numpy as np
import scipy, scipy.stats as stat
import moby2
from moby2.tod.cuts import CutsVector
from cutslib import TODCuts, util
from cutslib.glitch import SnippetInfo


def analyze_scan(tod, qlim=1, vlim=0.01, n_smooth=0):
    """Gather some information about the scan

    Parameters:
    -----------
    qlim: percentile of az to find turnaround (0-100)
    vlim: lower-limit of scan speed to consider as scanning
    n_smooth: window size to smooth scan speed

    Returns:
    --------
    scan_params (dict)

    """
    # first: find useful scan parameters
    scan_params = {}
    # get turnaround
    az = np.unwrap(tod.az)
    lo, hi = np.percentile(az, [qlim,100-qlim*1])
    scan_params['az_lo'] = lo
    scan_params['az_hi'] = hi
    # get scan speed
    t = tod.ctime
    daz, dt = np.diff(az), np.diff(t)
    vaz = daz / dt
    vaz = np.r_[vaz[0], vaz]  # add the missing point
    # smooth if that's what we want
    if n_smooth > 0:
        # the first/last n_smooth/2 points are no longer trustable
        kernel = np.ones(n_smooth)/n_smooth
        vaz = np.convolve(kernel, mode='same')
    v_typ, dt_typ = np.median(np.abs(vaz)), np.median(dt)
    scan_params['vscan'] = v_typ
    scan_params['dt'] = dt_typ
    scan_params['srate'] = 1/dt_typ
    # get quick est of scan freq and period
    # (2 times faster than np.fft.rfft, off by 0.1%)
    pivots_l = np.array(CutsVector.from_mask(az < lo).mean(axis=1), dtype=int)
    pivots_h = np.array(CutsVector.from_mask(az > hi).mean(axis=1), dtype=int)
    tscan = np.median(np.diff(pivots_l)*dt_typ)
    scan_params['tscan'] = tscan
    scan_params['fscan'] = 1/tscan
    scan_params['pivots_l'] = pivots_l
    scan_params['pivots_h'] = pivots_h
    scan_params['n_scan'] = len(pivots_l)
    # estimate az min/max
    scan_params['az_min'] = np.median(az[az<lo])
    scan_params['az_max'] = np.median(az[az>hi])
    # get some flags
    scan_f = (az > lo) * (az < hi)
    turn_f = ~scan_f
    stop_f = np.abs(vaz) < v_typ * vlim
    pick_f = np.abs(vaz) > v_typ * 2
    # generate flag objects and return flags manager
    scan_params['scan_flags'] = (stop_f + pick_f)*scan_f
    scan_params['turn_flags'] = turn_f
    return scan_params


def analyze_common_mode(fdata, nsamps=1, preselector=None, pman=None):
    """perform a simple common mode analysis in a given frequency range.
    It requires the fft signal (fsignal) to be available in the tod.

    Parameters
    ----------
    tod: moby2 tod object
    nsamps: use to normalize, default to 1
    preselector: mask or a function to generate preselection mask

    Returns
    -------
    dict with common mode stats

    """
    # find correlation matrix from fourior modes
    c = fdata @ fdata.T.conj()
    a = np.linalg.norm(fdata,axis=1)
    aa = np.outer(a,a)
    aa[aa==0] = 1
    cc = c/aa
    # assuming the majority of dets are good and couple to the atmosphere,
    # preselect this set of good dets to estimate common mode.
    presel = preselect_dets(cc, preselector)
    # now we proceed with common mode analysis
    # get common mode using the pre-selected det list
    # svd in scipy is much faster than numpy
    u, s, v = scipy.linalg.svd(fdata[presel], full_matrices=False)
    cm = v[0]
    # get gains for all data through cm (not limited to pre-selected data)
    gain = np.abs(fdata @ cm.conj())/s[0]
    # get norm
    fnorm = np.sqrt(np.abs(np.diag(c)))
    norm = fnorm*np.sqrt(2./nsamps)
    # get correlations
    # note that the s[0] here is from the pre-selected data which
    # might be different to the actual s[0] using un-preselected data
    with util.nowarn():
        corr = gain * s[0] / fnorm
    cm_params = {}
    cm_params['presel'] = presel
    cm_params['cc'] = cc
    cm_params['cm'] = cm / s[0]
    cm_params['s'] = s
    cm_params['gain'] = gain
    cm_params['norm'] = norm
    cm_params['corr'] = corr
    if pman:
        for k, v in cm_params.items(): pman.add(k, v)
    return cm_params

def analyze_detector_noise(fdata, preselector=None, n_deproject=0,
                           nsamps=1, pman=None):
    """Perform a simple analysis of noise property of dets in a given
    frequency range. In particular, we look for noise level (rms) and
    deviation from gaussian statistics (skew, kurt). Results will be
    returnned in term of a dictionary. Also accept a preselector.

    Parasmeters:
    -----------
    tod: tod object of AxisManager class
    preselector: mask or function that generate a mask
    n_deproject: number of common modes to deproject before
      estimate noise properties

    """
    noises = {}
    # find det-to-det covariance and deproject first few common modes
    # before estimating noise properties
    c = fdata @ fdata.T.conj()
    presel = preselect_dets(c, preselector)
    # deproject `n_deproject` common modes from data
    if n_deproject > 0:
        # svd to covariance matrix gives U S^2 V*
        u, s2, v = scipy.linalg.svd(c[np.ix_(presel, presel)], full_matrices=False)
        # get first `n_deproject` common mode
        # kernel => (nmode, ndet)
        kernel = v[:n_deproject]/np.sqrt(s2[:n_deproject])[:,None]
        # common mode representation in freq-space
        # (nmode, ndet) . (ndet, nfreq) => (nmode, nfreq)
        modes = kernel @ fdata[presel]
        # deproject these modes
        coeff = modes @ fdata.T.conj()
        fdata -= coeff.T.conj() @ modes
        noises['s2'] = s2
    # compute the rms for the detectors
    # use: var(real space) * Nsamps = var(freq space)
    rms = np.sqrt(np.var(fdata, axis=1)/nsamps)
    # compute higher order statistics back in time domain
    tdata = fmodes_to_tmodes(fdata)
    skew = stat.skewtest(tdata,axis=1)
    kurt = stat.kurtosistest(tdata,axis=1)
    noises['c'] = c
    noises['rms'] = rms
    noises['kurt'] = kurt[0]
    noises['kurt_pval'] = kurt.pvalue
    noises['skew'] = skew[0]
    noises['skew_pval'] = skew.pvalue
    if pman:
        for k,v in noises.items():
            pman.add(k, v)
    return noises

def deproject_modes(fdata, n_modes=0, preselector=None, inplace=False):
    c = fdata @ fdata.T.conj()
    presel = preselect_dets(c, preselector)
    to_deproj = None
    # deproject `n_deproject` common modes from data
    if n_modes > 0:
        # svd to covariance matrix gives U S^2 V*
        u, s2, v = scipy.linalg.svd(c[np.ix_(presel,presel)], full_matrices=False)
        # get first `n_deproject` common mode
        # kernel => (nmode, ndet)
        kernel = v[:n_modes]/np.sqrt(s2[:n_modes])[:,None]
        # common mode representation in freq-space
        # (nmode, ndet) . (ndet, nfreq) => (nmode, nfreq)
        modes = kernel @ fdata[presel]
        # deproject these modes
        coeff = modes @ fdata.T.conj()
        if not inplace:
            fdata = fdata.copy()
        to_deproj = coeff.T.conj() @ modes
        fdata -= to_deproj
    return fdata, to_deproj

def fmodes_to_tmodes(fmodes, pad_left=1, pad_right=0, nsamps=1):
    """represent frequency domain modes as time domain modes. Note
    that it does not necessarily convert to physical frequencies.
    It simply represent the f-modes in a zero-padded array.
    modes to time domain modes. By default we left-pad with 1 samp
    to avoid freq=0 element.

    Parameters
    ----------
    fmodes: fourier space data with shape (ndets, nfreqs)
    pad_left/pad_right: zeros to pad on each side
    nsamps: a normalization number

    Returns
    -------
    time domain modes with shape (ndets, samps)

    """
    # irfft requires fmodes[:,-1] = real if nmodes is odd, so pad more
    if (pad_left+fmodes.shape[-1]+pad_right)%2==1:
        pad_right += 1
    with_pad = np.pad(fmodes,((0,0),(pad_left,pad_right)))
    normalization = np.sqrt(2*fmodes.shape[-1]/nsamps)
    modes = np.fft.irfft(with_pad) * normalization
    return modes

def corrmat(fmodes):
    """calculate correlation matrix from fourior modes"""
    c = fmodes @ fmodes.T.conj()
    a = np.linalg.norm(fmodes,axis=1)
    aa = np.outer(a,a)
    aa[aa==0] = 1
    return c/aa

def covmat(fmodes):
    return fmodes @ fmodes.T.conj()

def preselect_dets(corrmat, preselector):
    """Pre-select a set good of detectors to extract commmon mode,
    no actual work is done here, it's an interface to various
    preselection method defined in preselect module. Here we allow
    three types of preselector argument:

    1. boolean mask / det ids
    2. preselector function that takes in corrmat and generate mask
    3. tuple of preselector that define fall-back mechanisms

    Parameters:
    -----------
    corrmat (np.ndarray): n x n correlation / cov matrix of dets
    preselector: preselector

    Returns:
    --------
    det mask showing which dets are pre-selected

    """
    if preselector is None:
        return np.ones(corrmat.shape[0], dtype=bool)
    if isinstance(preselector, np.ndarray):
        if preselector.dtype == np.int:  # a det ids
            presel = np.zeros(fdata.shape[0]).astype('bool')
            presel[preselector] = 1
        elif preselector.dtype == np.bool:  # det mask
            presel = preselector
        else:
            raise ValueError(f"Unsupported preselector type")
    elif callable(preselector):  # a function
        presel = preselector(corrmat)
    elif isinstance(preselector, tuple):  # fallback list
        presel_success = False
        for selector in preselector:
            try:
                presel = preselector(corrmat)
                presel_success = True
                break  # if successful: skip the rest
            except PreselectionError:
                continue  # if failed, try next
        # raise error if nothing worked
        if not presel_success:
            raise PreselectionError
    else: raise ValueError(f"Unsupported preselector type")
    return presel


def analyze_calibration(tod, cutparams, cuts=None, write=False, **kwargs):
    import moby2
    from moby2 import products
    from moby2 import TODCuts
    # load parameters
    params = moby2.util.MobyDict.from_file(cutparams)
    cutParams = moby2.util.MobyDict.from_file(cutparams.replace('cutp','cutP'))
    pathop = cutParams['pathologyParams']
    # initialize depot
    depot = moby2.util.Depot(params.get("depot"))
    if cuts is None: cuts = depot.read_object(TODCuts, tag=params.get('tag_out'), tod=tod)
    name = tod.info.name
    # get calibration
    flatfield = pathop["calibration"]["flatfield"]
    resp = products.get_calibration(pathop["calibration"]["config"], tod.info)
    resp_sel = (resp.cal != 0.0)
    flatfield_object = moby2.detectors.RelCal.from_dict(flatfield)
    dets = tod.info.array_data['det_uid']
    _, ff = flatfield_object.get_property('cal', det_uid=dets, default=1.)
    _, stable = flatfield_object.get_property('stable', det_uid=dets, default=False)
    if flatfield_object.calRMS is not None:
        _, ffRMS = flatfield_object.get_property('calRMS', det_uid=dets, default=1.)
    else: ffRMS = np.zeros_like(dets)

    # Select detectors that passed the cuts and with valid responsivities
    # and which are stable. Also make sure they are nonzero
    sel = cuts.get_mask()*resp_sel*stable #*(gains != 0)
    if sel.sum() == 0:
        raise RuntimeError("Unable to calibrate")

    # generate calibration
    calib = np.zeros(dets.shape)
    freqs = np.unique(tod.info.array_data.get('nom_freq',["single"]))
    for freq in freqs:
        if freq == "single":
            sf = np.ones(pa.dets.size,dtype=bool)
        else:
            sf = tod.info.array_data['nom_freq'] == freq
        calib[sf] = resp.cal[sf]*ff[sf]
    calib *= tod.info.array_data["optical_sign"]

    # Store results
    s = cuts.get_uncut()
    calObj = moby2.Calibration(det_uid = dets[s])
    calObj.set_property(["cal", "calRMS"], [calib[s], ffRMS[s]])
    if write:
        depot.write_object(calObj,
                           tag=params.get("tag_cal"),
                           tod=tod,
                           force=True)
    return calObj

class TODWrapper:
    def __init__(self, tod):
        """This avoids storing a full tod inside the pathology
        object. How it works is that it parses the relevant
        information that downstream functions may need and
        act as a thin wrapper"""
        self.det_uid = tod.det_uid
        self.nsamps = tod.nsamps
        self.info = SnippetInfo.from_todinfo(tod.info)

class CutsManager:
    def __init__(self, tod=None):
        """Manage different cuts"""
        self.tod = TODWrapper(tod)
        self.cuts = {}
    def add(self, name, cuts):
        if name in self.cuts:
            raise ValueError("naming conflict")
        # we allow cuts to be either TODCuts or mask or detid
        # here i check which one it belongs to
        if isinstance(cuts, TODCuts) or isinstance(cuts, CutsVector):
            cuts_ = cuts
        elif isinstance(cuts, np.ndarray):
            cuts_ = TODCuts.for_tod(self.tod)
            if cuts.dtype == np.bool_:  # if a mask is given
                cuts_.set_always_cut(np.where(cuts==False)[0])
            else:  # assume it's detid
                cuts_.set_always_cut(cuts)
        else: raise ValueError("Unrecognized cuts format")
        self.cuts[name] = cuts_
        return self
    @classmethod
    def for_tod(cls, tod):
        cm = cls(tod)
        return cm
    def combine_cuts(self, fields=[], exclude=[]):
        # default to combine all cuts into 'final'
        final = TODCuts.for_tod(self.tod)
        if len(fields) == 0:
            fields = [k for k in self.cuts.keys()
                      if k not in exclude]
        if len(fields) != 0:
            for f in fields:
                fcuts = self.cuts[f]
                if isinstance(fcuts, CutsVector):
                    for d in final.det_uid:
                        final.add_cuts(d, fcuts)
                elif isinstance(fcuts, TODCuts):
                    final.merge_tod_cuts(fcuts)
                else: raise ValueError
        return final
    def merge(self, other):
        assert isinstance(other, CutsManager)
        overlap = [k for k in other.cuts if k in self.cuts]
        assert len(overlap) <= 1, f"Naming conflicts, rename {overlap}"
        if len(overlap) == 1: assert overlap[0] == 'final'  # obsolete
        for k in other.cuts:
            self.add(k, other.cuts[k])
        return self
    def move(self, name, newname):
        if name not in self.cuts: return self
        if newname is not None:
            self.cuts[newname] = self.cuts[name].copy()
        del self.cuts[name]
        return self

class PathologyManager:

    _depot_structure = '{class}/{tag}/{first_five}/{tod_name}.pickle'

    def __init__(self, tod=None):
        """Manage different cuts"""
        self.tod = TODWrapper(tod)
        self.patho = {}
        self.crits = {}
        self.dets = tod.info.array_data['det_uid']
    def add(self, name, patho, dets=None, default=0):
        if name in self.patho:
            raise ValueError("naming conflict")
        if dets is not None:
            patho_ = np.ones_like(self.dets, dtype=np.float)
            patho_ *= defaut
            patho_[np.asarray(dets)] = patho
            patho = patho_
        self.patho[name] = patho
        return self
    def drop(self, name):
        if name in self.patho:
            del self.patho[name]
        if name in self.crits:
            del self.crits[name]
        return self
    def restrict_dets(self, dets):
        """restrict to a list of detectors, this step discards
        the data associated with the other dets"""
        self.dets = dets
        return self
    @classmethod
    def for_tod(cls, tod):
        cm = cls(tod)
        return cm
    def add_crit(self, field, limits=[5, 95], method='rel'):
        """add a criteria to find cuts

        Parameters
        ----------
        field: field name to apply crit on
        limits: lengh 2 list with lower and higher limits, if no limits
          is needed put None accordingly, i.e. [2, None]
        method: can be 'rel' or 'abs'. When 'rel', limits refer to percentile
          and when 'abs', limits refer to absolute values
        """
        # TODO: support metadata
        assert len(limits) == 2
        lo, hi = limits
        self.crits.update({field: {
            'lo': lo, 'hi': hi,
            'method': method
        }})
        return self
    def clear_crit(self, fields=None):
        """remove crit for a given list of fields, defaults to all"""
        if not fields:
            fields = list(self.patho.keys())
        for f in fields:
            if f in self.crits:
                del self.crits[f]
        return self
    def apply_crit(self):
        """apply crit to get a cut"""
        crits = [crit for crit in self.crits if crit in self.patho]
        cman = CutsManager.for_tod(self.tod)
        for f in crits:
            v = self.patho[f]
            lo, hi = self.crits[f]['lo'], self.crits[f]['hi']
            method = self.crits[f]['method']
            m = np.ones(self.tod.det_uid.shape[0], dtype=bool)
            if method == 'rel':
                if lo: lo = np.percentile(v[self.dets], lo)
                if hi: hi = np.percentile(v[self.dets], hi)
            with util.nowarn():
                if lo: m *= (v >= lo)  # open
                if hi: m *= (v < hi)   # close
            # add masks to CutsManager
            cman.add(f, m)
        # merge the flags into a cuts field and keep the origin fields
        combined = cman.combine_cuts(fields=crits)
        cman.add('all_crits', combined)
        return cman
    def __getitem__(self, attr):
        if attr in self.patho:
            return self.patho[attr]
        else: raise ValueError
    def write_to_path(self, path):
        with open(path, "wb") as f:
            pickle.dump(self, f)
    @classmethod
    def read_from_path(cls, path, tod=None, params=None):
        with open(path, "rb") as f:
            patho = pickle.load(f)
        return patho
    def get_calibration(self, dets=None):
        # generate calibration
        cal = self.patho['ff'] * self.patho['resp']
        cal *= self.tod.info.array_data["optical_sign"]
        # Store results
        if dets is None: dets = self.dets
        cal_obj = moby2.Calibration(det_uid=dets)
        cal_obj.set_property("cal", cal[dets])
        return cal_obj
