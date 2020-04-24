"""Grab some of the codes written for SO"""
import scipy, numpy as np
from moby2.tod.cuts import CutsVector


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
    az = tod.az
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


def analyze_common_mode(fdata, nsamps, preselector=None):
    """perform a simple common mode analysis in a given frequency range.
    It requires the fft signal (fsignal) to be available in the tod.

    Parameters
    ----------
    tod: moby2 tod object
    frange: frequency range to look for common mode (in Hz)
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
    corr = gain * s[0] / fnorm
    cm_params = {}
    cm_params['presel'] = presel
    cm_params['cc'] = cc
    cm_params['cm'] = cm / s[0]
    cm_params['s'] = s
    cm_params['gain'] = gain
    cm_params['norm'] = norm
    cm_params['corr'] = corr
    return cm_params

def analyze_detector_noise(fdata, preselector=None, n_deproject=0):
    """Perform a simple analysis of noise property of dets in a given
    frequency range. In particular, we look for noise level (rms) and
    deviation from gaussian statistics (skew, kurt). Results will be
    returnned in term of a dictionary. Also accept a preselector.

    Parasmeters:
    -----------
    tod: tod object of AxisManager class
    frange: range of frequency to consider
    preselector: mask or function that generate a mask
    n_deproject: number of common modes to deproject before
      estimate noise properties

    """
    # find det-to-det covariance and deproject first few common modes
    # before estimating noise properties
    c = fdata @ fdata.T.conj()
    presel = preselect_dets(c, preselector)
    # deproject `n_deproject` common modes from data
    if n_deproject > 0:
        # svd to covariance matrix gives U S^2 V*
        u, s2, v = scipy.linalg.svd(c[presel][presel], full_matrices=False)
        # get first `n_deproject` common mode
        # kernel => (nmode, ndet)
        kernel = v[:n_deproject]/np.sqrt(s2[:n_deproject])[:,None]
        # common mode representation in freq-space
        # (nmode, ndet) . (ndet, nfreq) => (nmode, nfreq)
        modes = kernel @ fdata[presel]
        # deproject these modes
        coeff = modes @ fdata.T.conj()
        fdata -= coeff.T.conj() @ modes
    # compute the rms for the detectors
    # use: var(real space) * Nsamps = var(freq space)
    rms = np.sqrt(np.var(fdata, axis=1)/tod.samps.count)
    # compute higher order statistics back in time domain
    tdata = fmodes_to_tmodes(fdata)
    skew = stat.skewtest(tdata,axis=1)
    kurt = stat.kurtosistest(tdata,axis=1)
    noises = {}
    noises['c'] = c
    noises['rms'] = rms
    noises['s2'] = s2
    noises['kurt'] = kurt.statistics
    noises['kurt_pval'] = kurt.pvalue
    noises['skew'] = skew.statistics
    noises['skew_pval'] = skew.pvalue

def deproject_modes(fdata, n_modes=0, preselector=None):
    c = fdata @ fdata.T.conj()
    presel = preselect_dets(c, preselector)
    # deproject `n_deproject` common modes from data
    if n_modes > 0:
        # svd to covariance matrix gives U S^2 V*
        u, s2, v = scipy.linalg.svd(c[presel][presel], full_matrices=False)
        # get first `n_deproject` common mode
        # kernel => (nmode, ndet)
        kernel = v[:n_modes]/np.sqrt(s2[:n_modes])[:,None]
        # common mode representation in freq-space
        # (nmode, ndet) . (ndet, nfreq) => (nmode, nfreq)
        modes = kernel @ fdata[presel]
        # deproject these modes
        coeff = modes @ fdata.T.conj()
        fdata_deproj = fdata - coeff.T.conj() @ modes
    return fdata_deproj

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
    modes = np.fft.irfft(fcm) * normalization
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
