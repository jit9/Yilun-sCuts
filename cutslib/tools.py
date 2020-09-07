from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from past.builtins import basestring

import numpy as np
from cutslib.errors import PreselectionError
# import pyfftw

def power( x, dt = 1., nbin = 1, binsize = 0, detrend = False,
           Bartlett = False, Welch = False, Hann = False,
           useRegular = False):
    """
    Take the power spectrum of input x.

    dt is the time sample spacing (seconds).
    Can break a long input up into sections or "bins" before transforming.
    Frankly, bin seems like a crap name for this, but it's in the interface,
    so I guess we are stuck with it.  Read "section" wherever you see "bin".
    nbin is the number of sections to take.
    binsize is the number of samples per section.
    Caller can specify either nbin or binsize but not both.
    When you ask for 2 sections, actually average power over the 1st 50%
    of data, the middle 50%, and the last 50%; in general, averages power
    over (2nbin-1) half-offset ranges.

    Returns (power, nu, window) where power is the power spectrum in
    units of x^2 per Hz, nu is the frequency sample vector in Hz,
    and window is the window function weight applied to the timestream.
    """

    if binsize and nbin != 1:
        raise ValueError("Please specify either binsize or nbin, but not both")

    nsamp = np.shape(x)[-1]

    if binsize == 0:
        binsize = int(nsamp/nbin)             # length of a bin
    else:
        nbin = int(nsamp/binsize)

    if nbin <= 0:
        raise ValueError("You have requested %d bins (len = %d)" % \
            (nbin, nsamp))

    if detrend: detrendData(x, window = 200)

    if Bartlett + Hann + Welch > 1:
        raise ValueError("Please choose at most one type of window")

    if Bartlett:
        window = 1 - abs((np.arange(binsize) - binsize/2.0)/(binsize/2.0))
    elif Hann:
        window = 0.5*(1-np.cos(2*math.pi*(np.arange(binsize))/binsize))
    elif Welch:
        window = 1 - pow((np.arange(binsize) - binsize/2.0)/(binsize/2.0), 2)
    else:
        window = 1.0*np.ones(binsize)


    one_d = x.ndim == 1
    if one_d:
        x.shape = (1,-1)

    if useRegular:
        nt = nextregular(binsize)
    else:
        nt = binsize

    power = 0
    if nbin != 1:
        for b in xrange(2*nbin - 1):
            y = x[:,b*binsize/2 : b*binsize/2 + binsize].copy()
            detrendData(y, window = 200)
            fx = np.fft.rfft(window[np.newaxis,:]*y,nt)
            power += (fx.real*fx.real + fx.imag*fx.imag)
            #fx = np.fft.fft(window[np.newaxis,:]*y)
            #fx = pyfftw.interfaces.numpy_fft.fft(window*y)
            #power += (fx.real*fx.real + fx.imag*fx.imag)[:,:binsize/2+1]
    else:
        y = x.copy()
        detrendData(y, window = 200)
        fx = np.fft.rfft(window[np.newaxis,:]*y,nt)
        power = (fx.real*fx.real + fx.imag*fx.imag)
        #fx = np.fft.fft(window[np.newaxis,:]*x)
        #fx = pyfftw.interfaces.numpy_fft.fft(window*x)
        #power = (fx.real*fx.real + fx.imag*fx.imag)[:,:binsize/2+1]


    # Normalizations
    power_scale = 2.0                     # 1-sided power
    power_scale /= (2.*nbin - 1.)         # Allow for multiple 'bins'
    power_scale /= pow(sum(window),2)     # Allow for window
    power_scale *= float(nt)*dt      # per unit time (total time in a bin)
    #power_scale *= float(binsize)*dt      # per unit time (total time in a bin)
    power *= power_scale

    # 1-sided power double-counts at 0 and Nyquist frequency.  Correct that.
    #power[0] /= 2.
    #power[-1] /= 2.

    nf = power.shape[-1]
    nu = np.arange(nf) / float(nf) / (2*dt)
    nu[0] = 0.5*nu[1]  # Make sure the x=0 isn't killing power

    if one_d:
        x.shape = -1
        power.shape = -1

    return power, nu, window

def detrendData(y, window = 1000):
    """
    @brief Remove the trend and mean from a data vector
    @param y        Data to detrend
    @param window   Number of elements to consider at each end of vector
    """
    n = y.shape[-1]
    one_d = y.ndim == 1
    if one_d: y.shape = (1,-1)
    if window > n/2: window = n/2
    y0 = np.mean(y[:,:window],axis=1)
    y1 = np.mean(y[:,-window:],axis=1)
    m1 = (y1+y0)/2.0
    m2 = np.mean(y,axis=1)
    slope = (y1-y0)/(n-1)
    x = np.arange(n)
    y -= (y0 - m1 + m2)[:,np.newaxis].repeat(n,1) + slope[:,np.newaxis] * x[np.newaxis,:]
    if one_d: y.shape = -1


def nextregular(n):
    while not checksize(n): n+=1
    return n

def checksize(n):
    while not (n%16): n//=16
    while not (n%13): n//=13
    while not (n%11): n//=11
    while not (n%9): n//=9
    while not (n%7): n//=7
    while not (n%5): n//=5
    while not (n%3): n//=3
    while not (n%2): n//=2
    return (1 if n == 1 else 0)

def presel_by_median(cc, sel=None, **kwargs):
    """
    minCorr: minimum correlation requiered for preselection
    superMinCorr: minimum correlation requiered in case minCorr produces less than
        max(<<minSel>>,numberOfDetectors/<<minFrac>> detectors
    minSel: minimum number of detectors preselected
    minFrac: determines the minimum number of detectors preselected by determining a
        fraction of the number of detectors available.
    Note: to go back to c9 you can set:
        superMinCorr = 0.5
        minSel = 0
        minFrac = 10000
    """
    if sel is None:
        sel = np.ones(cc.shape[0],dtype=bool)

    minCorr = kwargs.get("minCorr", 0.6)
    superMinCorr = kwargs.get("superMinCorr", 0.3)
    minSel = kwargs.get("minSel", 10)
    minFrac = kwargs.get("minFrac", 10)

    # select those detectors whose medium are above a specified threshold
    sl = (np.median(abs(cc),axis=0) > minCorr)*sel

    if kwargs.get("forceSel") is not None:
        sl *= kwargs.get("forceSel") # NOT PRETTY

    if sl.sum() < np.max([cc.shape[0]//minFrac,minSel]):
        print(f"sl.sum()={sl.sum()},minSel={minSel}")
        print("ERROR: did not find any valid detectors for low frequency analysis.")
        sl = (np.median(abs(cc),axis=0) > superMinCorr)*sel
        if sl.sum() < minSel:
            print(f"-> sl.sum()={sl.sum()},minSel={minSel}")
            raise PreselectionError("PRESELECTION FAILED, did not find any valid detectors for low frequency analysis.")
    else:
        sl = ((abs(cc[sl]).mean(axis=0)-1./len(cc[sl]))*len(cc[sl])/(len(cc[sl])-1) > minCorr)*sel
    return sl


def group_detectors(cc, sel = None, **kwargs):
    """
    Groups detectors according to their correlation.
    Returns:
        G: list of lists of detectors containing the indexes of the detectors in each group
        ind: index of the last group included in the live detector preselection
        ld: indexes of detectors from the main correlated groups
    Note: Indexes are provided according to the correlation matrix given
    """
    thr0 = kwargs.get("initCorr",0.99)
    thr1 = kwargs.get("groupCorr",0.8)
    thrg = kwargs.get("minCorr",0.6)
    dthr = kwargs.get("deltaCorr", 0.005)
    Nmin = kwargs.get("Nmin",20)
    Gmax = kwargs.get("Gmax",5)

    if sel is None: sel = np.ones(cc.shape[0],dtype=bool)
    smap = np.where(sel)[0]
    scc = cc[sel][:,sel]

    G = []
    g0 = []
    allind = np.arange(scc.shape[0])
    ss = np.zeros(scc.shape[0],dtype=bool)
    thr = thr0
    while ss.sum() < len(allind):
        if np.sum(~ss) <= Nmin or len(G) >= Gmax:
            G.append(smap[np.where(~ss)[0]])
            break

        ind = allind[~ss]
        N = np.sum(~ss)
        cco = scc[~ss][:,~ss]

        # Find reference mode
        n0 = np.sum(np.abs(cco)>thr,axis=0)
        imax = np.argmax(n0)
        if n0[imax] < np.min([Nmin,N//2]):
            thr -= dthr
            continue

        # Find initial set of strongly correlated modes
        gg = np.where(np.abs(cco[imax])>thr)[0]
        s = np.argsort(cco[imax][gg])
        g = ind[gg[s]].tolist()

        # Extend set until thr1
        while thr > thr1:
            thr -= dthr
            sg = np.ones(scc.shape[0],dtype=bool)
            sg[g] = False
            sg[ss] = False
            ind = np.where(sg)[0]

            if np.sum(sg) <= Nmin:
                g0.extend(np.where(sg)[0])
                break

            cci = scc[g][:,sg]
            g0 = np.where(~np.any(np.abs(cci)<thr,axis=0))[0].tolist()
            g.extend(ind[g0])

        # Append new group result
        G.append(smap[g])
        ss[g] = True
        #print(len(g), thr)
        thr = thr0

    ind = 0
    ld = G[ind].tolist()
    while (ind < len(G)-1):
        if np.mean(np.abs(cc[G[0]][:,G[ind+1]])) > thrg:
            ind += 1
            ld.extend(G[ind].tolist())
        else:
            break

    return G, ind, ld, smap

def get_sine2_taper(frange, edge_factor = 6):
    # Generate a frequency space taper to reduce ringing in lowFreqAnal
    band = frange[1]-frange[0]
    edge = band//edge_factor
    x = np.arange(edge, dtype=float) / (edge-1)
    taper = np.ones(band)
    taper[:edge] = np.sin(x*np.pi/2)**2
    taper[-edge:] = np.sin((x+1)*np.pi/2)**2
    return taper

def get_iharm(frange, df, scan_freq, wide = False):
    # Get the harmonic mode of the scan frequency
    n_harm = int(np.ceil(frange[1]*df/scan_freq))
    f_harm = (np.arange(n_harm)+1)*scan_freq
    if wide:
        i_harm = np.array(np.sort(np.hstack(
                      [np.round(f_harm/df-1),
                       np.round(f_harm/df),
                       np.round(f_harm/df+1)])),dtype=int)
    else:
        i_harm = np.array(np.round(f_harm/df), dtype=int)
    i_harm = i_harm[(i_harm>=frange[0])*(i_harm<frange[1])] - frange[0]
    return i_harm


def get_time_domain_modes(fmodes, n_l, nsamps, df=1.,):
    """What's happening here is: supposed fmodes has shape (ndet, nfreq)
              n_l      nfreq-1       1
    fcm = [0 0 ... 0] [        ] [        ]
          [0 0 ... 0] [        ] [ real   ]
     ndet [. . ... 0] [ fmodes ] [ fmodes ]
          [. . ... 0] [        ] [  -1    ]
          [0 0 ... 0] [        ] [        ]

    """
    # Get modes in time domain
    if fmodes.ndim == 1:
        fmodes = fmodes[np.newaxis,:]
    fcm = np.hstack([np.zeros((len(fmodes),n_l)),
                     fmodes[:,:-1],
                     np.expand_dims(np.real(fmodes[:,-1]),1)])
    modes = np.fft.irfft(fcm)
    modes_dt = 1./modes.shape[1]/df
    modes *= np.sqrt(2.*fmodes.shape[1]/nsamps)
    return modes, modes_dt


def printDictionary( dictio, tabLevel ):
    """
    @brief Visualize dictionary contents
    """
    for k in dictio.keys():
        if isinstance(dictio[k], dict):
            print("%s%s:"%(tabLevel*4*' ', k))
            printDictionary(dictio[k], tabLevel+1)
        else:
            print("%s%s:"%(tabLevel*4*' ', k), dictio[k])


def compareDictionaries( testDict, refDict):
    """
    @brief Compare two dictionaries to make sure that they have the same entries
    """
    same = True
    for k in refDict.keys():
        same *= np.any(np.array(list(testDict.keys())) == k)
        if not(same): return False
        elif isinstance(refDict[k], dict):
            same *= compareDictionaries(testDict[k], refDict[k])
    return same

def moments(data):
    """
    @brief Find the variance, skewness and kurtosis of the data
    @param data  Data array with vectors to analize
    @return var, skew, kurt  Variance, Skewness and Kurtosis
    """
    n = len(data[0])
    dat = data*data
    m2 = np.sum(dat, axis = 1)/n
    dat = dat*data
    m3 = np.sum(dat, axis = 1)/n
    dat = dat*data
    m4 = np.sum(dat, axis = 1)/n
    del dat
    var = m2
    skew = m3/m2**(3./2.)
    kurt = m4/m2**2

    return var, skew, kurt


def skewTest(skew, n):
    """
    @brief Transforms the skewness to a normal distribution
    """
    n = float(n)
    Y = skew*np.sqrt((n+1.)*(n+3.)/6./(n-2.))
    b = 3.*(n**2+27.*n-70.)*(n+1.)*(n+3.)/(n-2.)/(n+5.)/(n+7.)/(n+9.)
    w2 = np.sqrt(2.*(b-1.)) - 1.
    delta = 1./np.sqrt(np.log(w2)/2.)
    alfa = np.sqrt(2./(w2-1.))
    return delta*np.log(Y/alfa + np.sqrt((Y/alfa)**2+1.))

def kurtosisTest(kurt, n):
    """
    @brief Transforms the kurtosis to a normal distribution
    """
    n = float(n)
    e = 3.*(n-1.)/(n+1.)
    v = 24.*n*(n-2.)*(n-3.)/(n+1.)**2/(n+3.)/(n+5.)
    mask = kurt != 0
    x = (kurt[mask] - e)/np.sqrt(v)
    b = 6.*(n**2-5.*n+2.)/(n+7.)/(n+9.)*np.sqrt(6.*(n+3.)*(n+5.)/n/(n-2.)/(n-3.))
    A = 6.+8./b*(2./b+np.sqrt(1.+4./b**2))
    kt = np.zeros(kurt.shape)
    kt[mask] = ((1.-2./9./A)-((1.-2./A)/(1.+x*np.sqrt(2./(A-4.))))**(1./3.))/np.sqrt(2./9./A)
    kt[~mask] = -1000
    return kt

# This must be called after removing the common mode
def findNoiseLevel(data, nwin = 10, winsize = 1000):
    """
    @brief  Find the RMS noise of the vectors in data. It samples the vectors in order
            to estimate the white noise level.
    @param nwin      number of windows to use.
    @param winsize   number of elements in each window.
    """
    ndet = np.size(data,0)
    ndata = np.size(data,1)
    step = int(ndata/nwin)
    if nwin*winsize > ndata:
        winsize = step
    noises = np.zeros(ndet)
    for i in range(ndet):
        ns = []
        for j in range(nwin):
            ns.append(np.var(data[i][j*step:j*step+winsize]))
        noises[i] = np.median(ns)
    return noises


class cutsSummary(object):
    """
    @brief Class to visualize cut criteria selections
    """
    def __init__(self, pa, type = "live"):
        if type == "live": self.keys = pa.activeLiveKeys
        else: self.keys = pa.activeDarkKeys
        self.sel = []
        for k in self.keys:
            if pa.crit[k]["apply"]:
                self.sel.append(pa.crit[k]["sel"])
        self.sel = np.array(self.sel,dtype=int)
    def plot(self, sel = None):
        import pylab
        if sel is None: sel = np.ones(self.sel.shape[1],dtype=bool)
        pylab.figure(333,figsize=[15,5])
        m=pylab.matshow(self.sel[:,sel],fignum=0)
        m.axes.set_aspect(30)
        m.axes.set_yticks(list(range(len(self.keys))))
        m.axes.set_yticklabels(self.keys)
        pylab.draw()


# function to find glitches
def find_glitches(tod, gp):
    """find glitch cuts of a tod, the different between this
    and the find glitch cut function in moby2 is that it can
    work with MCE file more easily"""
    import moby2.libactpol as libactpol
    # define filter function
    filtvec = moby2.tod.filters.sine2highPass(tod, fc=gp['highPassFc']) * \
        moby2.tod.filters.gaussianFilter(tod, timeSigma=gp['tGlitch'])
    # get glitch cuts
    glitch_cuts = libactpol.get_glitch_cuts(
        tod.data,
        np.asarray(tod.det_uid, dtype='int32'),
        np.asarray(filtvec, dtype='float32'),
        gp['nSig'], gp['maxGlitch'], gp['minSeparation'])
    return glitch_cuts

# function to load mce data
def load_mce_data(infile, dt=1):
    """Wrap the MCeFile with a thin layer of TOD structure
    so some of the tools for TOD works on MCEFile as well
    """
    mce = MCEFile(filename=infile)
    tod = mce.Read()
    tod.nsamps = tod.data.shape[-1]
    tod.det_uid = np.arange(tod.data.shape[0])
    # give a pseudo time
    tod.ctime = np.arange(tod.nsamps,dtype=np.float32)*dt
    tod.data = tod.data.astype(np.float32)
    return tod

################
# cuts related #
################

# find cuts that affect multiple detectors
def find_common_cuts(cuts, min_dets=5, nsamps=None):
    """
    Args:
        cuts: list of cuts for each det
        min_dets: minimum dets required for common cuts
        nsamps (int): number of samples, default to max in cuts
        """
    # find the maximum sample in cuts if nsamps is not speficied
    if not nsamps:
        cuts_ends = [c[-1][-1] for c in cuts if len(c)!=0]
        if len(cuts_ends) == 0:
            return [],[]
        nsamps = np.max([c[-1][-1] for c in cuts if len(c)!=0])

    # make a mask
    mask = np.zeros(nsamps, dtype=int)
    for i, c in enumerate(cuts):
        for dc in c:
            mask[dc[0]:dc[1]] += 1
    # find common cuts enforce min dets
    t_ = np.where(mask > min_dets)[0]
    ranges = find_ranges_from_list(t_)
    dets = [find_dets_with_cut(cuts, r) for r in ranges]
    return dets, ranges

def find_ranges_from_list(lst, min_sep=1):
    """Find ranges in the list, for example, the input
    list can be [1,2,3,4,8,9,10]. This function will return
    [[1,4],[8,10]].
    """
    ranges = []
    for i in range(len(lst)):
        if i == 0:
            start_ = lst[i]
        else:
            new_ = lst[i]
            if (new_ - prev_) > min_sep:
                end_ = prev_
                ranges.append([start_,end_])
                start_ = new_
            elif i == len(lst)-1:
                end_ = new_
                ranges.append([start_,end_])
        prev_ = lst[i]
    return ranges


def find_dets_with_cut(cuts, target_cut):
    """Find the dets that has certain target_cuts
    Args:
        cuts: cuts for all dets
        target_cut: a cut to find among dets"""
    dets = []
    for i,c in enumerate(cuts):
        if len(c) > 0:
            for dc in c:
                start_ = max([dc[0],target_cut[0]])
                end_ = min([dc[-1],target_cut[-1]])
                if end_ - start_ > 0:
                    dets.append(i)
                    break
    return dets


################
# cuts and cal #
################

# calibration and sel tools
def quick_cuts_and_cal(data, dt=1):
    """Quick tool to find the cuts and cal for tod-like data
    Args:
        data: input tod-like data (ndets, nsamps)
        dt: unit time
    Rets:
        cut, cal
    """
    ndets = data.shape[0]
    nsamps = data.shape[1]
    # remove mean and detrend
    data_ = data.copy()
    data_ -= np.mean(data_,axis=1)[:,None]
    tools.detrendData(data_)
    # fft
    nu, fdata = cuts_rfft(data, dt)
    # low freq
    freq_sel = (nu > 0.01) * (nu < 0.05)
    lf_data = fdata[:,freq_sel]
    # get correlation matrix
    cc = corr_mat(lf_data)
    # preselection based on median
    # presel = np.ones(ndets)
    presel = tools.presel_by_median(cc)
    # find common mode
    u, s, v = np.linalg.svd(lf_data, full_matrices=False)
    # gain
    gain = u[:,0]
    # gain_sel = np.ones(ndets)
    gain_sel = (np.real(gain) > 0.01)
    # high freq
    freq_sel = (nu > 10) * (nu < 20)
    hf_data = fdata[:,freq_sel]
    rms = np.sqrt(np.sum(abs(hf_data)**2,axis=1)/hf_data.shape[1]/nsamps)
    rms_sel = np.ones(ndets)
    # rms_sel = rms < 8e3
    # return np.array(presel*gain_sel*rms_sel,dtype=bool), 1./gain
    return np.array(presel*gain_sel*rms_sel,dtype=bool), 1./rms

def cuts_rfft(data, dt=1, use_regular=True):
    """Quick tool to get rfft for tod-like data"""
    if use_regular:
        nf = tools.nextregular(data.shape[-1])
    else:
        nf = data.shape[-1]
    fdata = np.fft.rfft(data, nf)
    nf = fdata.shape[-1]
    nu = np.arange(nf) / nf / (2*dt)
    return nu, fdata

# Get correlation matrix
def corr_mat(fdata):
    """calculate correlation matrix for tod-like data
    in frequency space"""
    c = np.dot(fdata,fdata.T.conjugate())
    a = np.linalg.norm(fdata,axis=1)
    aa = np.outer(a,a)
    aa[aa==0.] = 1.
    cc = c/aa
    return cc
