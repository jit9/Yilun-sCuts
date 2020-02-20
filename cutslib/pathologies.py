from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from past.builtins import basestring

import numpy as np, time, pickle, os, matplotlib
import scipy.stats as stat

pylab = None
def import_pylab(interactive = True):
    global pylab
    import pylab as pl
    pylab = pl
    if interactive: pylab.ion()
    else: pylab.ioff()

import moby2
import moby2.util.log as psLib
import moby2.analysis.tod_ana.correlations
import scipy.stats.mstats as ms
from moby2.scripting import products
from scipy.cluster.vq import kmeans2

from cutslib.tools import *


class Pathologies( object ):
    """
    @brief The pathologies object calculates a series os statistics about a TOD to
    determine a set of functioning live detectors and a set of dark and dead detectors.
    It also reports plots with the various statistics, as well as statistical indicators.
    """
    # Depot support
    _depot_structure = '{class}/{tag}/{first_five}/{tod_name}.pickle'


    def __init__( self, tod, params, calibrated = False, calibratedTOD = False, noExclude = False):
        """
        @brief Initialize the pathologies object.
        @param   tod          TOD object to analyze.
        @param   calibrated   set True if the data has been previously calibrated
                              BEFORE FINDING pathology statistics.
        @param   noExclude    if set, the previously cut detectors are not excluded.
        @param   flatfield      address to calibration dictionary (flat-field)
        @param   params       Pathology parameter dictionary defined as follows:
            findPathoParams: (define what needs to be calculated)
                useCorr:      Find correlation with common mode
                useNoiseFit:  Fit power law model to data to obtain knees and noise
                useNormality: Find RMS, skewness and kurtosis in chunks
                cutJumps:     Find small jumps
                jumpParams:   Small jump parameters
                    threshold1:   Signal/Noise needed to declare a detection
                    window1:      Window size for mean calculation
                    window2:      Window size for median calculation
                    nSamples:     Number of samples for median comparison
                    skip:         Number of times to skip "window2" between samples
            makeSelParamsLive: (decide how to make (generate) live detector selections)
                selType: Select between relative or absolute selections. The options
                         are 'relative', 'absolute', 'or' or 'and'. Live detectors only.
                         The dictionary fields are:
                    correlation:  selection type for correlation
                    driftError:   selection type for drift error
                    norm:         selection type for norm
                    gain:         selection type for gain
                    knee:         selection type for knee (kneeHardLimit)
                    rms:          selection type for rms (noise from chunks)
                    skewness:     selection type for skewness
                    kurtosis:     selection type for kurtosis
                relSigma: (Criteria for relative selections)
                    sigmaCorr:   dispertion of correlation with "live common mode"
                    sigmaDE:     dispertion of "driftError".
                    sigmaNorm:   dispertion of "data vector norm".
                    sigmaGain:   dispertion of "gain".
                    sigmaKnee:   dispertion of "1/f knee".
                    sigmaRMS:    dispertion of "white noise rms" from chunks
                    sigmaSKEW:   dispertion of "skewness" for from chunks.
                    sigmaKURT:   dispertion of "kurtosis" for from chunks.
                absCrit: (Criteria for absolute selections)

                    corrLimit:   lower limit for correlation with common mode
                    DELimit:     upper limit for drift error
                    normLimit:   upper limit for data norm
                    gainLimit:   maximum gain factor versus the common mode
                    kneeLimit:   upper limit for knee (kneeHardLimit)
                    rmsLimit:    [lower, upper] limits of RMS noise from chunks
                    skewLimit:   [lower, upper] limits of Skewness from chunks
                    kurtLimit:   [lower, upper] limits of Kurtosis from chunks
                gainCrit: Maximum allowed dispersion of gains with respect to the
                          common modes before considering that the calibration is
                          unreliable.
            liveSelParams: (decide how to combine live selections)
                correlation:   Use live common mode correlations
                driftError:    Use drift error selection
                norm:          Use norm of data vector
                gain:          Use gain versus live common mode
                knees:         Use 1/f knees
                rms:           Use noise rms from chunks
                skew:          Use skewness from chunks
                kurtosis:      Use kurtosis from chunks
                partialRMS:    Cut detectors with more that certain fraction of RMS chunks cut
                partialSKEW:   Same for SKEW chunks
                partialKURT:   Same for KURT chunks
                jump:          Exclude detectors that suffered jumps
                forceDark:     Force live detectors to NOT belong with the original set
                               of dark detectors.
            makeSelParamsDark: (decide how to make (generate) dark detector selections)
                selType: Select between relative or absolute selections. The options
                         are 'relative', 'absolute', 'or' or 'and'. Live detectors only.
                         The dictionary fields are:
                    correlation:  selection type for correlation
                    driftError:   selection type for drift error
                    norm:         selection type for norm
                    gain:         selection type for gain versus dark common mode
                    rms:          selection type for rms (noise from chunks)
                relSigma: (Criteria for relative selections)
                    sigmaCorr:   dispertion of correlation with "dark common mode"
                    sigmaDE:     dispertion of "driftError".
                    sigmaNorm:   dispertion of "data vector norm".
                    sigmaGain:   dispertion of "gain".
                    sigmaRMS:    dispertion of "white noise rms" from chunks
                absCrit: (Criteria for absolute selections)

                    corrLimit:   lower limit for correlation with common mode
                    DELimit:     upper limit for drift error
                    normLimit:   upper limit for data norm
                    gainLimit:   maximum gain factor versus the common mode
                    rmsLimit:    [lower, upper] limits of RMS noise from chunks
            darkSelParams: (decide how to combine dark selections)
                correlation:   Use common mode correlations
                driftError:    Use drift error selection
                norm:          Use norm of data vector
                gain:          Use gain with respect to dark common mode
                rms:           Use white noise rms
                forceDark:     Force live detectors to NOT belong with the original set
                               of dark detectors.
            partialCuts: (parameters related to partial cuts from chunks)

                RMS:   wether to perform partial cuts based on RMS chunks (True or False)
                SKEW:  wether to perform partial cuts based on SKEW chunks (True or False)
                KURT:  wether to perform partial cuts based on KURT chunks (True or False)
                maxFracCut: maximum fraction of chunks cut in detector before cutting it.
        """

        self.tod = tod
        self.name = tod.info.name
        self.dets = tod.info.det_uid.copy()
        self.ndet = len(self.dets)
        self.ndata = tod.nsamps
        self.Nrows = np.unique(tod.info.array_data['row']).size
        self.Ncols = np.unique(tod.info.array_data['col']).size
        self.rows, self.cols = tod.info.array_data.get_property(
            ['row','col'], det_uid=self.dets)
        self.sampleTime = (tod.ctime[-1] - tod.ctime[0]) / (tod.ctime.shape[0]-1)
        dsl = tod.info.downsample_level
        self.offsets = (tod.info.sample_index*dsl, tod.nsamps*dsl)
        self.downsample4rms = dsl

        self.calibrated = calibrated
        self.calibratedTOD = calibratedTOD
        self.calData = None

        # Read Parameters
        self.params = params
        self.setParams(noExclude)

        self.crit = {}
        self.initializeCriteria()

        # Flags
        self.gainCut = False
        self.temperatureCut = False
        self.removeTemp = False
        self.azErrorCode = None
        self.Temp = -1
        self.dTemp = -1

        # Az analysis results
        self.scan = None

        # STATISTICS
        self.cm = None
        self.cmd = None
        self.noises = None
        self.powlaw = None
        self.totalRMS = None
        self.cmodes = None
        self.fcmodes = None

        # SELECTIONS
        self.liveSel = None
        self.darkSel = None
        self.zeroSel = None
        self.extraSel = {}

        self.chunkParams = None
        self.report = None

    def setParams( self, noExclude = False):
        """
        Update parameters in pathologies object
        """
        # Get lists of detectors
        exclude, dark, liveCandidates = get_detector_params(self.params["detectorLists"])
        if 'det_uid' in exclude:
            self.exclude = self.tod.info.array_data.select_inner(
                {'det_uid': exclude['det_uid']}, mask = True, det_uid = self.dets)
        else:
            self.exclude = self.tod.info.array_data.select_inner(
            {'row': exclude['rows'], 'col': exclude['cols']}, mask = True, det_uid = self.dets)
        if not(noExclude): self.exclude[list(self.tod.cuts.get_cut())] = True
        if 'det_uid' in dark:
            self.origDark = self.tod.info.array_data.select_inner(
                {'det_uid': dark['det_uid']}, mask = True, det_uid = self.dets)
        else:
            self.origDark = self.tod.info.array_data.select_inner(
            {'row': dark['rows'], 'col': dark['cols']}, mask = True, det_uid = self.dets)
        if 'det_uid' in liveCandidates:
            self.liveCandidates = self.tod.info.array_data.select_inner(
                {'det_uid': liveCandidates['det_uid']}, mask = True, det_uid = self.dets)
        else:
            self.liveCandidates = self.tod.info.array_data.select_inner(
            {'row': liveCandidates['rows'], 'col': liveCandidates['cols']}, mask = True, det_uid = self.dets)


        # Temperature params
        if self.params['findPathoParams']['thermParams']['autoTmax']:
            T_set = tod.get_hk(self.params['findPathoParams']['thermParams']["srv_setp"])
            dT_max = self.params['findPathoParams']['thermParams']['dT_max']
            self.params['findPathoParams']['thermParams']['T_max'] = T_set + 5 * dT_max


    def initializeCriteria( self ):
        p = self.params
        self.liveKeys = []
        self.darkKeys = []

        #Live Criteria
        for C in p["liveSelParams"].keys():
            K = "%sLive"%C
            self.liveKeys.append(K)
            if not(K in self.crit): self.crit[K] = {}
            self.crit[K].update(p["liveSelParams"][C])
            self.crit[K]["proc"] = True
            if "values" not in self.crit[K]:
                self.crit[K]["values"] = None

        for C in p["darkSelParams"].keys():
            K = "%sDark"%C
            self.darkKeys.append(K)
            if not(K in self.crit): self.crit[K] = {}
            self.crit[K].update(p["darkSelParams"][C])
            self.crit[K]["proc"] = True
            if "values" not in self.crit[K]:
                self.crit[K]["values"] = None

        self.activeLiveKeys = []
        for k in self.liveKeys:
            if (k.find("partial") == 0) and not(p["findPathoParams"].get("getPartial",False)):
                self.crit[k]["apply"] = False
            if self.crit[k]["apply"]:
                self.activeLiveKeys.append(k)

        self.activeDarkKeys = []
        for k in self.darkKeys:
            if self.crit[k]["apply"]:
                self.activeDarkKeys.append(k)

    def addCriterion( self, key, selection, apply = True ):
        self.liveKeys.append(key)
        self.crit[key] = {"sel": selection, "apply": apply, "proc": False}
        if apply: self.activeLiveKeys.append(key)

    def findPathologies( self, retrend = False, verbose = False,):
        """
        @brief Finds the common mode for both live and dark detectors and calculates the
        standard deviation of the detectors around the common modes, normalized by their
        white noise; it finds the correlation of all detectors to both common modes; and
        it finds fits the noise and finds the 1/f knees for all of them. It finally selects
        the live and dark detectors as isotated normal distributions in the previous
        calculations.
        @param  kneeHardLimit Maximum allowable 1/f knee frequency.
        @param  retrend       Retrend TOD after analysis
        @param  verbose       Show resulting number of selected detectors
        """
        assert self.tod.data is not None
        tictic = time.time()
        par = self.params['findPathoParams']

        # ANALYZE SCAN
        self.scan = analyzeScan( np.unwrap(self.tod.az), self.sampleTime,
                                **self.params.get("scanParams",{}) )
        self.scan_freq = self.scan["scan_freq"]
        self.chunkParams = {'T': self.scan["T"]*self.tod.info.downsample_level,
                            'pivot': self.scan["pivot"]*self.tod.info.downsample_level,
                            'N': self.scan["N"]}

        # ANALYZE TEMPERATURE
        self.Temp, self.dTemp, self.temperatureCut = checkThermalDrift(self.tod, par["thermParams"])

        # FIND ZERO DETECTORS
        if self.zeroSel is None:
            tic = time.time(); psLib.trace('moby', 2, "Finding Zero Detectors")
            self.zeroSel = ~self.tod.data[:,::100].any(axis=1)
            toc = time.time();
            psLib.trace('moby', 2, "It took %f seconds to find zero detectors"%(toc-tic))

        # GET CANDIDATE DETECTORS
        fullRMSlim = par.get("fullRMSlim",1e8)
        self.fullRMSsel = np.std(self.tod.data,axis=1) < fullRMSlim
        live = self.liveCandidates * ~self.zeroSel * self.fullRMSsel
        dark = self.origDark * ~self.zeroSel * self.fullRMSsel

        # Calibrate TOD to pW
        tic = time.time(); psLib.trace('moby', 2, "Calibrating")
        self.calibrate2pW()
        resp = self.calData["resp"]; ff = self.calData["ff"]
        cal = resp*ff
        if not(numpy.any(self.calData["calSel"])):
            psLib.trace('moby', 0, "ERROR: no calibration for this TOD")
            return 1
        toc = time.time();
        psLib.trace('moby', 2, "It took %f seconds to calibrate"%(toc-tic))


        # FIND JUMPS
        self.crit["jumpLive"]["values"] = moby2.libactpol.find_jumps(self.tod.data,
                            par['jumpParams']['dsStep'],
                            par["jumpParams"]["window"])
        self.crit["jumpDark"]["values"] = self.crit["jumpLive"]["values"]

        # FIND AMPLITUDE
        if 'ampLive' in self.crit:
            self.crit["ampLive"]["values"] = self.tod.data.max(axis=1) - self.tod.data.min(axis=1)

        # FREQUENCY SPACE ANALYSIS
        trend = moby2.tod.detrend_tod(self.tod)
        nf = nextregular(self.tod.nsamps)
        fdata = numpy.fft.rfft(self.tod.data, nf)
        dt = (self.tod.ctime[-1]-self.tod.ctime[0])/(self.tod.nsamps-1)
        df = 1./(dt*nf)

        # Low-frequency dark analysis
        res = multiFreqCorrAnal(fdata, dark, df, nf, self.ndata, self.scan_freq, par,
                          "darkCorrPar")
        self.preDarkSel = res["preSel"]
        self.crit["corrDark"]["values"] = res["corr"]
        self.crit["gainDark"]["values"] = res["gain"]
        self.crit["normDark"]["values"] = res["norm"]
        self.darkSel = self.preDarkSel.copy()

        # Low-frequency live analysis
        self.fbandSel = []
        self.fbands = []
        if par["liveCorrPar"].get("separateFreqs",False):
            fbs = np.array(list(set(self.tod.info.array_data["nom_freq"])))
            fbs = fbs[fbs != 0]
            for fb in fbs:
                self.fbandSel.append((self.tod.info.array_data["nom_freq"] == fb)*live)
                self.fbands.append(str(int(fb)))
        else:
            self.fbandSel.append(live)
            self.fbands.append("all")

        self.preLiveSel = np.zeros(self.ndet,dtype=bool)
        self.liveSel = np.zeros(self.ndet,dtype=bool)

        self.crit["darkRatioLive"]["values"] = np.zeros(self.ndet,dtype=float)
        self.crit["corrLive"]["values"] = np.zeros(self.ndet,dtype=float)
        self.crit["gainLive"]["values"] = np.zeros(self.ndet,dtype=float)
        self.crit["normLive"]["values"] = np.zeros(self.ndet,dtype=float)
        #self.darkCoeff = np.zeros(self.ndet,dtype=complex)
        #self.darkRatio = np.zeros(self.ndet,dtype=float)
        #self.preSelGroups = {}
        #self.preCC = {}
        self.multiFreqData = {}
        for fbs,fbn in zip(self.fbandSel,self.fbands):
            res = multiFreqCorrAnal(fdata, fbs, df, nf, self.ndata, self.scan_freq, par,
                              "liveCorrPar", darkSel=self.darkSel, tod=self.tod,
                              respSel = self.calData["respSel"], flatfield = self.flatfield_object)
            self.preLiveSel[fbs] = res["preSel"][fbs]
            self.liveSel[fbs] = res["preSel"][fbs]
            if 'darkRatio' in res:
                self.crit["darkRatioLive"]["values"][fbs] = res["darkRatio"][fbs]
            self.crit["corrLive"]["values"][fbs] = res["corr"][fbs]
            self.crit["gainLive"]["values"][fbs] = res["gain"][fbs]
            self.crit["normLive"]["values"][fbs] = res["norm"][fbs]
            # self.darkCoeff[fbs] = res["dcoeff"].squeeze()[fbs]
            # self.darkRatio[fbs] = res["ratio"][fbs]
            # self.preSelGroups[fbn] = res["groups"]
            # self.preCC[fbn] = res["cc"]
            self.multiFreqData[fbn] = res["all_data"]
        self.res = res
        # Undo flatfield correction
        self.crit["gainLive"]["values"] /= numpy.abs(ff)

        # Get slow Common Mode
        n_l = 1
        n_h = nextregular(int(round(par['driftFilter']/df))) + 1
        lf_data = fdata[:,n_l:n_h]
        fcmL = lf_data[self.preLiveSel].mean(axis = 0)
        fcmD = lf_data[self.preDarkSel].mean(axis = 0)
        self.dsCM, self.dsCM_dt = get_time_domain_modes(fcmL,n_l, self.ndata, df)
        self.dsDCM, _ = get_time_domain_modes(fcmD,n_l, self.tod.nsamps, df)

        # Add trend to CM
        trL = numpy.array(trend).T[self.preLiveSel].mean(axis=0)
        trLt = trL[:,np.newaxis]
        moby2.tod.retrend_tod(trLt, data = self.dsCM)
        trD = numpy.array(trend).T[self.preDarkSel].mean(axis=0)
        trDt = trD[:,np.newaxis]
        moby2.tod.retrend_tod(trDt, data = self.dsDCM)

        # Get Drift-Error
        DE = highFreqAnal(fdata, live, [n_l,n_h], self.ndata, nmodes = par["DEModes"],
                             highOrder = False, preSel = self.preLiveSel)
        self.crit["DELive"]["values"] = DE

        # Mid-frequency Analysis
        n_l = int(round(par["midFreqFilter"][0]/df))
        n_h = int(round(par["midFreqFilter"][1]/df))
        MFE = highFreqAnal(fdata, live, [n_l,n_h], self.ndata, nmodes = par["MFEModes"],
                             highOrder = False, preSel = self.preLiveSel)
        self.crit["MFELive"]["values"] = MFE

        # High-frequency analysis Live
        n_l = int(round(par["highFreqFilter"][0]/df))
        n_h = int(round(par["highFreqFilter"][1]/df))
        n_h = nextregular(n_h-n_l) + n_l
        if not(par["getPartial"]):
            rms, skewt, kurtt = highFreqAnal(fdata, live, [n_l,n_h], self.ndata,
                                           nmodes = par["HFLiveModes"],
                                           highOrder = True, preSel = self.preLiveSel)
        else:
            rms, skewt, kurtt, prms, pskewt, pkurtt = highFreqAnal(fdata, live,
                                           [n_l,n_h], self.ndata, nmodes = par["HFLiveModes"],
                                           highOrder = True, preSel = self.preLiveSel,
                                           scanParams = self.scan)
            self.crit["partialRMSLive"]["values"] = numpy.zeros([self.ndet,self.chunkParams["N"]])
            self.crit["partialSKEWLive"]["values"] = numpy.zeros([self.ndet,self.chunkParams["N"]])
            self.crit["partialKURTLive"]["values"] = numpy.zeros([self.ndet,self.chunkParams["N"]])
            self.crit["partialSKEWPLive"]["values"] = numpy.zeros([self.ndet,self.chunkParams["N"]])
            self.crit["partialKURTPLive"]["values"] = numpy.zeros([self.ndet,self.chunkParams["N"]])
            self.crit["partialRMSLive"]["values"][live] =  prms
            self.crit["partialSKEWLive"]["values"][live] = pskewt.T[:,0]
            self.crit["partialKURTLive"]["values"][live] = pkurtt.T[:,0]
            self.crit["partialSKEWPLive"]["values"][live] = pskewt.T[:,1]
            self.crit["partialKURTPLive"]["values"][live] = pkurtt.T[:,1]
        self.crit["rmsLive"]["values"] = rms
        self.crit["skewLive"]["values"] = numpy.zeros(self.ndet)
        self.crit["kurtLive"]["values"] = numpy.zeros(self.ndet)
        self.crit["skewpLive"]["values"] = numpy.zeros(self.ndet)
        self.crit["kurtpLive"]["values"] = numpy.zeros(self.ndet)
        self.crit["skewLive"]["values"][live] = skewt[0]
        self.crit["kurtLive"]["values"][live] = kurtt[0]
        self.crit["skewpLive"]["values"][live] = skewt[1]
        self.crit["kurtpLive"]["values"][live] = kurtt[1]

        # High-frequency analysis Dark
        rms = highFreqAnal(fdata, dark, [n_l,n_h], self.ndata, nmodes = par["HFDarkModes"],
                             highOrder = False, preSel = self.preDarkSel)
        self.crit["rmsDark"]["values"] = rms

        # Atmosphere -- 1/f analysis
        if par.get("fitPowerLaw",False):
            sel = self.preLiveSel + self.preDarkSel
            rms = self.crit["rmsLive"]["values"] + self.crit["rmsDark"]["values"]
            powLaw, level, knee = fit_atm(fdata, sel, dt, df, rms, self.scan_freq,
                                          **par.get("atmFit",{}))
            self.crit.update({"atmPowLaw":{"values":powLaw}})
            self.crit.update({"atmLevel":{"values":level}})
            self.crit.update({"atmKnee":{"values":knee}})

        # Bring dimensional observables back to RAW units
        recal = numpy.abs(cal)
        recal[cal==0] = 1.
        self.crit["normLive"]["values"] /= recal
        self.crit["rmsLive"]["values"] /= recal
        self.crit["DELive"]["values"] /= recal
        self.crit["MFELive"]["values"] /= recal
        if self.crit["partialRMSLive"]["values"] is not None:
            self.crit["partialRMSLive"]["values"] /= \
                numpy.resize(recal,self.crit["partialRMSLive"]["values"].T.shape).T

        # Retrend TOD
        if retrend:
            moby2.tod.retrend_tod(trend,self.tod)

        toctoc = time.time()
        dtime = (toctoc-tictic)/60
        psLib.trace('moby', 1, "It took %4.3f minutes to find pathologies." % dtime)
        return 0


    def makeNewSelections( self, params = None,  verbose = False):
        """
        @brief  Make detector selections (whole and partial).
        @param  verbose       Show resulting number of selected detectors
        @param  params        Dictionary with specific parameters (see __init__ for details.
        """
        # Update parameters
        if params is not None:
            self.params.update(params)
            self.initializeCriteria()

        if not(self.calibrated): self.calibrateValues()
        normalizeGains(self.crit["gainLive"]["values"],
                       sel = self.calData["stable"]*self.preLiveSel,
                       rejectOutliers = True, outlierSigma = 1.)

        # Make Live Selections
        for k in self.liveKeys:
            if self.crit[k]["proc"]: self._processSelection(k, self.preLiveSel)

        # Make Dark Selections
        for k in self.darkKeys:
            if self.crit[k]["proc"]: self._processSelection(k, self.preDarkSel)

        # SELECT LIVE DETECTORS
        self.liveSel = ~self.zeroSel*~self.exclude
        self.cutCounter = numpy.zeros(len(self.dets))
        self.cutCounter[self.zeroSel] += 1
        self.cutCounter[self.exclude] += 1

        live = numpy.ones(len(self.dets), dtype = bool)
        live[1024:] = False
        if verbose:
            psLib.trace('moby', 0, 'Number of live detectors cut out of 1024:')
            psLib.trace('moby', 0, '%-20s = %d'%('zero_cuts',len(self.dets[self.zeroSel*live])))
            psLib.trace('moby', 0, '%-20s = %d'%('exclude',len(self.dets[self.exclude*live])))

        for k in self.activeLiveKeys:
            self.liveSel *= self.crit[k]["sel"]
            self.cutCounter[~self.crit[k]["sel"]] += 1
            if verbose: psLib.trace('moby', 0,
                  '%-15s_cuts = %d'%(k,len(self.dets[~self.crit[k]["sel"]*live])))

        if self.crit["gainLive"]["apply"]:
            if self.crit["gainLive"]["sigma"] > self.params['otherParams']['gainCrit']:
                self.gainCut = True

        if self.params["otherParams"]["forceCalib"]:
            self.liveSel *= self.calData["calSel"]
            if verbose: psLib.trace('moby', 0,
                  '%-20s = %d'%("forceCalib_liveCuts",len(self.dets[self.calData["calSel"]])))

        # SELECT DARK DETECTORS
        self.darkSel = ~self.zeroSel

        for k in self.activeDarkKeys:
            sel = self.crit[k]["sel"]
            self.darkSel *= sel
            if verbose: psLib.trace('moby', 0,
                  '%-15s_darkCuts = %d'%(k,len(self.dets[~sel*self.origDark])))

        if self.params["otherParams"]['forceDark']:
            self.liveSel *= ~self.origDark
            self.darkSel *= self.origDark
            if verbose: psLib.trace('moby', 0,
                  '%-20s = %d'%("forceDark_Cuts",len(self.dets[self.origDark])))

        if verbose:
            psLib.trace('moby', 0, '%-20s = %d'%('# Live Detectors',len(self.dets[self.liveSel])))
            psLib.trace('moby', 0, '%-20s = %d'%('# Dark Detectors', len(self.dets[self.darkSel])))


    def _processSelection( self, key, initSel=None ):
        if initSel is None:
            initSel = numpy.ones(self.ndet,dtype = bool)
        p = self.crit[key]
        if p["values"] is None: return

        relSel, m, s = selectBySigma(p["values"], initSel, p['relSigma'])
        p.update({'median':m, 'sigma':s, 'relLims': (m-s*p['relSigma'], m+s*p['relSigma'])})
        if p["normalize"]:
            mm = numpy.median(p["values"][initSel])
            absSel = (p["values"]/mm >= p['absCrit'][0])*(p["values"]/mm <= p['absCrit'][1])
            p["abs_median"] = mm
        else:
            absSel = (p["values"] >= p['absCrit'][0])*(p["values"] <= p['absCrit'][1])

        if (key.find("partial") == 0):
            p["pSel"] = _combineSelTypes( relSel, absSel, p['selType'])
            p["sel"] = _findExcessPartial( p["pSel"], self.params['otherParams']['maxFracCut'])
        else:
            p["sel"] = _combineSelTypes( relSel, absSel, p['selType'])


    def getpWCalibration(self, flatfield = None):
        """
        @brief  Get responsivities and flatfield for calibration
        """
        if flatfield is None:
            flatfield = self.params["calibration"]["flatfield"]
        # Get responsivity
        resp = products.get_calibration(self.params["calibration"]["config"], self.tod.info)
        respSel = (resp.cal != 0.0)
        # Get flatfield
        self.flatfield_object = moby2.detectors.RelCal.from_dict(flatfield)
        ffSel, ff = self.flatfield_object.get_property('cal', det_uid=self.dets, default=1.)
        # Default resposibity to median of stable detectors
        _, stable = self.flatfield_object.get_property( 'stable', det_uid = self.dets, default=False)
        if self.params["calibration"].get("forceNoResp",False):
            rm = numpy.median(resp.cal[stable*respSel])
            resp.cal[~respSel] = rm
        if self.flatfield_object.calRMS is not None:
            _, ffRMS = self.flatfield_object.get_property(
                'calRMS', det_uid = self.dets, default=1.)
        else:
            ffRMS = np.zeros_like(self.dets)
        self.calData = {"resp": resp.cal,
                        "ff": ff,
                        "ffRMS": ffRMS,
                        "ffSel": ffSel,
                        "respSel": respSel,
                        "calSel": ffSel*respSel,
                        "stable": stable}
        return resp.cal, ff, ffRMS, respSel, ffSel, stable

    def calibrate2pW( self, flatfield = None, full = False):
        """
        @brief  Calibrate TOD to pW using responsivity and flatfield
        """
        if self.calData is None or flatfield is not None:
            self.getpWCalibration(flatfield = flatfield)
        factor = self.calData["resp"] * self.calData["ff"]
        # Apply to all except for original dark detectors
        s = ~self.origDark
        if not(self.calibratedTOD):
            moby2.libactpol.apply_calibration(self.tod.data,
                                          s.nonzero()[0].astype('int32'),
                                          factor[s].astype('float32'))
            self.calibratedTOD = True

    def calibrateValues( self, flatfield = None):
        """
        @brief   Apply calibration to those statistics that depend on data units.
        @param   flatfield  Calibration filename
        """
        if not(self.calibrated):
            resp, ff = self.getpWCalibration(flatfield = flatfield)[:2]
            dcal = numpy.where(~self.origDark)[0]

            # Downsampling effect
            aff = abs(ff[dcal])
            cal = resp[dcal]*aff
            if self.params['liveSelParams']['gain'].get('calibrate',True):
                self.crit["gainLive"]["values"][dcal] *= aff
            # if self.params['liveSelParams']['amp'].get('calibrate',True):
            #     self.crit["ampLive"]["values"][dcal] *= aff
            if self.params['liveSelParams']['norm'].get('calibrate',True):
                self.crit["normLive"]["values"][dcal] *= cal
            if self.params['liveSelParams']['rms'].get('calibrate',True):
                self.crit["rmsLive"]["values"][dcal] *= cal
            if self.params['liveSelParams']['DE'].get('calibrate',True):
                self.crit["DELive"]["values"][dcal] *= cal
            if self.params['liveSelParams']['MFE'].get('calibrate',True):
                self.crit["MFELive"]["values"][dcal] *= cal
            if self.params['liveSelParams']['partialRMS'].get('calibrate',True):
                if self.crit["partialRMSLive"]["values"] is not None:
                    self.crit["partialRMSLive"]["values"][dcal] *= \
                        numpy.resize(cal,self.crit["partialRMSLive"]["values"][dcal].T.shape).T
            # for i in range(len(dcal)):
            #     aff = abs(ff[dcal[i]])
            #     cal = resp[dcal[i]]*aff
            #     self.crit["gainLive"]["values"][dcal[i]] *= aff
            #     self.crit["normLive"]["values"][dcal[i]] *= cal
            #     self.crit["rmsLive"]["values"][dcal[i]] *= cal
            #     self.crit["DELive"]["values"][dcal[i]] *= cal
            #     self.crit["MFELive"]["values"][dcal[i]] *= cal
            #     if self.crit["partialRMSLive"]["values"] is not None:
            #         self.crit["partialRMSLive"]["values"][dcal[i]] *= \
            #             numpy.resize(cal,self.crit["partialRMSLive"]["values"][dcal[i]].T.shape).T
            self.calibrated = True
        else: psLib.trace('moby', 0, 'ERROR: pathologies already calibrated')


    def makeAzCuts( self, applyToTOD = False ):
        """
        @brief  Apply azimuth cuts to TOD
        """
#        if self.scan is None:
        sampleTime = ( (self.tod.ctime[-1]-self.tod.ctime[0])
                       / self.tod.nsamps )
        self.scan = analyzeScan(np.unwrap(self.tod.az), sampleTime,
                                **self.params.get("scanParams",{}) )

        c_obj = moby2.TODCuts.for_tod(self.tod)
        for d in self.tod.det_uid:
            c_obj.add_cuts(d, self.scan["az_cuts"])

        if applyToTOD:
            self.tod.cuts.merge_tod_cuts(c_obj)

        return c_obj

    def makeCuts( self, applyToTOD = False ):
        """
        @brief  Apply cuts to tod.
        """
        print('tod.dark not updated!')
        #self.tod.dark = self.dets[self.darkSel]
        #self.tod.ndark = len(self.tod.dark)
        c_obj = moby2.TODCuts.for_tod(self.tod, assign=False)
        c_obj.set_always_cut( numpy.nonzero(~self.liveSel)[0] )
        if applyToTOD:
            self.tod.cuts.merge_tod_cuts(c_obj)
        return c_obj

    def makePartialCuts( self, applyToTOD = False ):
        """
        @brief  Apply partial cuts from RMS, SKEW and KURT
        """
        SEL = numpy.ones([self.ndet, self.chunkParams['N']], dtype = bool)
        for k in self.activeLiveKeys:
            if (k.find("partial")==0):
                SEL *= self.crit[k]["pSel"]
        T = self.chunkParams['T']/self.tod.info.downsample_level
        pivot = self.chunkParams['pivot']/self.tod.info.downsample_level
        N = self.chunkParams['N']
        assert N == numpy.shape(SEL)[1]

        c_obj = moby2.TODCuts.for_tod(self.tod, assign=False)
        for d in self.dets[self.liveSel]:
            clist = []
            for i in range(N):
                if not(SEL[d][i]):
                    clist.append((pivot+i*T, pivot+(i+1)*T))
            if len(clist) > 0:
                c_obj.cuts[d] = moby2.tod.cuts.CutsVector(clist, c_obj.nsamps)

        if applyToTOD:
            self.tod.cuts.merge_tod_cuts(c_obj)
        return c_obj

    def showParams( self ):
        """
        @brief   Show parameters
        """
        printDictionary(self.params, 0)


    def viewChunksData( self, DATA, title = None, vmin = None, vmax = None, filename = None,
                        selection = None, display = True, units = None, interactive = True):
        """
        @brief   Produce 2D plot of chunk data
        """
        import_pylab(interactive = interactive)
        if DATA.dtype == 'bool':
            DATA = numpy.array(~DATA, dtype = int)
            vmin = 0; vmax = 1

        if selection is None: selection = numpy.ones(self.ndet, dtype = bool)
        m = pylab.matshow(DATA[selection].transpose(), vmin = vmin, vmax = vmax)
        b = pylab.colorbar(shrink=0.8)
        if units is not None: b.set_label(units)
        m.axes.set_aspect(len(self.dets[selection])/numpy.shape(DATA)[1])
        m.figure.set_size_inches(9.5, 8, forward = True)
        if title is not None: pylab.title(title)
        pylab.xlabel('Detector Number')
        pylab.ylabel('Chunk Number')
        if filename is not None: pylab.savefig(filename)
        if display: pylab.show()
        else: pylab.clf()


    def viewSelection( self, selection, title = None, filename = None, display = True,
                       interactive = True):
        """
        @brief  Visualize the location in the array of the detectors that would be cut or not
                given a selection of good detectors.
        @param  selection  Bool numpy array with False for detectors to be cut and True otherwise.
        @param  title      Specify the title to use in the plot.
        @param  filename   Specify the file name where the figure should be saved as a ".png"
        @param  display    True for showing and False for not showing the plot on screen.
        """
        import_pylab(interactive = interactive)
        fig = pylab.figure(1, figsize=(6,6))
        ax = pylab.subplot(111)
        if self.ndet < self.tod.info.array_data['det_uid'].size:
            x = self.tod.info.array_data['row']#numpy.arange(1056)/32
            y = self.tod.info.array_data['col']#numpy.arange(1056)%32
            l0 = pylab.plot(x,y,'s', label = 'Unread detectors')
            pylab.setp(l0, 'markeredgecolor', 'k', \
                           'markeredgewidth', 1.5, \
                           'markerfacecolor', 'w', \
                           'markersize', 6)
        l1 = pylab.plot(self.rows[selection],self.cols[selection], \
                        's', label = 'Uncut detectors')
        pylab.setp(l1, 'markeredgecolor', 'g', \
                       'markeredgewidth', 1.5, \
                       'markerfacecolor', 'w', \
                       'markersize', 6)
        l2 = pylab.plot(self.rows[~selection],self.cols[~selection], \
                        'rs', label = 'Cut detectors')
        pylab.setp(l2, 'markeredgecolor', 'r', \
                       'markeredgewidth', 1.5, \
                       'markerfacecolor', 'w', \
                       'markersize', 6)
        pylab.xlim([-1,self.Nrows])
        pylab.ylim([-1,self.Ncols])
        pylab.subplots_adjust(left = 0.13, right = 0.9, top = 0.9, bottom = 0.2)
        pylab.xlabel('rows')
        pylab.ylabel('cols')
        fo = matplotlib.font_manager.FontProperties(size=10)
        pylab.legend(loc = (0.3, -0.25), markerscale = 1, numpoints = 1, prop = fo)
        if title is not None: pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if display: pylab.show()
        else: pylab.clf()


    def comparativePlot( self, data, semilogy = False, xlim = None, ylim = None, \
                         title = None, ylabel = None, filename = None, display = True,
                         interactive = True):
        """
        @brief   Plot detector # versus value for various statistics.
        @param   data      data vector (ndet long) to plot.
        @param  semilogy  use logaritmic scale in Y axis.
        @param  xlim      specify range in X axis: [xmin,xmax]
        @param  ylim      specify range in Y axis: [ymin,ymax]
        @param  filename  string with name of the file where to save the figure.
        @param  title     string with title to use in figure.
        @param  ylabel    string with ylabel to use in figure.
        @param  display   display or not the image.
        """
        import_pylab(interactive = interactive)

        if semilogy:
            plot = pylab.semilogy
        else:
            plot = pylab.plot

        plot(self.dets, data, 'b.')
        plot(self.dets[self.liveSel], data[self.liveSel], 'g.')
        plot(self.dets[self.darkSel], data[self.darkSel], 'r.')
        pylab.xlabel('Detector #')
        if ylabel is not None:
            pylab.ylabel(ylabel)
        if xlim is not None:
            pylab.xlim(xlim)
        if ylim is not None:
            pylab.ylim(ylim)
        if title is not None:
            pylab.title(title)
        if filename is not None:
            pylab.savefig(filename)
        if display:
            pylab.show()
        else:
            pylab.clf()


    def plotArrayVector(self, data, selection = None, vmin = None, vmax = None,
                        title = None, filename = None, display = True, units = None,
                        interactive = True):
        """
        @brief   Plot the value for various statistics across the array.
        @param   data      data vector (ndet long) to plot.
        @param  selection bool array with selection of detectors to include in the plot.
        @param  vmin      minimum value in color range.
        @param  vmax      maximum value in color range.
        @param  title     string with title to use in figure.
        @param  filename  string with name of the file where to save the figure.
        @param  display   display or not the image.
        """
        import_pylab(interactive = interactive)

        if selection is None:
            selection = numpy.ones(len(self.dets), dtype = 'bool')
        data[~selection] = 0.0
        mat = numpy.reshape(data, [self.Nrows,self.Ncols]).transpose()
        m = pylab.matshow(mat, vmin = vmin, vmax = vmax)
        pylab.xlabel('Columns')
        m.axes.xaxis.set_label_position( 'top' )
        pylab.ylabel('Rows')
        b = pylab.colorbar(shrink = 0.8)
        if units is not None: b.set_label(units)
        if title is not None:
            pylab.title(title)
        if filename is not None:
            pylab.savefig(filename)
        if display:
            pylab.show()
        else:
            pylab.clf()


    def histogram(self, data, selection = None, bins = None, drange = None, \
                  title = None, xlabel = None, filename = None, display = True,
                  interactive = True):
        """
        @brief   Make histogram for the value for various statistics.
        @param   data      data vector (ndet long) to plot.
        @param  selection bool array with selection of detectors to include in the plot.
        @param  bins      number of bins to use.
        @param  title     string with title to use in figure.
        @param  filename  string with name of the file where to save the figure.
        @param  display   display or not the image.
        """
        import_pylab(interactive = interactive)

        if selection is None:
            selection = ~self.zeroSel
        dat = numpy.array(data)[selection].flatten()

        if drange is not None:
            dat = dat[(dat >= drange[0])*(dat <= drange[1])*~numpy.isnan(dat)]

        if bins is None:
            bins = min(max(int(len(dat)/20),10),100)

        a = pylab.hist(dat, bins = bins)
        if xlabel is not None: pylab.xlabel(xlabel)
        if title is not None: pylab.title(title)
        if filename is not None: pylab.savefig(filename)
        if display:
            if not(interactive): pylab.show()
        else: pylab.clf()

    def makeHistograms(self, outdir = "./", live = True, dark = True, nsig = 5):
        """
        @brief Makes histograms of all statistical values of active criteria
        """
        if not(self.calibrated): self.makeNewSelections()
        if live:
            sel = ~self.zeroSel*~self.origDark*~self.exclude
            for k in self.activeLiveKeys:
                if (not self.crit[k]["proc"]) or ("values" not in self.crit[k]): continue
                dat = self.crit[k]
                r = numpy.array([dat["median"]-nsig*dat["sigma"],dat["median"]+nsig*dat["sigma"]])
                fn = os.path.join(outdir,"hist_%s_%s.png"%(self.tod.info.name, k))
                if dat["normalize"]: mm = dat["abs_median"]
                else: mm = 1.
                self.histogram(dat["values"]/mm, drange=r/mm, filename = fn, display = False,
                               title = k, selection=sel)
        if dark:
            for k in self.activeDarkKeys:
                if not(self.crit[k]["proc"]): continue
                dat = self.crit[k]
                r = (dat["median"]-nsig*dat["sigma"],dat["median"]+nsig*dat["sigma"])
                fn = os.path.join(outdir,"hist_%s_%s.png"%(self.tod.info.name, k))
                if dat["normalize"]: mm = dat["abs_median"]
                else: mm = 1.
                self.histogram(dat["values"]/mm, drange=r, filename = fn, display = False,
                               title = k, selection = self.origDark)

    def quickReport( self, verbose = True ):
        """
        @brief  Give: Number of live detectors.
                      Fraction of uncut live time.
                      Number of dark detectors.
        """
        nLive = len(self.dets[self.liveSel])
        nDark = len(self.dets[self.darkSel])
        T = self.chunkParams['T']/self.tod.info.downsample_level
        N = self.chunkParams['N']
        tLiveTot = float(self.ndata)*len(self.dets[self.liveSel])
        partialSel = numpy.ones([self.ndet, N], dtype = bool)
        for k in self.activeLiveKeys:
            if (k.find("partial")==0):
                partialSel *= self.crit[k]["pSel"]
#        if self.params['partialCuts']['partialRMS']:
#            partialSel *= self.crit["partialRMS"]["pSel"]
#        if self.params['partialCuts']['partialSKEW']:
#            partialSel *= self.crit["partialSKEW"]["pSel"]
#        if self.params['partialCuts']['partialKURT']:
#            partialSel *= self.crit["partialKURT"]["pSel"]
        tDead = float(numpy.sum(numpy.array(~partialSel[self.liveSel], dtype = int))*T)
        if tLiveTot > 0:
            fLive = (tLiveTot - tDead)/tLiveTot
        else: fLive = 0.0
        if verbose:
            psLib.trace('moby', 0, 'Quick Cut Report:\n'
                                   '     # Live Detectors:    %4d\n'
                                   '     # Dark Detectors:    %4d\n'
                                   '     Fraction Live Time: %5.3f\n'%(nLive, nDark, fLive))
        return nLive, nDark, fLive



    def makeReport( self, filename ):
        """
        @brief   Generates a report with various statistical indicators.
        """
        nLive, nDark, fLive = self.quickReport( verbose = False )
        f = open(filename, "w")
        f.write("# Pathologies Report:   %5s\n\n"%self.name)
        f.write("# Live detectors:       %5d\n"%nLive)
        f.write("# Percentage time live: %5.1f\n"%(fLive*100))
        f.write("# Dark detectors:       %5d\n\n"%nDark)
        f.write("# Live candidates:      %5d\n"%len(self.dets[self.liveCandidates]))
        f.write("# Original dark:        %5d\n"%len(self.dets[self.origDark]))
        f.write("# Zero dectectors:      %5d\n"%len(self.dets[self.zeroSel]))
        f.write("# Good calibration:     %5d\n\n"%len(self.dets[self.calData["calSel"]]))
        f.write("# Criteria Live | Passed |  Abs(-)  |  Rel(-)  |  MeanV  |  Rel(+)  |  Abs(+)  |  Sigma  | Type \n")
        for k in self.activeLiveKeys:
            if not(self.crit[k]["proc"]): continue
            f.write("%14s  "%k)
            f.write("%7d  "%len(self.dets[self.crit[k]["sel"]]))
            if self.crit[k]["normalize"]: mm = self.crit[k]["abs_median"]
            else: mm = 1
            f.write("%10.3e  "%(self.crit[k]["absCrit"][0]*mm))
            f.write("%10.3e  "%(self.crit[k]["relLims"][0]))
            f.write("%10.3e  "%(self.crit[k]["median"]))
            f.write("%10.3e  "%(self.crit[k]["relLims"][1]))
            f.write("%10.3e  "%(self.crit[k]["absCrit"][1]*mm))
            f.write("%10.3e  "%(self.crit[k]["sigma"]))
            f.write("%s\n"%self.crit[k]["selType"])
        f.write("\n# Criteria Dark | Passed | Abs(-) | Rel(-) | MeanV | Rel(+) | Abs(+) | Sigma | Type \n")
        for k in self.activeDarkKeys:
            if not(self.crit[k]["proc"]): continue
            f.write("%14s  "%k)
            f.write("%7d  "%len(self.dets[self.crit[k]["sel"]]))
            if self.crit[k]["normalize"]: mm = self.crit[k]["abs_median"]
            else: mm = 1
            f.write("%10.3e  "%(self.crit[k]["absCrit"][0]*mm))
            f.write("%10.3e  "%(self.crit[k]["relLims"][0]))
            f.write("%10.3e  "%(self.crit[k]["median"]))
            f.write("%10.3e  "%(self.crit[k]["relLims"][1]))
            f.write("%10.3e  "%(self.crit[k]["absCrit"][1]*mm))
            f.write("%10.3e  "%(self.crit[k]["sigma"]))
            f.write("%s\n"%self.crit[k]["selType"])
        f.close()
        self.report = os.path.abspath(filename)

    def reportDetector(self, det):
        """
        @brief Quick report to see which criteria was applied to a certain detector
        @param det   detector to report (only one).
        """
        if not(isinstance(det, int)): det = int(self.tod.dets[det[0]][det[1]])
        psLib.trace('moby', 0, "Selection veredicts:")
        if self.origDark[det]:
            psLib.trace('moby', 0, "WARNING: This is a dark detector")
            for k in self.activeDarkKeys:
                if self.crit[k]["sel"][det]: psLib.trace('moby', 0, '    %11s: live'%k)
                else:   psLib.trace('moby', 0, '    %11s: dead'%k)
        else:
            for k in self.activeLiveKeys:
                if self.crit[k]["sel"][det]: psLib.trace('moby', 0, '    %11s: live'%k)
                else:   psLib.trace('moby', 0, '    %11s: dead'%k)
        if self.zeroSel[det]:      psLib.trace('moby', 0, '        zeroSel: dead')
        else:                      psLib.trace('moby', 0, '        zeroSel: live')


    def writeToPath( self, path ):
        """
        @brief Stores the pathologies object in Path.
        @param  path   directory path where to store the data object.
        """
        # Remove TOD from pathologies object

        tod = self.tod
        del self.tod
        data = self.__dict__.copy()
        data["calibratedTOD"] = data["calibrated"]

        if path[-1] == '/':
            path = path[0:-1]
        # filename = "%s.pickle" % path
        f = file( path, 'w' )
        p = pickle.Pickler( f, 2 )
        p.dump( data )
        f.close()

        self.tod = tod
        del tod

    @staticmethod
    def readFromPath( path, tod = None, params = None ):
        """
        @brief  Reloads a stored pathologies object.
        @param  path   Path to the directory where the data is stored.
        @param  tod    TOD that corresponds to the object you want to read.
        @return pa     Pathologies object with the data read.
        """
        if path[-1] == '/':
            path = path[0:-1]
        f = file( path, 'rb' )
        p = pickle.Unpickler( f )
        data = p.load()
        f.close()
        if "todName" in data:
            old = True
            todName = data.get("todName")
        else:
            old = False
            todName = data.get("name")
        assert tod.info.basename == todName, \
            "ERROR: TOD %s and stored name %s don't match" % \
            ( tod.info.name, data['todName'] )
        if params is None: params = data.get('params')
        else: data["params"] = params
        pa = Pathologies( tod, params )
        if old:
            pa.parseOldFile(data)
        else:
            pa.__dict__.update(data)
        pa.initializeCriteria()
        return pa


    @classmethod
    def read_from_path(cls, filename, tod=None, params = None):
        return cls.readFromPath(filename, tod, params)

    write_to_path = writeToPath


def selectBySigma( data, preSelection, thrld ):
    """
    @brief  Select the detectors that fall inside a normal distribution within a given standard deviation from the mean. The distribution mean and
            standard deviation are calculated given an initial set of detectors.
    @param  preSelection  bool array with preselection of detectors in distribution.
    @param  thrld         number of sigmas away from the mean to include in selection.
    """
    if data is None:
        return numpy.ones(len(preSelection), dtype = bool), 0.0, 0.0
    if len(data[preSelection]) == 0:
        return preSelection, 0., 0.
    datas = numpy.sort(data[preSelection].flatten())
    n = len(datas)
    if n == 0:
        return preSelection, 0., 0.
    q25 = datas[n//4]
    q75 = datas[n*3//4]
    s = 0.741*(q75-q25)
    m = datas[n//2]
    if s == 0: return numpy.ones(len(data), dtype = bool), m, s
    else: return (data >= m-thrld*s)*(data <= m+thrld*s), m, s


def _combineSelTypes( relSel, absSel, selType):
    """
    @brief Combine relative and partial cuts
    """
    if selType == 'relative':
        return relSel
    elif selType == 'absolute':
        return absSel
    elif selType == 'or':
        return relSel + absSel
    elif selType == 'and':
        return relSel*absSel
    else:
        psLib.trace("moby", 0, "ERROR: Unknown selection type")
        return numpy.ones(len(relSel), dtype = bool)

def _findExcessPartial( selection, threshold ):
    """
    @brief Find detectors with too many partial chunk cuts and select them to be cut.
    """
    n = numpy.shape(selection)[1]
    n_part = numpy.sum(numpy.array(~selection, dtype = int), axis = 1)
    sel = (numpy.array(n_part, dtype = float)/n <= threshold)
    return sel


class darkDets(object):
    """
    @brief  Class intended to output dark detector results in text format
    """
    _depot_structure = '{class}/{tag}/{night_name}/{tod_name}.txt'

    def __init__(self, patho = None):
        if patho is not None:
            self.darkDets = numpy.array(patho.dets[patho.darkSel], dtype = int)
        else: self.darkDets = None

    def applyToTOD( self, tod ):
        """
        brief   Updates tod.dark with the dark detectors from the object
        """
        tod.dark = list(self.darkDets)

    def writeToPath( self, path ):
        """
        @brief  Writes the list of dark detectors in the specified directory
        """
        if path[-1] == '/':
            path = path[0:-1]
#        filename = "%s/darkDets.txt" % path
        f = open(path, 'w')
        for d in self.darkDets:
            f.write('%d\n'%d)
        f.close()

    @staticmethod
    def readFromPath( path, tod ):
        """
        @brief  Read a list of dark detectors from the specified directory
        """
        if path[-1] == '/':
            path = path[0:-1]
#        filename = "%s/darkDets.txt" % path
        f = open(path, 'r')
        dd = darkDets()
        dd.darkDets = []
        for l in f:
            dd.darkDets.append(int(l.split('\n')[0].strip()))
        dd.darkDets = numpy.array(dd.darkDets, dtype = int)
        f.close()
        return dd


    @classmethod
    def read_from_path(cls, filename, tod=None):
        return cls.readFromPath(filename, tod)

    write_to_path = writeToPath

def get_detector_params(params):
    source =  params.get("source")
    if source == "matrix":
        mat = np.loadtxt(params.get("filename"), dtype = str)
        r,c = np.where(mat == "d") # Dark detectors
        dark = {}.fromkeys(("rows", "cols"))
        dark["rows"] = r; dark["cols"] = c
        r,c = np.where(mat == "b") # Bad detectors
        exclude = {}.fromkeys(("rows", "cols"))
        exclude["rows"] = r; exclude["cols"] = c
        r, c = np.where(mat == "s") # Stable detectors
        r2, c2 = np.where(mat == "c") # Candidate detectors
        liveCandidates = {}.fromkeys(("rows", "cols"))
        liveCandidates["rows"] = np.hstack([r,r2]); liveCandidates["cols"] = np.hstack([c,c2])
    elif source == "individual":
        # Excludable
        exclude = moby2.util.MobyDict.from_file(params.get("exclude"))
        # Dark
        dark = moby2.util.MobyDict.from_file(params.get("dark"))
        # liveCandidates
        liveCandidates = moby2.util.MobyDict.from_file(params.get("live"))
    else: raise ValueError("Unknown detector params source")
    return (exclude, dark, liveCandidates)


class darkModes(object):
    def __init__(self, pa = None):
        self.modes = None   # Downsamples and normalized
        self.dt = None      # Time between mode samples in seconds
        self.offset = 0.0   # Time until first sample in seconds
        self.emodes = None
        if pa is not None: self.fromPathologies(pa)

    def fromPathologies(self, pa):
        self.modes = pa.cmodes
        self.dt = pa.cmodes_dt
        self.offset = pa.offsets[0]*pa.sampleTime/pa.downsample4rms

    def expandModes(self, tod):
        nmod, smod = self.modes.shape
        tmod = np.arange(smod+1)*self.dt
        mod = numpy.hstack([self.modes,self.modes[:,0:1]])

        st = (tod.ctime[-1] - tod.ctime[0]) / (tod.nsamps-1)
        tt = tod.ctime - tod.ctime[0]

        # first index of reconstructed mode
        ind0 = numpy.max([int(self.offset/st) - tod.info.sample_index,0])
        # last sample of reconstructed mode
        ind1 = numpy.min([int((smod)*self.dt/st) + ind0, tod.nsamps])

        self.emodes = numpy.zeros((nmod,tod.nsamps))
        from scipy.interpolate import interp1d
        intermode = interp1d(tmod,mod,"cubic")
        self.emodes[:,ind0:ind1] = intermode(tt[ind0:ind1])
        self.emodes /= np.linalg.norm(self.emodes, axis=1)[:,numpy.newaxis]

    def deproject(self, tod, dets = None):
        if dets is None: dets = tod.det_uid
        if self.emodes is None: self.expandModes(tod)
        self.coeff = moby2.libactpol.data_dot_modes(tod.data,
                                  numpy.array(dets,dtype="int32"),
                                  numpy.asarray(self.emodes,dtype="float32"),
                                  moby2.util.encode_c(tod.cuts))
        moby2.libactpol.remove_modes( tod.data,
                                  numpy.array(dets, dtype='int32'),
                                  self.emodes.astype('float32'),
                                  self.coeff.astype('float64'))

        #for m in self.emodes:
        #    coeff = numpy.dot(tod.data[dets],m)
        #    tod.data[dets] -= numpy.outer(coeff,m)

    def write_to_path( self, path ):
        """
        @brief  Writes object in the specified directory
        """
        emodes = None
        if self.emodes is not None:
            emodes = self.emodes.copy()
            self.emodes = None
        if path[-1] == '/':
            path = path[0:-1]
        f = open(path, 'wb')
        numpy.save(f,self)
        f.close()
        self.emodes = emodes

    @staticmethod
    def readFromPath( path ):
        """
        @brief  Reads object from the specified directory
        """
        if path[-1] == '/':
            path = path[0:-1]
#        filename = "%s/???" % path
        f = open(path, 'rb')
        dm = numpy.load(f)
        f.close()
        return dm

    @classmethod
    def read_from_path(cls, filename):
        return cls.readFromPath(filename)



def multiFreqCorrAnal(fdata, sel, df, nf, nsamps, scan_freq, par, parTag,
                      darkSel=None, tod=None, respSel = None, flatfield = None):
    """
    Note that gainLive can be defined as the factor by which to multiply the common mode to fit a
    given detector signal.
    For this reason the calibration process divides each detector by this gain to equalize
    the array.
    Note that stable detectors here should be the same as the calibration fiducial detectors, but in
    this code they are defined as the common mode stable detectors.
    """
    from numpy import ma
    fmin   = par[parTag]["freqRange"].get("fmin",  0.017)
    fshift = par[parTag]["freqRange"].get("fshift",0.009)
    band   = par[parTag]["freqRange"].get("band",  0.070)
    Nwin   = par[parTag]["freqRange"].get("Nwin",      1)
    full   = par[parTag].get("fullReport", False)
    if not(par[parTag].get("forceResp",True)): respSel = None
    all_data = []
    psel=[]; corr=[]; gain=[]; norm=[]; darkRatio=[]
    fcm=[]; cm=[]; cmdt=[];
    fcmi=None; cmi = None; cmdti=None
    minFreqElem = 16
    for i in range(Nwin):
        n_l = int(round((fmin + i*fshift)/df))
        n_h = int(round((fmin + i*fshift + band)/df))
        if n_h - n_l < minFreqElem: n_h = n_l+minFreqElem
        if par[parTag].get("removeDark",False):
            if darkSel is None:
                print("ERROR: no dark selection supplied")
                return 0
            fcmi, cmi, cmdti = getDarkModes(fdata, darkSel, [n_l,n_h],
                                            df, nf, nsamps, par, tod)
            fcm.append(fcmi); cm.append(cmi); cmdt.append(cmdti)
        r = lowFreqAnal( fdata, sel, [n_l,n_h], df, nsamps, scan_freq, par.get(parTag,{}),
                    fcmodes = fcmi, respSel = respSel, flatfield = flatfield )
        psel.append(r["preSel"]); corr.append(r["corr"]); gain.append(np.abs(r["gain"]))
        norm.append(r["norm"]); darkRatio.append(r["ratio"])
        if full: all_data.append(r)
    spsel = np.sum(psel,axis=0)
    Nmax = spsel.max()
    #apsel = spsel == Nwin
    #psel50 = spsel >= Nwin/2.
    #for g in gain: g /= np.mean(g[apsel])
    psel50 = spsel >= Nmax/2.
    for g,s in zip(gain,psel): g /= np.mean(g[psel50*s])
    gain = np.array(gain)
    gain[np.isnan(gain)] = 0.
    #sg = np.sort(gain,axis=0)
    #gm = sg[Nwin/2]
    #gs = 0.741*(sg[(3*Nwin)/4] - sg[Nwin/4])
    #gs[gs==0] = 1.
    #gain_mask = np.abs(gain - np.repeat([gm],Nwin,axis=0)) > np.repeat([gs],Nwin,axis=0)
    #mgain = ma.MaskedArray(gain,gain_mask)
    mgain = ma.MaskedArray(gain,~np.array(psel))
    mgain_mean = mgain.mean(axis=0)
    mcorr = ma.MaskedArray(corr,~np.array(psel))
    #mcorr_mean = mcorr.mean(axis=0)
    mcorr_max = mcorr.max(axis=0)
    mnorm = ma.MaskedArray(norm,~np.array(psel))
    mnorm_mean = mnorm.mean(axis=0)
    res = {"preSel": psel50, "corr": mcorr_max.data, "gain": mgain_mean.data,
          "norm": mnorm_mean.data, "all_data": all_data}
    if par[parTag].get("removeDark",False):
        mdarkRatio = ma.MaskedArray(darkRatio,~np.array(psel))
        mdarkRatio_mean = mdarkRatio.mean(axis=0)
        res['darkRatio'] = mdarkRatio_mean.data
    return res

def getDarkModes(fdata,darkSel,frange,df,nf,nsamps,par,tod = None):
    """
    @brief Get dark or thermal modes from dark detectors and thermometer
           data.
    @return correlated modes in frequency and time domain, plus thermometer
            info.
    """
    n_l,n_h=frange
    fc_inputs = []
    # Dark detector drift
    if par["darkModesParams"].get("useDarks", False):
        dark_signal = fdata[darkSel,n_l:n_h].copy()
        fc_inputs.extend(list(dark_signal))
    # TEST CRYOSTAT TEMPERATURE
    if par["darkModesParams"].get("useTherm", False) and tod is not None:
        thermometers = []
        for channel in par['thermParams']['channel']:
            thermometer = tod.get_hk( channel, fix_gaps=True)
            if len(np.diff(thermometer).nonzero()[0]) > 0:
                thermometers.append(thermometer)
        if len(thermometers) > 0:
            thermometers = numpy.array(thermometers)
            fth = numpy.fft.rfft( thermometers, nf )[:,n_l:n_h]
            fc_inputs.extend(list(fth))
    elif par["darkModesParams"].get("useTherm", False) and tod is None:
            print("WARNING: TOD requiered to obtain thermometer data")
    fc_inputs = np.array(fc_inputs)
    if par.get("useTaper",False):
        taper = get_sine2_taper(frange, edge_factor = 6)
        fc_inputs *= np.repeat([taper],len(fc_inputs),axis=0)
    # Normalize modes
    fc_inputs /= np.linalg.norm(fc_inputs, axis=1)[:,np.newaxis]
    # Obtain main svd modes to deproject from data
    if par["darkModesParams"].get("useSVD",False):
        Nmodes = par["darkModesParams"].get("Nmodes",None)
        u, s, v = np.linalg.svd( fc_inputs, full_matrices=False )
        if Nmodes is None: fcmodes = v[s > s.max()/10]
        else: fcmodes = v[:Nmodes]
    else:
        fcmodes = fc_inputs

    # Get modes in time domain
    cmodes, cmodes_dt = get_time_domain_modes(
        fcmodes, n_l, nsamps, df)
    cmodes /= np.linalg.norm(cmodes, axis=1)[:,np.newaxis]
    return fcmodes, cmodes, cmodes_dt


def lowFreqAnal(fdata, sel, frange, df, nsamps, scan_freq, par,
                fcmodes=None, respSel = None, flatfield = None):
    """
    @brief Find correlations and gains to the main common mode over a frequency range
    """
    lf_data = fdata[sel,frange[0]:frange[1]].copy()
    ndet = len(sel)
    dcoeff = None
    ratio = None
    res = {}

    # Apply sine^2 taper to data
    if par.get("useTaper",False):
        taper = get_sine2_taper(frange, edge_factor = 6)
        lf_data *= np.repeat([taper],len(lf_data),axis=0)
    else:
        taper = np.ones(frange[1]-frange[0])

    # Deproject correlated modes
    if fcmodes is not None:
        data_norm = np.linalg.norm(lf_data,axis=1)
        dark_coeff = []

        for m in fcmodes:
            coeff = numpy.dot(lf_data.conj(),m)
            lf_data -= numpy.outer(coeff.conj(),m)
            dark_coeff.append(coeff)

        # Reformat dark coefficients
        if len(dark_coeff) > 0:
            dcoeff = numpy.zeros([len(dark_coeff),ndet],dtype=complex)
            dcoeff[:,sel] = np.array(dark_coeff)

        # Get Ratio
        ratio = numpy.zeros(ndet,dtype=float)
        data_norm[data_norm==0.] = 1.
        ratio[sel] = np.linalg.norm(lf_data,axis=1)/data_norm

    # Scan frequency rejection
    if par.get("cancelSync",False) and (scan_freq/df > 7):
        i_harm = get_iharm(frange, df, scan_freq, wide = par.get("wide",True))
        lf_data[:,i_harm] = 0.0

    # Get correlation matrix
    c = numpy.dot(lf_data,lf_data.T.conjugate())
    a = numpy.linalg.norm(lf_data,axis=1)
    aa = numpy.outer(a,a)
    aa[aa==0.] = 1.
    cc = c/aa

    # Get Norm
    ppar = par.get("presel",{})
    norm = numpy.zeros(ndet,dtype=float)
    fnorm = np.sqrt(np.abs(np.diag(c)))
    norm[sel] = fnorm*np.sqrt(2./nsamps)
    nnorm = norm/np.sqrt(nsamps)
    nlim = ppar.get("normLimit",[0.,1e15])
    if np.ndim(nlim) == 0: nlim = [0,nlim]
    normSel = (nnorm > nlim[0])*(nnorm < nlim[1])

    # Check if norms are divided in 2 groups. Use the higher norms
    sigs = ppar.get("sigmaSep",None)
    if sigs is not None:
        cent, lab = kmeans2(nnorm[normSel],2)
        frac = 0.2
        if lab.sum() > len(lab)*frac and lab.sum() < len(lab)*(1-frac):
            c0 = np.sort(nnorm[normSel][lab==0])
            c1 = np.sort(nnorm[normSel][lab==1])
            mc0 = c0[len(c0)//2]; mc1 = c1[len(c1)//2];
            sc0 = 0.741*(c0[(3*len(c0))//4] - c0[len(c0)//4])
            sc1 = 0.741*(c1[(3*len(c1))//4] - c1[len(c1)//4])
            sep =  (mc0 + sigs*sc0 - (mc1 - sigs*sc1))*np.sign(mc1-mc0)
            if sep < 0.:
                if mc1 > mc0:
                    normSel[normSel] *= (lab==1)
                else:
                    normSel[normSel] *= (lab==0)
        elif lab.sum() > 0:
            if lab.sum() > len(lab)//2:
                normSel[normSel] *= (lab==1)
            else:
                normSel[normSel] *= (lab==0)

    # Get initial detectors
    if respSel is None:
        respSel = np.ones(sel.shape,dtype=bool)
    if par.get("presel",{}).get("method","median") is "median":
        sl = presel_by_median(cc,sel=normSel[sel],
                              forceSel=respSel[sel],**par.get("presel",{}))
        res["groups"] = None
    elif par.get("presel",{}).get("method") is "groups":
        G, ind, ld, smap = group_detectors(cc, sel=normSel[sel], **par.get("presel",{}))
        sl = np.zeros(cc.shape[1],dtype=bool)
        sl[ld] = True
        res["groups"] = {"G": G, "ind": ind, "ld": ld, "smap": smap}
    else:
        raise ValueError("ERROR: Unknown preselection method")
    #if respSel is not None: sl *= respSel[sel]
    preSel = sel.copy()
    preSel[sel] = sl

    # Apply gain ratio in case of multichroic
    if (flatfield is not None) and ("scale" in flatfield.fields):
        scl = flatfield.get_property("scale",det_uid=np.where(sel)[0],
                                                 default = 1.)
        lf_data *= np.repeat([scl],lf_data.shape[1],axis=0).T

    # Get Correlations
    u, s, v = np.linalg.svd( lf_data[sl], full_matrices=False )
    corr = numpy.zeros(ndet)
    if par.get("doubleMode",False):
        corr[preSel] = np.sqrt(abs(u[:,0]*s[0])**2+abs(u[:,1]*s[1])**2)/fnorm[sl]
    else:
        corr[preSel] = np.abs(u[:,0])*s[0]/fnorm[sl]

    # Get Gains
    #
    # data = CM * gain
    #
    gain = numpy.zeros(ndet,dtype=complex)
    gain[preSel] = np.abs(u[:,0])  #/np.mean(np.abs(u[:,0]))
    res.update({"preSel": preSel, "corr": corr, "gain": gain, "norm": norm,
                "dcoeff": dcoeff, "ratio": ratio, "cc": cc, "normSel": normSel})
    return res


def highFreqAnal(fdata, sel, range, nsamps,
                 nmodes=0, highOrder = False, preSel = None,
                 scanParams = None):
    """
    @brief Find noise RMS, skewness and kurtosis over a frequency band
    """
    ndet = len(sel)
    hf_data = fdata[sel,range[0]:range[1]]
    if nmodes > 0:
        if preSel is None: preSel = numpy.ones(sel.sum(),dtype=bool)
        else: preSel = preSel[sel]
        c = numpy.dot(hf_data[preSel],hf_data[preSel].T.conjugate())
        u, w, v = numpy.linalg.svd(c, full_matrices = 0)
        kernel = v[:nmodes]/numpy.repeat([numpy.sqrt(w[:nmodes])],len(c),axis=0).T
        modes = numpy.dot(kernel,hf_data[preSel])
        coeff = numpy.dot(modes,hf_data.T.conj())
        hf_data -= numpy.dot(coeff.T.conj(),modes)
    rms = numpy.zeros(ndet)
    rms[sel] = numpy.sqrt(numpy.sum(abs(hf_data)**2,axis=1)/hf_data.shape[1]/nsamps)
    if highOrder:
        #hfd = numpy.hstack([numpy.expand_dims(numpy.zeros(len(hf_data)),1),
        #                    hf_data[:,:-1],
        #                    numpy.expand_dims(numpy.real(hf_data[:,-1]),1)])
        #hfd = numpy.fft.irfft(hfd)*numpy.sqrt(2.*hf_data.shape[1]/nsamps)
        hfd, _ = get_time_domain_modes( hf_data, 1, nsamps)
        skewt = stat.skewtest(hfd,axis=1)
        kurtt = stat.kurtosistest(hfd,axis=1)
        if scanParams is not None:
            T = scanParams["T"]
            pivot = scanParams["pivot"]
            N = scanParams["N"]
            f = float(hfd.shape[1])/nsamps
            t = int(T*f); p = int(pivot*f)
            prms = []; pskewt = []; pkurtt = []
            for c in range(N):
               prms.append(hfd[:,c*t+p:(c+1)*t+p].std(axis=1))
               pskewt.append(stat.skewtest(hfd[:,c*t+p:(c+1)*t+p],axis=1))
               pkurtt.append(stat.kurtosistest(hfd[:,c*t+p:(c+1)*t+p],axis=1))
            prms = numpy.array(prms).T
            pskewt = numpy.array(pskewt)
            pkurtt = numpy.array(pkurtt)
            return (rms, skewt, kurtt, prms, pskewt, pkurtt)
        else:
            return (rms, skewt, kurtt)
    else:
        return rms


def fit_atm(fdata, sel, dt, df, noise, scanf,
            fminA=0.2, fmaxA=3., fmaxT = 10.,
            its = 1, width = 0.005, **kwargs):
    """
    Fit a power law to the atmosphere signal in a range of frequencies.
    """
    scale = 2*dt**2*df
    delta = 0.7
    kneeM = fmaxA+delta
    ind_ini = int(fminA/df)
    ps = np.power(np.abs(fdata[sel,:int(fmaxT/df)]),2)*scale
    # Get scan harmonics
    n_harm = int(np.ceil(fmaxT/scanf))
    i_harm = np.array(np.round(np.arange(1,n_harm+1)*scanf/df), dtype=int)
    i_harmw = i_harm.copy()
    di = int(width/df)
    for i in range(di):
        i_harmw = np.hstack([i_harmw,i_harm-(i+1)])
        i_harmw = np.hstack([i_harmw,i_harm+(i+1)])
    i_harmw.sort()
    # Iterate range fit
    for it in range(its):
        fmax = kneeM-delta-fminA
        imin = int(fminA/df); imax = int(fmax/df)
        psr = ps[:,imin:imax]
        log_ps = np.log(psr)
        freq = np.arange(imin,imax)*df
        log_freq = np.log(freq)
        w = np.diff(log_freq)
        w = np.hstack([w,w[-1]])
        iharm = i_harmw[(i_harmw>imin)*(i_harmw<imax)]-imin
        s = np.ones(len(freq),dtype=bool)
        s[iharm] = False
        m,n = np.polyfit(log_freq[s], log_ps[:,s].T, 1, w=w[s])
        pl = np.power(freq[np.newaxis,s].repeat(ps.shape[0],0),
                      m[:,np.newaxis].repeat(len(freq[s]),1))*\
                      np.exp(n[:,np.newaxis].repeat(len(freq[s]),1))
        c = np.sum(psr[:,s]*pl,axis=1)/np.sum(pl*pl,axis=1)
        level = np.exp(n)*c
        knee = np.power(noise[sel]/level,1./m)
        kneeM = np.median(knee)
    mA = np.zeros(fdata.shape[0])
    levelA = np.zeros(fdata.shape[0])
    kneeA = np.zeros(fdata.shape[0])
    mA[sel] = m
    levelA[sel] = level
    kneeA[sel] = knee
    return mA, levelA, kneeA

def analyzeScan(az, dt = 0.002508,  N=50, vlim=0.01, qlim=0.01):
    """
    @brief Find scan parameters and cuts
    """
    # Find no motion
    lo, hi = ms.mquantiles(az, (qlim, 1-qlim))
    daz = np.diff(az)
    v_scan = np.r_[2*daz[0] - daz[1], daz]
    v_smooth = np.convolve(v_scan, np.ones((N,))/N)[(N-1)//2:-(N-1)//2]
    speed = numpy.median(abs(v_smooth))
    stop = abs(v_scan) < vlim*speed
    pick = abs(v_smooth) > 2*speed

    # Exit now in case of a stare TOD
    if all(stop):
        scan={"az_max":hi,"az_min":lo,"az_speed":speed,
              "scan_freq": 0.0, "az_cuts":None,
              "T":len(az), "pivot":0, "N":1}
        return scan

    # Find sections with no scanning
    noscan = stop*(az>lo)*(az<hi)

    # Get scan frequency
    faz = numpy.fft.rfft(az-az.mean())
    fscan = numpy.where(abs(faz)==abs(faz).max())[0][0]/dt/len(az)

    # Find turnarounds
    az_min = numpy.median(az[stop*(az<lo)])
    az_max = numpy.median(az[stop*(az>hi)])
    td = (abs(lo-az_min)+abs((az_max-hi)))

    # Find scan period parameters
    T_scan = int(1./fscan/dt) # number of samples in scan period
    T_ex = int(1.2*T_scan)
    onescan = az[~noscan][:T_ex]
    az_min0 = numpy.min(onescan[ stop[~noscan][:T_ex] * (onescan<lo) * (onescan>lo-td) ])
    az_max0 = numpy.max(onescan[ stop[~noscan][:T_ex] * (onescan>hi) * (onescan<hi+td) ])
    imin0 = numpy.where(az==az_min0)[0][0]
    imax0 = numpy.where(az==az_max0)[0][0]
    pivot = numpy.min([imin0,imax0]) # Index of first scan minima or maxima
    N_scan = (len(az)-pivot)/T_scan # Number of complete scan periods

    # Find cuts
    if hi - lo < 1. * np.pi / 180:
        flag = np.ones_like(az, dtype=bool)
    else:
        flag = (pick+stop)*(az>lo)*(az<hi)+(az<lo-td)+(az>hi+td)

    c_vect = moby2.tod.cuts.CutsVector.from_mask(flag).get_buffered(100)
    scan={"az_max":az_max,"az_min":az_min,"az_speed":speed,
          "scan_freq":fscan,"az_cuts":c_vect,
          "T":T_scan, "pivot":pivot, "N":N_scan}
    return scan


def checkThermalDrift(tod, thermParams):
    """
    @brief Measure the mean temperature and thermal drift, and
           suggest a thermalCut
    @return mean temperature, thermal drift and thermal cut flag
    """
    Temp = None
    dTemp = None
    temperatureCut = False
    channels = thermParams.get("channel",None)
    T_max = thermParams.get('T_max',None)
    dT_max = thermParams.get('dT_max',None)
    if channels is None or T_max is None or dT_max is None:
        return Temp, dTemp, temperatureCut
    thermometers = []
    for ch in channels:
        thermometer = tod.get_hk( ch, fix_gaps=True)
        if len(np.diff(thermometer).nonzero()[0]) > 0:
            thermometers.append(thermometer)
    if len(thermometers) > 0:
        thermometers = numpy.array(thermometers)
        # Get thermometer statistics
        th_mean = moby2.tod.remove_mean(data=thermometers)
        th_trend = moby2.tod.detrend_tod(data=thermometers)
        Temp = th_mean[0]
        dTemp = th_trend[1][0] - th_trend[0][0]
        if (Temp > T_max) or (abs(dTemp) > dT_max):
            temperatureCut = True
    return Temp, dTemp, temperatureCut


# Gain normalization code
def normalizeGains(gains, sel = None, useMedian = False,
                   rejectOutliers = False, outlierSigma = 2.,
                   weights = None, min_sel = 1, **kwargs):
    """
    Normalize a 1-d vector of gains as a function of the gains of a subset of detectors
    Parameters:
        sel: selection of detectors used to reference the normalization
        useMedian: use the median instead of the mean to set the normalization
        rejectOutliers: reject detectors with gains <outlierSigma> sigmas away from
                        the median gain of the reference detectors
        outlierSigma: number of sigmas to reject outliers
        weights: linear weights to average the gains of the reference detectors. This
                 option is not used if <useMedian> is True.
        min_sel: minimum number of valid detectors.
    Notes:
        MR1 cuts used weights = responsivity * flatfield, with the gains already multipied
        by the flatfield.
        Pathologies uses <rejectOutliers> = True and <sel> = None.
    """
    # Select reference detectors
    s = gains != 0
    if sel is not None:
        if sel.sum() < min_sel: return 1, False
        s *= sel
    if rejectOutliers:
        sgains = np.sort(gains[s])
        N = len(sgains)
        med = sgains[N//2]
        sig = 0.741*(sgains[(3*N)//4] - sgains[N//4])
        s *= np.abs(gains - med) < sig*outlierSigma
    if s.sum() < min_sel: return 1, False
    if useMedian: norm = np.median(gains[s])
    else:
        if weights is None: w = 1./s.sum()
        else: w = weights[s]/weights[s].sum()
        norm = np.sum(w*gains[s])
    gains /= norm
    return norm, True
