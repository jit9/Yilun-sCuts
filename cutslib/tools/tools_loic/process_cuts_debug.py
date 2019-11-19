import sys, os
import time
import numpy as np

import moby2.util.noninteractive_plots

import moby2
from moby2.scripting import products
from moby2.analysis.tod_ana import pathologies
from moby2.analysis import hwp
from optparse import OptionParser

o = OptionParser()
o.add_option('-v','--verbosity',type=int,default=1)
o.add_option('-o','--output-prefix',default=None)
o.add_option('-i','--interactive-errors',action='store_true')
opts, args = o.parse_args()




#
# Load parameters, update them from command line as needed
#

obs = args[0]
params = moby2.util.MobyDict.from_file(args[1])
IV_wait = params.get("IV_wait",0)
rebias_wait = params.get("rebias_wait",0)
end = params.get("sample_end",None)
if end is not None:
    end = -end


depot = moby2.util.Depot(params.get('depot'))
cutParams = params.get('cutParams')
cpar = moby2.util.MobyDict.from_file(cutParams)
pathop = cpar['pathologyParams']
glitchp = cpar['glitchParams']

# Fill cuts parameter
no_noise = not(cpar.get("fillWithNoise",True))

user_config = moby2.util.get_user_config()
moby2.pointing.set_bulletin_A(params=user_config.get('bulletin_A_settings'))

# Find IV/rebias gap time
ct = int(obs.split("/")[-1].split(".")[0])
ctimes = (ct-IV_wait,ct)
from moby2.instruments import actpol as inst
config_file = params.get('manifest_conf', None)
db = inst.TODDatabase(config_file=config_file)
recs = db.select_acqs(suffix='iv', ctime=ctimes)
if len(recs) > 0:
    offset_IV = (IV_wait - ct + recs[-1].ctime)*400
else:
    offset_IV = 0
ctimes = (ct-rebias_wait,ct)
recs = db.select_acqs(suffix='bc1_step', ctime=ctimes)
if len(recs) > 0:
    offset_rebias = (rebias_wait - ct + recs[-1].ctime)*400
else:
    offset_rebias = 0
offset = max( offset_IV, offset_rebias, params.get("offset",0)*400 )
print "Offset set to %is" %offset



# READ TOD
tic1 = time.time()
params_loadtod = {
    'filename':obs,
    'fix_sign':True,
    'start':offset,
    'end':end,
    'repair_pointing':True
}
params_loadtod.update(params.get("params_loadtod",{}))
tod = moby2.scripting.get_tod(params_loadtod)
dt = (time.time() - tic1) / 60
minFreq = pathop["findPathoParams"]['liveCorrPar']['freqRange']['fmin']
minPeriods = pathop["findPathoParams"].get("minPeriods", 1.)
minTODlength = 1 / minFreq * minPeriods
print "It took %6.2f minutes to load TOD"%dt
print "Processing TOD %s" % tod.info.name
print "ndata = %d"%tod.data.shape[1]
if minTODlength > tod.ctime[-1] - tod.ctime[0]:
    raise RuntimeError, "TOD too short to perform LF analysis"
tic2 = time.time()

# CUT MCE FRAME ERROR
mce_cuts = moby2.tod.get_mce_cuts(tod)
moby2.tod.fill_cuts(tod, mce_cuts, no_noise = no_noise)

# CUT SOURCES
if cpar.get_deep(('source_cuts','source_list'),None) is not None:
    print "Finding new source cuts"
    tod.fplane = products.get_focal_plane(cpar['pointing'], tod.info)
    f = open(cpar['source_cuts']['source_list'], 'r')
    source_list = f.readlines()
    source_list = [(s.strip('\n'), 'source') for s in source_list]
    f.close()
    matched_sources = moby2.ephem.get_sources_in_patch(
        tod=tod, source_list=source_list)
    mask_params = cpar.get_deep(('source_cuts','mask_params'),{})
    shift_params = cpar.get_deep(('source_cuts', 'mask_shift_generator'))
    if shift_params is not None:
        mask_params['offset'] = products.get_pointing_offset(
            shift_params, tod=tod, source_offset=True)
        if mask_params['offset'] is None:
            mask_params['offset'] = (0.,0.)
        if max(mask_params['offset']) > 20./60:
            mask_params['map_size'] = max(mask_params['offset']) + 10./60
    pos_cuts_sources = moby2.TODCuts.for_tod(tod, assign=False)
    pos_cut_dict = {}
    print matched_sources
    for source in matched_sources:
        pos_cut_dict[source[0]] = moby2.tod.get_source_cuts(
            tod, source[1], source[2], **mask_params)
        pos_cuts_sources.merge_tod_cuts(pos_cut_dict[source[0]])
    moby2.tod.fill_cuts(tod, pos_cuts_sources, no_noise = no_noise)

# CUT PLANETS
if params.get("cut_planets",False):
    print "Finding new planet cuts"
    if not hasattr(tod,'fplane'):
        tod.fplane = products.get_focal_plane(cpar['pointing'], tod.info)
    matched_sources = moby2.ephem.get_sources_in_patch(
        tod=tod, source_list=None)
    mask_params = cpar.get_deep(('planet_cuts','mask_params'),{})
    shift_params = cpar.get_deep(('planet_cuts', 'mask_shift_generator'))
    if shift_params is not None:
        mask_params['offset'] = products.get_pointing_offset(
            shift_params, tod=tod, source_offset=True)
        if mask_params['offset'] is None:
            mask_params['offset'] = (0.,0.)
        if max(mask_params['offset']) > 20./60:
            mask_params['map_size'] = max(mask_params['offset']) + 10./60
    pos_cuts_planets = moby2.TODCuts.for_tod(tod, assign=False)
    pos_cut_dict = {}
    print matched_sources
    for source in matched_sources:
        pos_cut_dict[source[0]] = moby2.tod.get_source_cuts(
            tod, source[1], source[2], **mask_params)
        pos_cuts_planets.merge_tod_cuts(pos_cut_dict[source[0]])
    moby2.tod.fill_cuts(tod, pos_cuts_planets, no_noise = no_noise)


# PARTIAL CUTS
tic5 = time.time()
# Generate and save new glitch cuts (note calbol may not be implemented...)
cuts_partial = moby2.tod.get_glitch_cuts(tod=tod, params=glitchp)
cuts_partial.merge_tod_cuts(mce_cuts)
moby2.tod.fill_cuts(tod, cuts_partial, extrapolate = False, no_noise = no_noise)
tod.cuts = cuts_partial

# REMOVE MEAN
tic4 = time.time()
if params.get("remove_median", True): moby2.tod.remove_median(tod)
else: moby2.tod.remove_mean(tod)
if params.get("detrend", False): moby2.tod.detrend_tod(tod)
if params.get("remove_filter_gain", False): moby2.tod.remove_filter_gain(tod)

# DOWNSAMPLE
tic6 = time.time()
n_downsample = params.get("n_downsample")
if n_downsample is not None:
    tod = tod.copy(resample=2**n_downsample, resample_offset=1)
    
# FIND PATHOLOGIES
tic7 = time.time()
pa = pathologies.Pathologies(tod, pathop, noExclude = True)




####### Here starts findpathologies ########
############################################

par = pa.params['findPathoParams']

# ANALYZE SCAN                                                                                    
pa.scan = pathologies.analyzeScan( np.unwrap(pa.tod.az), pa.sampleTime,
                         **pa.params.get("scanParams",{}) )
pa.scan_freq = pa.scan["scan_freq"]
pa.chunkParams = {'T': pa.scan["T"]*pa.tod.info.downsample_level,
                    'pivot': pa.scan["pivot"]*pa.tod.info.downsample_level,
                    'N': pa.scan["N"]}


# ANALYZE TEMPERATURE                                                                             
pa.Temp, pa.dTemp, pa.temperatureCut = pathologies.checkThermalDrift(pa.tod, par["thermParams"])


# FIND ZERO DETECTORS                                                                             
if pa.zeroSel is None:
    pa.zeroSel = ~pa.tod.data[:,::100].any(axis=1)


# GET CANDIDATE DETECTORS                                                                         
fullRMSlim = par.get("fullRMSlim",1e8)
pa.fullRMSsel = np.std(pa.tod.data,axis=1) < fullRMSlim
live = pa.liveCandidates * ~pa.zeroSel * pa.fullRMSsel
dark = pa.origDark * ~pa.zeroSel * pa.fullRMSsel


# Calibrate TOD to pW                                                                             
pa.calibrate2pW()
resp = pa.calData["resp"]; ff = pa.calData["ff"]
cal = resp*ff
if not(np.any(pa.calData["calSel"])):
    print "ERROR: no calibration for this TOD"
    stop


# FIND JUMPS                                                                                      
pa.crit["jumpLive"]["values"] = moby2.libactpol.find_jumps(
    pa.tod.data,
    par['jumpParams']['dsStep'],
    par["jumpParams"]["window"])
pa.crit["jumpDark"]["values"] = pa.crit["jumpLive"]["values"]


# FREQUENCY SPACE ANALYSIS                                                                        
trend = moby2.tod.detrend_tod(pa.tod)
nf = pathologies.nextregular(pa.tod.nsamps)
fdata = np.fft.rfft(pa.tod.data, nf)
dt = (pa.tod.ctime[-1]-pa.tod.ctime[0])/(pa.tod.nsamps-1)
df = 1./(dt*nf)



# Low-frequency dark analysis                                                                     
res = pathologies.multiFreqCorrAnal(fdata, dark, df, nf, pa.ndata, pa.scan_freq, par,
                        "darkCorrPar")
pa.preDarkSel = res["preSel"]
pa.crit["corrDark"]["values"] = res["corr"]
pa.crit["gainDark"]["values"] = res["gain"]
pa.crit["normDark"]["values"] = res["norm"]
pa.darkSel = pa.preDarkSel.copy()


# Low-frequency live analysis                                                                     
pa.fbandSel = []
pa.fbands = []
if par["liveCorrPar"].get("separateFreqs",False):
    fbs = np.array(list(set(pa.tod.info.array_data["nom_freq"])))
    fbs = fbs[fbs != 0]
    for fb in fbs:
        pa.fbandSel.append((pa.tod.info.array_data["nom_freq"] == fb)*live)
        pa.fbands.append(str(int(fb)))
else:
    pa.fbandSel.append(live)
    pa.fbands.append("all")

pa.preLiveSel = np.zeros(pa.ndet,dtype=bool)
pa.liveSel = np.zeros(pa.ndet,dtype=bool)
 
pa.crit["darkRatioLive"]["values"] = np.zeros(pa.ndet,dtype=float)
pa.crit["corrLive"]["values"] = np.zeros(pa.ndet,dtype=float)
pa.crit["gainLive"]["values"] = np.zeros(pa.ndet,dtype=float)
pa.crit["normLive"]["values"] = np.zeros(pa.ndet,dtype=float)
pa.multiFreqData = {}


fbs = pa.fbandSel[0]
fbn = pa.fbands[0]



#### Here starts Multifreqcorranal ####
#######################################


# fmin   = par["liveCorrPar"]["freqRange"].get("fmin",  0.017)
# fshift = par["liveCorrPar"]["freqRange"].get("fshift",0.009)
# band   = par["liveCorrPar"]["freqRange"].get("band",  0.070)
# Nwin   = par["liveCorrPar"]["freqRange"].get("Nwin",      1)
# full   = par["liveCorrPar"].get("fullReport", False)
nsamps = pa.ndata
fmin   = 0.017
fshift = 0.009
band   = 0.070
Nwin   = 1
full   = False

if not(par["liveCorrPar"].get("forceResp",True)): respSel = None
all_data = []
psel=[]; corr=[]; gain=[]; norm=[]; darkRatio=[]
fcm=[]; cm=[]; cmdt=[];
fcmi=None; cmi = None; cmdti=None
minFreqElem = 16

i = 0
n_l = int(round((fmin + i*fshift)/df))
n_h = int(round((fmin + i*fshift + band)/df))
if n_h - n_l < minFreqElem: n_h = n_l+minFreqElem
if par['liveCorrPar'].get("removeDark",False):
    if darkSel is None:
        print "ERROR: no dark selection supplied"
        stop
        fcmi, cmi, cmdti = pathologies.getDarkModes(fdata, darkSel, [n_l,n_h],
                                        df, nf, nsamps, par, tod)
    fcm.append(fcmi); cm.append(cmi); cmdt.append(cmdti)

#### Here starts lowfreqanal ####
#################################

par = par['liveCorrPar']
sel = fbs
lf_data = fdata[sel,n_l:n_h].copy()
ndet = len(sel)
dcoeff = None
ratio = None
res = {}

# Apply sine^2 taper to data                                                                          
if par.get("useTaper",False):
    taper = pathologies.get_sine2_taper([n_l,n_h], edge_factor = 6)
    lf_data *= np.repeat([taper],len(lf_data),axis=0)
else:
    taper = np.ones(n_h-n_l)
    

fcmodes = fcm
# Deproject correlated modes                                                                          
if fcmodes is not None:
    data_norm = np.linalg.norm(lf_data,axis=1)
    dark_coeff = []
    
    for m in fcmodes:
        coeff = np.dot(lf_data.conj(),m)
        lf_data -= np.outer(coeff.conj(),m)
        dark_coeff.append(coeff)

    # Reformat dark coefficients                                                                      
    if len(dark_coeff) > 0:
        dcoeff = np.zeros([len(dark_coeff),ndet],dtype=complex)
        dcoeff[:,sel] = np.array(dark_coeff)

    # Get Ratio                                                                                       
    ratio = np.zeros(ndet,dtype=float)
    data_norm[data_norm==0.] = 1.
    ratio[sel] = np.linalg.norm(lf_data,axis=1)/data_norm

# Scan frequency rejection                                                                            
if par.get("cancelSync",False) and (pa.scan_freq/df > 7):
    i_harm = pathologies.get_iharm([n_l,n_h], df, pa.scan_freq, wide = par.get("wide",True))
    lf_data[:,i_harm] = 0.0

# Get correlation matrix                                                                              
c = np.dot(lf_data,lf_data.T.conjugate())
a = np.linalg.norm(lf_data,axis=1)
aa = np.outer(a,a)
aa[aa==0.] = 1.
cc = c/aa

# Get Norm                                                                                            
ppar = par.get("presel",{})
norm = np.zeros(ndet,dtype=float)
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
        mc0 = c0[len(c0)/2]; mc1 = c1[len(c1)/2];
        sc0 = 0.741*(c0[(3*len(c0))/4] - c0[len(c0)/4])
        sc1 = 0.741*(c1[(3*len(c1))/4] - c1[len(c1)/4])
        sep =  (mc0 + sigs*sc0 - (mc1 - sigs*sc1))*np.sign(mc1-mc0)
        if sep < 0.:
            if mc1 > mc0:
                normSel[normSel] *= (lab==1)
            else:
                normSel[normSel] *= (lab==0)
    elif lab.sum() > 0:
        if lab.sum() > len(lab)/2:
            normSel[normSel] *= (lab==1)
        else:
            normSel[normSel] *= (lab==0)
        


respSel = pa.calData["respSel"]
# Get initial detectors                                                                               
if respSel is None:
    respSel = np.ones(sel.shape,dtype=bool)
if par.get("presel",{}).get("method","median") is "median":
    sl = pathologies.presel_by_median(cc,sel=normSel[sel],
                          forceSel=respSel[sel],**par.get("presel",{}))
    res["groups"] = None
elif par.get("presel",{}).get("method") is "groups":
    G, ind, ld, smap = pathologies.group_detectors(cc, sel=normSel[sel], **par.get("presel",{}))
    sl = np.zeros(cc.shape[1],dtype=bool)
    sl[ld] = True
    res["groups"] = {"G": G, "ind": ind, "ld": ld, "smap": smap}
else:
    raise "ERROR: Unknown preselection method"
    #if respSel is not None: sl *= respSel[sel]                                                           
preSel = sel.copy()
preSel[sel] = sl







