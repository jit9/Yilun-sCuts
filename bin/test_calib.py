import pandas as pd
import moby2
from moby2.scripting import products
from moby2.scripting.pathologies_tools import get_pwv
import os
import numpy as np
from scipy import stats
from moby2.util.database import TODList
from matplotlib import pyplot as plt
import ephem
import sys, cPickle
import pyfits as pf
T_URANUS = 177.6        #K CMB    


s15_split_time = 1441810800

def get_beam_solid_angle(array,season):
    bso = {'s13': {'PA1':195.3},
           's14': {'PA1':189.0, 'PA2':180.1},
           's15a': {'PA1':207.2, 'PA2':193.6, 'PA3_90': 525.0, 'PA3_150':291.8},
           's15b': {'PA1':195.1, 'PA2':182.8, 'PA3_90': 493.9, 'PA3_150':261.5},
           's16': {'PA2':179.8, 'PA3_90':471.4, 'PA3_150':241.5, 'PA4_150':236.1, 'PA4_220':112.8},
           's17': {'PA4_150':218.2, 'PA4_220':116.8, 'PA5_90':428.1, 'PA5_150':201.1, 'PA6_90':463.7, 'PA6_150':215.5},
       }
    return bso[season][array] * 1e-9


def get_planet_solid_angle(ctime):
    if np.asarray(ctime).ndim > 0:
        return np.array([get_planet_solid_angle(t) for t in ctime])
    U = ephem.Uranus(moby2.util.ctime.to_string(ctime,format='%Y/%m/%d'))
    return np.pi*U.radius**2


df = pd.DataFrame()

# Initialize parameters
params = moby2.util.MobyDict.from_file(sys.argv[1])
fb = moby2.scripting.get_filebase()
depot = moby2.util.Depot(params.get('depot'))
if params.has_key('tag'):
    tag_calib = params.get('tag')
    tag_cuts = params.get('tag')
    tag_selected = params.get('tag')
else:
    tag_calib = params.get('tag_calib')
    tag_cuts = params.get('tag_cuts')
    tag_selected = params.get('tag_selected')
if params.has_key('flatfield'):
    FF = moby2.util.MobyDict.from_file(params.get('flatfield'))
tods, flag = np.loadtxt(
    '/data/actpol/depot_yguan/SelectedTODs/%s/selectedTODs_uranus.txt' %tag_selected, dtype=str, usecols=[0,5]).T
df['tods'] = np.asarray( tods, dtype=str )
df['flag'] = np.asarray( flag, dtype=int )
df = df[df.flag==2]

print "Start with %i selected TODs" %df.tods.size 
path = params.get('planet_path')
array_name = params.get('array')
array_data =  products.get_array_data(
    {'instrument':'actpol','array_name':'ar%s' %array_name[2],'season':'20%s' %params.get('season')[-2:]})
ctime_start = params.get('ctime_start', 0)
loading_max = params.get('loading_max', 10)
alma_pwv = moby2.aux_data.alma.WeatherChannel('pwv')
apex_pwv = moby2.aux_data.apex.WeatherChannel('radiometer')


# Read the data
if not os.path.exists(os.path.join('PICKLE', tag_calib+'.pickle')):
    cal = np.zeros((2048,len(df.tods)))
    peak_DAC = np.zeros((2048,len(df.tods)))
    alt = np.zeros(len(df.tods))
    pwv = np.zeros(len(df.tods))
    hour_utc = np.zeros(len(df.tods))
    ctime = np.zeros(len(df.tods), dtype=int)
    for idx, obs in enumerate(df.tods):
        obs = (obs.split('/')[-1])
        if ".zip" in obs: obs = obs[:-4]
        print ""
        print "%s - %i/%i" %(obs, idx, len(df.tods))
        try:
            filename = path+obs+'_fp.fits'
            print filename
            data = pf.getdata(filename)
            det_uid = data['det_uid']
            peak_DAC[det_uid,idx] = data['peak_dac']
        except:
            print 'No peak measurement'
            pass
        filename = fb.filename_from_name(obs, single=True)
        try:
            tod = moby2.scripting.get_tod(
                {'filename':filename, 'read_data':False})
        except: 
            print "Can't load tod %s" %filename
            continue
        alt[idx] = tod.alt.mean()
        pwv_, tdiff = apex_pwv.get_nearest(tod.info.ctime)
        print "Apex tdiff = %i" %tdiff
        if np.abs(tdiff) < 600.: 
            pwv[idx] = pwv_
        else:
            pwv_, tdiff = alma_pwv.get_nearest(tod.info.ctime)
            print "ALMA tdiff = %i" %tdiff
            if np.abs(tdiff) < 600.:
                pwv[idx] = pwv_
        hour_utc[idx] = ( tod.info.ctime % 86400. ) / 3600
        ctime[idx] = tod.info.ctime
        try:
            Cal = depot.read_object(
                moby2.Calibration, tag=tag_calib, tod=tod)
            cal[Cal.det_uid,idx] = np.abs(Cal.cal) * tod.info.runfile.ReadoutFilter().gain()
            if params.has_key('flatfield'):
                cal[FF['det_uid']] *= FF['cal']
        except: 
            print 'No cal for this TOD'
        try:
            cuts = depot.read_object(moby2.TODCuts,
                                     tag=tag_cuts,
                                     tod=tod)
            if cuts.get_uncut().size < 0:
                cal[:,idx] = 0
            else:
                dc = cuts.get_cut()
                cal[dc,idx] = 0
        except:
            print 'No cuts for this TOD'

    data = {}
    data["cal"] = cal
    data["peak_DAC"] = peak_DAC
    data["alt"] = alt
    data["pwv"] = pwv
    data["hour_utc"] = hour_utc
    data["ctime"] = ctime
    data["tods"] = tods
    f = open("PICKLE/%s.pickle"%tag_calib,"w")
    p = cPickle.Pickler(f)
    p.dump(data)
    f.close()
    df['alt'] = data['alt']
    df['pwv'] = data['pwv']
    df['hour_utc'] = data['hour_utc']
    df['ctime'] = data['ctime']
    df.index = pd.to_datetime(df.ctime, unit='s')
else:
    print "Data already exist, load them from pickle"
    f = open("PICKLE/%s.pickle"%tag_calib,"r")
    data = cPickle.load(f)
    f.close()
    cal = data["cal"]
    peak_DAC = data["peak_DAC"]
    df['alt'] = data['alt']
    df['pwv'] = data['pwv']
    df['hour_utc'] = data['hour_utc']
    df['ctime'] = data['ctime']
    df.index = pd.to_datetime(df.ctime, unit='s')


if params.has_key('solid_angle'):
    tods, omega = np.loadtxt(
        'SOLID_ANGLES/%s.txt' %tag_calib, dtype=str, usecols=[0,1]).T
    df2 = pd.DataFrame()
    df2['tods'] = np.asarray(tods, dtype=str)
    df2['omega'] = np.asarray(omega, dtype=float)
    df = df.merge(df2, how='outer', on='tods')
    df.index = pd.to_datetime(df.ctime, unit='s')


if params.has_key('recal'):
    f = open(params.get('recal'), 'r')
    lines = f.readlines()
    f.close()
    tods_recal = np.array([l.split()[0] for l in lines], dtype='str')
    recal = np.array([l.split()[1] for l in lines], dtype='float')
    for t in zip(tods):
        cal[:,tods_recal==t] *= recal[tods_recal==t][0]

# Reject day TODs, loading > 2.7 and TODs from the commissioning phase before the season
df['loading'] = df.pwv / np.sin(df.alt)
sel_day = np.logical_and( df.hour_utc > 11, df.hour_utc<23)
print "Discard %i TODs during daytime" %(sel_day.sum())
print "Discard %i TODs for lack of peak measurement" %( peak_DAC.sum(axis=0) == 0 ).sum()
print "Discard %i TODs due to loading > %.1f" %((df.loading>loading_max).sum(), loading_max)
print "Discard %i TODs for lack of loading" %((df.loading<1).sum())
print "Discard %i TODs before the beginning of the season" %((df.ctime<ctime_start).sum())
idx = ( peak_DAC.sum(axis=0) != 0 ) * (df.loading>0) * (df.ctime > ctime_start) * (df.loading < loading_max)# * ~sel_day
df = df[idx]
print "Discard %i TODs in total, left with %i" %((~idx).sum(),df.tods.size)
cal = cal[:,idx]
peak_DAC = peak_DAC[:,idx]

# Compute the peak values
if params.get('season') == 's15':
    bso = []
    for ct in df.ctime:
        if ct < s15_split_time:
            bso.append(get_beam_solid_angle(params.get('array'),'s15a'))
        else:
            bso.append(get_beam_solid_angle(params.get('array'),'s15b'))
    beam_solid_angle = np.array(bso)
else:
    beam_solid_angle = get_beam_solid_angle(params.get('array'), params.get('season'))
    
print beam_solid_angle


T_diluted = T_URANUS * get_planet_solid_angle(df.ctime) / beam_solid_angle
if any([ar in array_name for ar in ['PA3', 'PA4', 'PA5', 'PA6']]):
    freq = np.float(array_name[4:])
    sel_freq = array_data['nom_freq'] == freq
    peak_masked = np.ma.masked_equal(peak_DAC*cal/T_diluted,0)
    peak_masked.mask[~sel_freq,:] = True
else:
    peak_masked = np.ma.masked_equal(peak_DAC*cal/T_diluted,0)[:]


# Reject TODs that have too few peak measurements
df['Ndets'] = (~peak_masked.mask).sum(axis=0)
live_min = params.get('live_min', 10)
sel = df.Ndets > live_min
print "Discard %i TODs for lack of good detectors" %((~sel).sum())
cal = cal[:,sel]
peak_DAC = peak_DAC[:,sel]
peak_masked = peak_masked[:,sel]
df = df[sel]


# Get the statistics (average and dispersion for each TOD) 
k = stats.mstats.kurtosistest(peak_masked, axis=0)
s = stats.mstats.skewtest(peak_masked, axis=0)
sel = (k.pvalue*s.pvalue) < params.get('gauss_test',1.)
# peak_mean = np.ma.median(peak_masked, axis=0)

p17, df['peak_mean'], p83 = stats.mstats.mquantiles(peak_masked,[0.17,0.5,0.83], axis=0)
df['peak_error'] = (p83-p17)/np.sqrt(df.Ndets)
#peak_error = np.ma.std(peak_masked, axis=0) / np.sqrt(Ndets)


# Fit
# Only points for which distribution was gaussian enough are considered
p = np.polyfit(df.loading[sel], np.log(df.peak_mean[sel]), deg=1, w=1./df.peak_error[sel])
x = np.linspace(0,5,100)
y = np.exp( np.polyval(p,x) )
p0 = np.exp( np.polyval(p,0) )
p1 = np.exp( np.polyval(p,0.5/np.sin(np.radians(50))) )
df['residual'] = df.peak_mean - np.exp( np.polyval(p,df.loading) )


# Make the plot
plt.ioff()
plt.figure(figsize=(18,9))
ax1 = plt.axes([.1,.1,.8,.8])
ax1.plot(x, y, 'b--')
ax1.errorbar(df.loading, df.peak_mean, df.peak_error, fmt=',', color='k', ecolor='k', elinewidth=2)
splot = ax1.scatter(df.loading, df.peak_mean, s=100,
#            c=df.hour_utc, vmin=0, vmax=24,
            c=df.Ndets,
            edgecolor='None', 
#            cmap='hsv')
            cmap='RdYlBu_r')
clb = plt.colorbar(splot)
clb.set_label('Ndets')
ax1.set_xlim((0,loading_max))
y_lim = params.get('y_lim', p0*2)
ax1.set_ylim((0,y_lim))
ax1.plot(df.loading[~sel], df.peak_mean[~sel], 'kx', ms=10, lw=5)
ax1.plot(df.loading[df.peak_mean>y_lim], np.ones_like(df.loading[df.peak_mean>y_lim])*y_lim, 'bx', ms=30)


ax1.text(2.5, 0.8*y_lim, '%s observations' %(df.loading.size), ha='right')
ax1.text(2.5, 0.77*y_lim, 'slope = %.4f / mm' %p[0], ha='right')
ax1.text(2.5, 0.74*y_lim, 'peak(loading=0) = %.5e' %p0, ha='right')
ax1.text(2.5, 0.71*y_lim, 'peak(pwv=.5,alt=50) = %.5e' %p1, ha='right')
ax1.text(2.5, 0.68*y_lim, 'scatter = %.3f' %(
        df.residual.std()/df.peak_mean.mean()), ha='right')

ax1.set_xlabel('PWV / sin(alt) [mm]')
ax1.set_ylabel('pW / [delta K CMB]')
ax1.set_title('%s' %tag_calib)

ax2 = plt.axes([.15,.15,.15,.15])
ax2.hist(df.residual/df.peak_mean.mean(), bins=np.linspace(-0.22,0.22,21))
ax2.set_xticks([-0.2,0.,0.2])
ax2.set_title('Residual')

plt.savefig('PLOTS/peak_vs_loading_%s_%s.png' %(tag_calib,array_name))


# apex_temperature = moby2.aux_data.apex.WeatherChannel('temperature')
# df['temperature'], dt = apex_temperature.get_nearest(df.ctime)
# df['residual_rel'] = df.residual / df.peak_mean.mean()
# df.hour_utc[df.hour_utc>20] -= 24
# df.plot('hour_utc', 'residual_rel', c='temperature', s=50, kind='scatter', cmap='RdYlBu_r')
# plt.ylim(-0.2,0.2)
# plt.xlim(-4,20)
# plt.xlabel('APEX temperature')
# plt.ylabel('Residual')
# plt.title(tag_calib)
# plt.savefig('PLOTS/residual_vs_hour_%s.png' %tag_calib)

df.to_csv('%s.csv' %tag_calib)
