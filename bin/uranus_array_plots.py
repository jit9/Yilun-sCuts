import os
import moby2
import sys
import numpy as np
import pickle
import pandas as pd
from moby2.analysis.tod_ana import visual as v

# Initialize parameters
params = moby2.util.MobyDict.from_file(sys.argv[1])
fb = moby2.scripting.get_filebase()
depot = moby2.util.Depot(params.get('depot'))
if 'tag' in params:
    tag_calib = params.get('tag')
    tag_cuts = params.get('tag')
    tag_selected = params.get('tag')
else:
    tag_calib = params.get('tag_calib')
    tag_cuts = params.get('tag_cuts')
    tag_selected = params.get('tag_selected')
if 'flatfield' in params:
    FF = moby2.util.MobyDict.from_file(params.get('flatfield'))
tods, flag = np.loadtxt(
    os.path.join(params.get('depot'), 'SelectedTODs/%s/selectedTODs_uranus.txt' %tag_selected), dtype=str, usecols=[0,5]).T

df = pd.DataFrame()
df['tods'] = np.asarray( tods, dtype=str )
df['flag'] = np.asarray( flag, dtype=int )
df = df[df.flag==2]


f = open("PICKLE/%s.pickle"%tag_calib,"r")
data = pickle.load(f)
f.close()
cal = data["cal"]
peak_DAC = data["peak_DAC"]
df['alt'] = data['alt']
df['pwv'] = data['pwv']
df['hour_utc'] = data['hour_utc']
df['ctime'] = data['ctime']
df.index = pd.to_datetime(df.ctime, unit='s')

print("Start with %i selected TODs" %df.tods.size) 
path = params.get('planet_path')
array_name = params.get('array')
array_data =  moby2.scripting.products.get_array_data(
    {'instrument':'actpol','array_name':'ar%s' %array_name[2],'season':'20%s' %params.get('season')[-2:]})
ctime_start = params.get('ctime_start', 0)
loading_max = params.get('loading_max', 10)
alma_pwv = moby2.aux_data.alma.WeatherChannel('pwv')
apex_pwv = moby2.aux_data.apex.WeatherChannel('radiometer')


# Reject day TODs, loading > 2.7 and TODs from the commissioning phase before the season
df['loading'] = df.pwv / np.sin(df.alt)
sel_day = np.logical_and( df.hour_utc > 11, df.hour_utc<23)
print("Discard %i TODs during daytime" %(sel_day.sum()))
print("Discard %i TODs for lack of peak measurement" %( peak_DAC.sum(axis=0) == 0 ).sum())
print("Discard %i TODs due to loading > %.1f" %((df.loading>loading_max).sum(), loading_max))
print("Discard %i TODs for lack of loading" %((df.loading<1).sum()))
print("Discard %i TODs before the beginning of the season" %((df.ctime<ctime_start).sum()))
idx = ( peak_DAC.sum(axis=0) != 0 ) * (df.loading>0) * (df.ctime > ctime_start) * (df.loading < loading_max)# * ~sel_day
df = df[idx]
print("Discard %i TODs in total, left with %i" %((~idx).sum(),df.tods.size))
cal = cal[:,idx]
peak_DAC = peak_DAC[:,idx]

sel = array_data['nom_freq'] == params.get('freq', 150)
peak_masked = np.ma.masked_equal(peak_DAC*cal,0)[:]
for i in range(peak_masked.shape[1]):
    peak = peak_masked[:,i]
    v.array_plots(
        peak.data[sel],
        np.arange(peak_masked.shape[0])[sel],
        array='ar4', season=2017,
        pmin=2e-5, pmax=3e-2,
        title = '{} - {}'.format(df.tods.iloc[i], tag_calib),
        display='save',
        save_name='PLOTS/uranus_array_plots/{}_{}.png'.format(df.tods.iloc[i], tag_calib))
