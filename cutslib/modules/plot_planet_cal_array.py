import os
import moby2
import sys
import numpy as np
import cPickle
import pandas as pd
from moby2.analysis.tod_ana import visual as v

def init(config):
    global ctime_start, loading_max, path
    ctime_start = config.getint('ctime_start', 0)
    loading_max = config.getfloat('loading_max', 10)    
    path = config.get('planet_path')

def run(proj):
    global ctime_start, loading_max, path
    # Initialize parameters
    params = moby2.util.MobyDict.from_file(proj.i.cutparam)
    fb = moby2.scripting.get_filebase()
    depot = moby2.util.Depot(params.get('depot'))
    tag_calib = proj.tag
    tag_cuts = proj.tag
    tag_selected = proj.tag

    if params.has_key('flatfield'):
        FF = moby2.util.MobyDict.from_file(params.get('flatfield'))
    tods, flag = np.loadtxt(
        os.path.join(params.get('depot'), 'SelectedTODs/%s/selectedTODs_uranus.txt' %tag_selected), dtype=str, usecols=[0,5]).T

    df = pd.DataFrame()
    df['tods'] = np.asarray( tods, dtype=str )
    df['flag'] = np.asarray( flag, dtype=int )
    df = df[df.flag==2]

    f = open(proj.o.cal.root+"/%s.pickle"%tag_calib,"r")
    data = cPickle.load(f)
    f.close()
    cal = data["cal"]
    peak_DAC = data["peak_DAC"]
    df['alt'] = data['alt']
    df['pwv'] = data['pwv']
    df['hour_utc'] = data['hour_utc']
    df['ctime'] = data['ctime']
    df.index = pd.to_datetime(df.ctime, unit='s')

    print "Start with %i selected TODs" %df.tods.size
    season = 's'+proj.i.season[-2:]
    array = 'PA'+proj.i.ar[-1]
    if int(proj.i.ar[-1]) >= 3:
        array_name = array+"_"+str(proj.i.freq)
    else:
        array_name = array    
    array_data =  moby2.scripting.products.get_array_data(
        {'instrument':'actpol','array_name':proj.i.ar,'season':proj.i.season})
    alma_pwv = moby2.aux_data.alma.WeatherChannel('pwv')
    apex_pwv = moby2.aux_data.apex.WeatherChannel('radiometer')

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

    sel = array_data['nom_freq'] == int(proj.i.freq)
    peak_masked = np.ma.masked_equal(peak_DAC*cal,0)[:]
    for i in xrange(peak_masked.shape[1]):
        outfile = proj.o.cal.array+'/{}_{}.png'.format(df.tods.iloc[i], tag_calib)
        print("Saving plot: %s" % outfile)
        peak = peak_masked[:,i]
        v.array_plots(
            peak.data[sel],
            np.arange(peak_masked.shape[0])[sel],
            array=proj.i.ar, season=int(proj.i.season),
            pmin=1e-15,# pmax=3e-2,
            title = '{} - {}'.format(df.tods.iloc[i], tag_calib),
            display='save',
            save_name=outfile)

    sa_dict = moby2.util.MobyDict.from_file("/home/yguan/work/cuts_analysis/data/pa4_s16_f150_solid_angles.dict")
    sa_det_uid = np.array(sa_dict['det_uid']).astype(int)
    sa_cal = np.array(sa_dict['omega']) / 218.2

    # debug plot: uncalibrated peak
    sa_sel = np.zeros(len(sel)).astype(bool)
    sa_sel[sa_det_uid] = True
    sel = np.logical_and(sa_sel, sel)
    peak_masked = np.ma.masked_equal(peak_DAC,0)[:]
    for i in range(len(df.index)):
        outfile = proj.o.cal.peak+"/{}_peak_corrected.png".format(df.tods.iloc[i])
        peak = peak_masked[:,i]
        peak[sa_sel] *= sa_cal
        print("Saving plot: %s" % outfile)
        v.array_plots(peak.data[sel], np.arange(peak_masked.shape[0])[sel],
                      season=proj.i.season, array=proj.i.ar,
                      title=df.tods.iloc[i] + " peak_DAC (uncalibrated) - %s"
                      % tag_calib, pmin=1e-10, display='save',
                      save_name=outfile)
        
