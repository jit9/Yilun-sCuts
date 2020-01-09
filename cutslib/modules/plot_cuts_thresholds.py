"""This module plots the histograms for each selCriteria.

Options:
    calibrate: True if rms, norm, MFE, DE, jump is to be calibrated
    {crit}_adj: adjust the histogram range, larger value means larger range
    corr_l: lower plot range for corr
"""

import moby2
import cPickle
import sys
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
import numpy as np

def init(config):
    global calibrate, gain_adj, corr_l, rms_adj, norm_adj, de_adj, mfe_adj, jump_adj
    calibrate = config.getboolean("calibrate", False)
    gain_adj = config.getfloat("gain_adj", 5)
    corr_l = config.getfloat("corr_l", 0.8)
    rms_adj = config.getfloat("rms_adj", 5)
    norm_adj = config.getfloat("norm_adj", 10)
    de_adj = config.getfloat("de_adj", 20)
    mfe_adj = config.getfloat("mfe_adj", 20)
    jump_adj = config.getfloat("jump_adj", 3)

def run(proj):
    global calibrate, gain_adj, corr_l, rms_adj, norm_adj, de_adj, mfe_adj, jump_adj
    filename = proj.i.pickle_file
    array_name = proj.i.ar
    season = proj.i.season
    tag = proj.tag
    data = cPickle.load(open(filename,'r'))
    if calibrate:
        data['rmsLive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['normLive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['MFELive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['DELive'] *= data['resp'] * data['ff'][:,np.newaxis]
        data['jumpLive'] *= data['resp'] * data['ff'][:,np.newaxis]

    array_data = moby2.scripting.get_array_data({'array_name':array_name,'season':season})
    # freqs = np.unique(array_data['nom_freq'])
    # freqs = freqs[freqs!=0]
    freqs = [proj.i.freq]
    sel_freqs = [array_data['nom_freq'] == f for f in freqs]

    keys = [k for k in data.keys() if 'Live' in k]

    plt.ioff()
    plt.figure(figsize=(20,12))
    # Gain
    plt.subplot(341)
    d = data['gainLive']
    p1,p5,p95,p99 = scoreatpercentile(d[d!=0], [1,5,95,99])
    for f, sel in zip(freqs,sel_freqs):
        d = data['gainLive'][sel]
        adj = gain_adj
        h, b, e = plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/adj),np.log10(p95*adj),100), alpha=0.5, label='%i' %f)
    plt.axvline(p5, color='g')
    plt.text(p5, h.max()/2, 'p5 = %.2f'%p5, color='g', ha='right')
    plt.axvline(p95, color='g')
    plt.text(p95, h.max()/2, 'p95 = %.2f'%p95, color='g')
    plt.xscale('log')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('Gain')

    # Correlation
    plt.subplot(342)
    d = data['corrLive']
    p1,p5,p10 = scoreatpercentile(d[d!=0], [1,5,10])
    for f, sel in zip(freqs,sel_freqs):
        d = data['corrLive'][sel]
        h,b,e=plt.hist(d[np.isfinite(d)], bins=np.linspace(corr_l,1,100), alpha=0.5, label='%i' %f)
    plt.text(p10, h.max()/2, 'p10 = %.4f'%p10, color='g', ha='right')
    plt.axvline(p10, color='g')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('Correlation')

    # Dark ratio
    #plt.subplot(343)
    #d = data['darkRatioLive']
    #p1,p5,p10 = scoreatpercentile(d[d!=0], [1,5,10])
    #for f, sel in zip(freqs,sel_freqs):
    #    d = data['darkRatioLive'][sel]
    #    h,b,e=plt.hist(d[np.isfinite(d)], bins=np.linspace(0.01,1,100), alpha=0.5, label='%i' %f)
    #plt.text(p10, h.max()/2, 'p10 = %.4f'%p10, color='g', ha='right')
    #plt.axvline(p10, color='g')
    #plt.legend(loc='best',frameon=False)
    #plt.yticks(visible=False)
    #plt.title('Dark Fraction')

    # RMS
    plt.subplot(344)
    d = data['rmsLive']
    p1,p5,p95,p99 = scoreatpercentile(d[d!=0], [1,5,95,99])
    for f, sel in zip(freqs,sel_freqs):
        d = data['rmsLive'][sel]
        adj = rms_adj
        h,b,e=plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/adj),np.log10(p95*adj),100), alpha=0.5, label='%i' %f)
    plt.axvline(p5, color='g')
    plt.text(p5, h.max()/2, 'p5 = %.2e'%p1, color='g', ha='right')
    plt.axvline(p95, color='g')
    plt.text(p95, h.max()/2, 'p95 = %.2e'%p95, color='g')
    plt.xscale('log')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('RMS')

    # Norm
    plt.subplot(345)
    d = data['normLive']
    p1,p5,p95,p99 = scoreatpercentile(d[d!=0], [1,5,95,99])
    for f, sel in zip(freqs,sel_freqs):
        d = data['normLive'][sel]
        adj = norm_adj
        h,b,e=plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/adj),np.log10(p95*adj),100), alpha=0.5, label='%i' %f)
    plt.axvline(p5, color='g')
    plt.text(p5, h.max()/2, 'p5 = %.2e'%p1, color='g', ha='right')
    plt.axvline(p95, color='g')
    plt.text(p95, h.max()/2, 'p95 = %.2e'%p95, color='g')
    plt.xscale('log')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('Norm')

    # Drift Error
    plt.subplot(346)
    d = data['DELive']
    p1,p5,p95,p99 = scoreatpercentile(d[d!=0], [1,5,95,99])
    for f, sel in zip(freqs,sel_freqs):
        d = data['DELive'][sel]
        adj = de_adj
        h,b,e=plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/adj),np.log10(p95*adj),100), alpha=0.5, label='%i' %f)
    plt.axvline(p5, color='g')
    plt.text(p5, h.max()/2, 'p5 = %.2e'%p1, color='g', ha='right')
    plt.axvline(p95, color='g')
    plt.text(p95, h.max()/2, 'p95 = %.2e'%p95, color='g')
    plt.xscale('log')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('Drift Error')

    # Mid-Frequency Error
    plt.subplot(347)
    d = data['MFELive']
    p1,p5,p95,p99 = scoreatpercentile(d[d!=0], [1,5,95,99])
    for f, sel in zip(freqs,sel_freqs):
        d = data['MFELive'][sel]
        adj = mfe_adj
        h,b,e=plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/adj),np.log10(p95*adj),100), alpha=0.5, label='%i' %f)
    plt.axvline(p5, color='g')
    plt.text(p5, h.max()/2, 'p5 = %.2e'%p1, color='g', ha='right')
    plt.axvline(p95, color='g')
    plt.text(p95, h.max()/2, 'p95 = %.2e'%p95, color='g')
    plt.xscale('log')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('Mid-Frequency Error')

    # Kurtosis
    plt.subplot(349)
    d = data['kurtLive']
    p1,p5,p95,p99 = scoreatpercentile(d[d!=0], [1,5,95,99])
    for f, sel in zip(freqs,sel_freqs):
        d = data['kurtLive'][sel]
        h,b,e=plt.hist(d[np.isfinite(d)*(d!=0)], bins=np.linspace(p1*5,p95*3,100), alpha=0.5, label='%i' %f)
    plt.axvline(p1, color='g')
    plt.text(p1, h.max()/2, 'p1 = %.2f'%p1, color='g', ha='right')
    plt.axvline(p95, color='g')
    plt.text(p95, h.max()/2, 'p95 = %.2f'%p95, color='g')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('Kurtosis')

    # Skewness
    plt.subplot(3,4,10)
    d = data['skewLive']
    p1,p5,p95,p99 = scoreatpercentile(d[d!=0], [1,5,95,99])
    for f, sel in zip(freqs,sel_freqs):
        d = data['skewLive'][sel]
        h,b,e=plt.hist(d[np.isfinite(d)*(d!=0)], bins=np.linspace(p1*3,p99*3,100), alpha=0.5, label='%i' %f)
    plt.axvline(p1, color='g')
    plt.text(p1, h.max()/2, 'p1 = %.2f'%p1, color='g', ha='right')
    plt.axvline(p99, color='g')
    plt.text(p99, h.max()/2, 'p99 = %.2f'%p99, color='g')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('Skewness')

    # Jump
    plt.subplot(3,4,11)
    d = data['jumpLive']
    p1,p5,p95,p99 = scoreatpercentile(d[d>0], [1,5,95,99])
    for f, sel in zip(freqs,sel_freqs):
        d = data['jumpLive'][sel]
        adj = jump_adj
        h,b,e=plt.hist(d[np.isfinite(d)*(d>0)], bins=np.logspace(np.log10(p1/adj),np.log10(p99*adj),100), alpha=0.5, label='%i' %f)
    plt.xscale('log')
    plt.axvline(p1, color='g')
    plt.text(p1, h.max()/2, 'p1 = %.2e'%p1, color='g', ha='right')
    plt.axvline(p99, color='g')
    plt.text(p99, h.max()/2, 'p99 = %.2e'%p99, color='g')
    plt.legend(loc='best',frameon=False)
    plt.yticks(visible=False)
    plt.title('Jump')

    if calibrate:
        postfix="calibrated"
    else:
        postfix = "uncalibrated"
    plt.figtext(0.5,0.95,tag+' '+postfix,ha='center',va='center',fontsize='xx-large')
    plt.savefig(proj.o.patho.root+'/seasoncrit_hist_%s_%s.png' %(tag, postfix))
    plt.close('all')
