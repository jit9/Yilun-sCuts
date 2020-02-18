import moby2
import pickle
import sys
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
import numpy as np

filename = '/home/yguan/cuts/s17/c11/pa4_f150/pa4_f150_s17_c11_v3_results.pickle'
array_name = 'ar4'
season = '2017'
calibrate = False

tag = filename.split('/')[-1][:-15]
data = pickle.load(open(filename,'r'))
if calibrate:
    data['rmsLive'] *= data['resp'] * data['ff'][:,np.newaxis]
    data['normLive'] *= data['resp'] * data['ff'][:,np.newaxis]
    data['MFELive'] *= data['resp'] * data['ff'][:,np.newaxis]
    data['DELive'] *= data['resp'] * data['ff'][:,np.newaxis]
    data['jumpLive'] *= data['resp'] * data['ff'][:,np.newaxis]

array_data = moby2.scripting.get_array_data({'array_name':array_name,'season':season})
freqs = np.unique(array_data['nom_freq'])
freqs = freqs[freqs!=0]
# sel_freqs = [array_data['nom_freq'] == f for f in freqs]
sel_freqs = [150]

keys = [k for k in list(data.keys()) if 'Live' in k]

plt.ioff()
plt.figure(figsize=(20,12))
# Gain
plt.subplot(341)
d = data['gainLive']
p1,p5,p95,p99 = scoreatpercentile(d[d!=0], [1,5,95,99])
for f, sel in zip(freqs,sel_freqs):
    d = data['gainLive'][sel]
    h, b, e = plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/5),np.log10(p95*5),100), alpha=0.5, label='%i' %f)
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
    h,b,e=plt.hist(d[np.isfinite(d)], bins=np.linspace(0.8,1,100), alpha=0.5, label='%i' %f)
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
    h,b,e=plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/5),np.log10(p95*5),100), alpha=0.5, label='%i' %f)
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
    h,b,e=plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/10),np.log10(p95*10),100), alpha=0.5, label='%i' %f)
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
    h,b,e=plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/20),np.log10(p95*20),100), alpha=0.5, label='%i' %f)
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
    h,b,e=plt.hist(d[np.isfinite(d)], bins=np.logspace(np.log10(p5/20),np.log10(p95*20),100), alpha=0.5, label='%i' %f)
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
    h,b,e=plt.hist(d[np.isfinite(d)*(d>0)], bins=np.logspace(np.log10(p1/3),np.log10(p99*3),100), alpha=0.5, label='%i' %f)
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
plt.savefig('PLOTS/seasoncrit_hist_%s_%s.png' %(tag, postfix))
plt.close('all')
