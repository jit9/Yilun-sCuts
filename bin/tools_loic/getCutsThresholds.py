import matplotlib
matplotlib.use('Agg')
import numpy as np
import pickle
import sys
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

filename = sys.argv[1]

data = pickle.load(open(filename, 'r'))
liveKeys = [k for k in list(data.keys()) if 'Live' in k]

plt.ioff()
plt.figure()
# for k in liveKeys:
# corrLive
k = 'corrLive'
d = data[k]
d = d[d!=0]
p = scoreatpercentile(d, [1, 3, 5, 10])
plt.hist(d, range=[0.7,1],bins=50 )
plt.axvline(p[0], color='r')
plt.axvline(p[1], color='r')
plt.axvline(p[2], color='r')
plt.text(p[0], 50000, '%.3f' %p[0])
plt.text(p[1], 50000, '%.3f' %p[1])
plt.text(p[2], 50000, '%.3f' %p[2])
plt.title(k)
plt.savefig('plots/%s_hist_%s.png' %(filename[:-15],k))
plt.clf()

# gainLive
k = 'gainLive'
d = np.asarray(data[k]) * np.asarray(data['ff'])[:,np.newaxis]
d = d[d!=0]
p = scoreatpercentile(d, [1., 5, 95., 99])
plt.hist(d, bins=np.logspace(-2,2,50) )
plt.xscale('log')
plt.axvline(p[0], color='r')
plt.axvline(p[1], color='r')
plt.axvline(p[-1], color='r')
plt.axvline(p[-2], color='r')
plt.text(p[0], 50000, '%.3f' %p[0])
plt.text(p[1], 50000, '%.3f' %p[1])
plt.text(p[-1], 50000, '%.3f' %p[-1])
plt.text(p[-2], 50000, '%.3f' %p[-2])
plt.title(k)
plt.savefig('plots/%s_hist_%s.png' %(filename[:-15],k))
plt.clf()

# rmsLive
k = 'rmsLive'
d = np.abs( np.asarray(data[k])*np.asarray(data['cal'])*np.asarray(data['ff'])[:,np.newaxis] )[np.asarray(data['sel'],dtype=bool)]
d = d[d!=0]
p = scoreatpercentile(d, [0.01, 1., 99., 99.99])
plt.hist(d, bins=np.logspace(np.log10(p[0]), np.log10(p[-1]), 50) )
plt.xscale('log')
plt.axvline(p[1], color='r')
plt.axvline(p[2], color='r')
plt.text(p[1], 50000, '%.3f' %p[1])
plt.text(p[2], 50000, '%.3f' %p[2])
plt.title(k)
plt.savefig('plots/%s_hist_%s.png' %(filename[:-15],k))
plt.clf()

