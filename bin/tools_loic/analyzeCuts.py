"""This script produces some performance plots for the cuts results.
The plots iclude a season plot of live detectors and a killed-by
plots that's the violin plot for each cut crit. It also generate
frequency-space and time-space waterfall plots for 10 random TODs
with calibrations.
"""

from matplotlib import pyplot as plt
import os, sys, random, numpy as np

import moby2
from cutslib import visual
from cutslib.pathologyReport import pathoReport


params = moby2.util.MobyDict.from_file( sys.argv[1] )
pl = pathoReport( os.path.join(
    params.get('outdir'), params.get('report')+'.db') )

pl.addPWV()

if not os.path.exists(params.get('dirplot')):
    os.makedirs(params.get('dirplot'))


pl.seasonplot('liveDets', dets_lim=(0,pl.data.liveDets.max()*1.1),
              filename=os.path.join(
                  params.get('dirplot'), 'seasonplot_%s.png' %params.get('tag_out'))
          )
pl.killedbyplot(dets_lim=(0,pl.data.liveDets.max()*1.1),
                filename=os.path.join(
                    params.get('dirplot'), 'killedbyplot_%s.png' %params.get('tag_out'))
                )


depot_path = params.get('depot')
tag_cal = params.get('tag_cal')
tag_cuts = params.get('tag_out')
for obs in random.sample(pl.data.todName, 10):
    tod = moby2.scripting.get_tod({'filename':obs, 'repair_pointing':True})
    moby2.tod.remove_median(tod)
    moby2.tod.detrend_tod(tod)

    cuts = moby2.scripting.get_cuts(
        {'depot':depot_path,
         'tag':tag_cuts},
        tod=tod)
    tod.cuts = cuts
    moby2.tod.fill_cuts(tod)
    ld = cuts.get_uncut()

    cal = moby2.scripting.get_calibration(
        {'type':'depot_cal',
         'depot': depot_path,
         'tag': tag_cal},
        tod=tod)
    tod.data[cal.det_uid] *= cal.cal[:,np.newaxis]

    wt = visual.timeSpaceWaterfall(tod)
    wt.plot(selection=cuts.get_mask(), units='pW',
            title='%s - Time domain' %obs,
            show = False,
            filename = os.path.join(
                params.get('dirplot'), 'waterfall_time_%s_%s.png' %(obs,params.get('tag_out')) )
        )
    wf = visual.freqSpaceWaterfall(tod)
    # logx=False,fmin=8., fmax=50., nfreq=84)
    wf.plot(selection=cuts.get_mask(), units='pW',
            title='%s - Frequency domain' %obs,
            show = False,
            filename = os.path.join(
                params.get('dirplot'), 'waterfall_freq_%s_%s.png' %(obs,params.get('tag_out')) )
        )
    plt.close('all')
