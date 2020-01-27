"""This script aims to generate a report of the plots, it depends on
emacs and latex dependencies.

"""

TEMPLATE="""#+TITLE: ={tag}=
* Changelog
{changelog}
* Run
- TOD list: ={source_scans}=
  |-------+-----------+--------|
  | Total | Processed | ld>100 |
  |-------+-----------+--------|
  | {ntod}| {nproc}   | {ngood}|
  |-------+-----------+--------|

  |----------+------------+--------------+--------------|
  | mean(ld) | total dets | frac         | total frac   |
  |----------+------------+--------------+--------------|
  | {mld:.1f}| {ndets}    |{ldfrac:.1f}\%|{tfrac:.1f}\% |
  |----------+------------+--------------+--------------|

- Cuts parameters:
{cuts_summary}
* Statistics
- Uncut vs. loading
#+ATTR_LATEX: :width 16cm
[[{ld_vs_loading}]]
- Histogram:
#+ATTR_LATEX: :width 16cm
[[{cuts_threshold}]]

Note that the vertical lines on the histograms are for reference
only, they are not the same as the cuts applied on them.
- Killed by:
Number of detectors that passes each criteria

[[{killed_by_plot}]]
- Live fraction:
Fraction of the time that a detector is uncut

#+ATTR_LATEX: :width 10cm
[[{live_fraction}]]

- rms vs. input gain
[[{rms_vs_gain}]]
* Flatfield
#+BEGIN_center
#+ATTR_LaTeX: :height 0.45\\textwidth :center
[[{flatfield_in}]]
#+ATTR_LaTeX: :height 0.45\\textwidth :center
[[{flatfield_out}]]
#+END_center
The plot on the left (input) refers to the input flatfield for the
cuts pipeline. For this case it refers to the planet flatfield. On
the right is the output flatfield estimated from the atmosphere gain.
* Array plots
#+BEGIN_center
#+ATTR_LaTeX: :height 0.45\\textwidth :center
[[{gain_plot}]]
#+ATTR_LaTeX: :height 0.45\\textwidth :center
[[{corr_plot}]]
#+ATTR_LaTeX: :height 0.45\\textwidth :center
[[{rms_plot}]]
#+ATTR_LaTeX: :height 0.45\\textwidth :center
[[{norm_plot}]]
#+END_center
* Calibration
#+ATTR_LATEX: :width 16cm
[[{planet_cal}]]
"""

import os, glob
import os.path as op
import moby2
from moby2.util.database import TODList
from cutslib.pathologyReport import pathoReport

def init(config):
    global hits_map
    hits_map = config.get("hits_map", None)

def run(p):
    global hits_map
    # load cut parameters
    cutParam = moby2.util.MobyDict.from_file(p.i.cutParam)
    cutparam = moby2.util.MobyDict.from_file(p.i.cutparam)
    res = {'tag': p.tag}

    #############
    # changelog #
    #############
    res['changelog'] = cutParam.get('changelog', '')

    ############
    # tod list #
    ############
    source_scans = cutparam.get('source_scans')
    tod_list = TODList.from_file(source_scans)
    pr = pathoReport(filename=str(p.i.db))
    res['ntod'] = len(tod_list)
    res['nproc'] = len(pr.data)
    res['ngood'] = len(pr.data[pr.data.liveDets >= 100])
    res['source_scans'] = source_scans
    res['mld'] = pr.data.liveDets.mean(skipna=False)
    # get number of detectors
    ld_dict = cutParam.get_deep(('pathologyParams','detectorLists','live'))
    exclude_dict = cutParam.get_deep(('pathologyParams','detectorLists','exclude'))
    ld_data = moby2.util.MobyDict.from_file(ld_dict)
    exclude_data = moby2.util.MobyDict.from_file(exclude_dict)
    dets = [d for d in ld_data['det_uid'] if d not in exclude_data['det_uid']]
    res['ndets'] = len(dets)
    res['ldfrac'] = res['mld'] * 100. / res['ndets']
    res['tfrac'] = res['ldfrac'] * res['nproc'] / res['ntod']

    ################
    # cuts summary #
    ################
    live_par = cutParam.get_deep(('pathologyParams', 'liveSelParams'))
    cuts_summary =  "|--------------+-------------------+-------|\n"
    cuts_summary += "| type         | crit              | apply |\n"
    cuts_summary += "|--------------+-------------------+-------|\n"
    for k, v in live_par.items():
        cuts_summary += "| %s | %s | %s |\n" % (k,
                                                v['absCrit'],
                                                v['apply'])
    cuts_summary += "|--------------+-------------------+-------|\n"
    res['cuts_summary'] = cuts_summary

    #################
    # ld vs loading #
    #################
    res['ld_vs_loading'] = op.join(p.o.root, 'ld_vs_loading.png')

    #######################
    # cut threshold plots #
    #######################
    cuts_threshold = glob.glob(op.join(p.o.patho.root, "seasoncrit_hist*.png"))[0]
    res['cuts_threshold'] = cuts_threshold

    ##################
    # killed_by_plot #
    ##################
    killed_by_plot = glob.glob(op.join(p.o.patho.root,"killedbyplot.png"))[0]
    res['killed_by_plot'] = killed_by_plot

    ######################
    # live fraction plot #
    ######################
    res['live_fraction'] = glob.glob(op.join(p.o.root, "live_frac.png"))[0]

    ###############
    # rms vs gain #
    ###############
    res['rms_vs_gain'] = glob.glob(op.join(p.o.root, "rms_gain.png"))[0]

    ##################
    # flatfield plot #
    ##################
    res['flatfield_in'] = op.join(p.o.ff,"ff_%s_cal_input.png" % p.tag)
    res['flatfield_out'] = op.join(p.o.ff,"ff_%s_cal_output.png" % p.tag)

    ###############
    # array plots #
    ###############
    res['gain_plot'] = op.join(p.o.patho.array.root, 'gainLive_mean.png')
    res['corr_plot'] = op.join(p.o.patho.array.root, 'corrLive_mean.png')
    res['norm_plot'] = op.join(p.o.patho.array.root, 'normLive_mean.png')
    res['rms_plot']  = op.join(p.o.patho.array.root, 'rmsLive_mean.png')

    ###############
    # calibration #
    ###############
    res['planet_cal'] = glob.glob(op.join(p.o.cal.root, "peak_vs_loading*.png"))[0]

    ###################
    # generate report #
    ###################
    report = TEMPLATE.format(**res)
    # save org report
    outfile = op.join(p.e.root, "%s.org" % p.tag)
    print("Writing report: %s" % outfile)
    with open(outfile, "w") as f:
        f.write(report)
    # generate pdf report by converting org to pdf
    cmd="emacs --batch %s -f org-latex-export-to-pdf" % outfile
    print(cmd)
    os.system(cmd)
