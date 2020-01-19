"""This script aims to generate a report of the plots, it depends on
emacs and latex dependencies.

"""

TEMPLATE="""#+TITLE: ={tag}=
* Run
- TOD list: ={source_scans}=
  |-------+-----------+--------|
  | Total | Processed | ld>100 |
  |-------+-----------+--------|
  | {ntod}| {nproc}   | {ngood}|
  |-------+-----------+--------|
- Hits map:
#+ATTR_LATEX: :width 10cm
[[{hits_map}]]
- Cuts parameters:
{cuts_summary}
* Statistics
- Histogram:
#+ATTR_LATEX: :width 15cm
[[{cuts_threshold}]]
- Killed by:
[[{killed_by_plot}]]
* Flatfield
#+BEGIN_center
#+ATTR_LaTeX: :height 0.45\\textwidth :center
[[{flatfield_in}]]
#+ATTR_LaTeX: :height 0.45\\textwidth :center
[[{flatfield_out}]]
#+END_center
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
[[{planet_cal}]]
"""

import os, glob
import os.path as op
import moby2
from moby2.util.database import TODList
from cutslib.pathologyReport import pathoReport

def init(config):
    pass

def run(p):
    # load cut parameters
    cutParam = moby2.util.MobyDict.from_file(p.i.cutParam)
    cutparam = moby2.util.MobyDict.from_file(p.i.cutparam)
    res = {'tag': p.tag}

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

    ############
    # hits map #
    ############
    res['hits_map'] = op.join(p.o.root, "hits.png")

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

    #######################
    # cut threshold plots #
    #######################
    cuts_threshold = glob.glob(op.join(p.o.patho.root, "seasoncrit_hist*.png"))[0]
    res['cuts_threshold'] = cuts_threshold

    ##################
    # killed_by_plot #
    ##################
    killed_by_plot = glob.glob(op.join(p.o.patho.root,"killedbyplot_*.png"))[0]
    res['killed_by_plot'] = killed_by_plot

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
