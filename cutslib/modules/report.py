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

  |-----------------------+-----------------------|
  | Type                  | ndets                 |
  |-----------------------+-----------------------|
  | f{freq}+tes           | {ndets_tes}           |
  | f{freq}+tes+ff        | {ndets_tes_ff}        |
  | f{freq}+tes+ff+stable | {ndets_tes_ff_stable} |
  |-----------------------+-----------------------|

The total number of detector in the table below also includes
the exclude list of unfavorable detectors based on previous
cuts run, so it is different than the number above.

  |-----------+------------+----------------+---------------|
  | mean(ld)  | total dets | frac           | total frac    |
  |-----------+------------+----------------+---------------|
  | {mld:.1f} | {ndets}    | {ldfrac:.1f}\% | {tfrac:.1f}\% |
  |-----------+------------+----------------+---------------|

The number of detectors with valid bias-step measurement is shown in
the figure below

[[{resp_frac}]]

Below is a histogram of the PWV of the TODs

[[{hist_pwv}]]

- Cuts parameters:
{cuts_summary}
* Statistics
- Uncut vs. loading
#+ATTR_LATEX: :width 16cm
[[{ld_vs_loading}]]

- live dets
#+ATTR_LATEX: :width 16cm
[[{view_sel}]]

- Histogram:
#+ATTR_LATEX: :width 16cm
[[{cuts_threshold}]]

Note that the vertical lines on the histograms are for reference
only, they are not the same as the cuts applied on them. The thresholds
are shown in the plots below

#+ATTR_LATEX: :width 16cm
[[{hist_cuts}]]

Triangular plot (if produced it will show up below)

#+ATTR_LATEX: :width 16cm
[[{tri}]]

- Killed by:
Number of detectors that passes each criteria

[[{killed_by_plot}]]

Number of dets that passes each criteria at various pwv

[[{view_cuts}]]

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
* Pathologies
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
import os.path as op, numpy as np
import moby2
from moby2.util.database import TODList
from cutslib.pathologyReport import pathoReport


class Module:
    def __init__(self, config):
        self.hits_map = config.get("hits_map", None)

    def run(self, p):
        hits_map = self.hits_map

        # load cut parameters
        cutParam = moby2.util.MobyDict.from_file(p.i.cutParam)
        cutparam = moby2.util.MobyDict.from_file(p.i.cutparam)
        res = {'tag': p.tag}

        #############
        # changelog #
        #############
        res['changelog'] = format_changelog(cutParam.get('changelog', ''))

        ############
        # tod list #
        ############
        source_scans = cutparam.get('source_scans')
        tod_list = TODList.from_file(source_scans)
        pr = pathoReport(filename=str(p.i.db))
        pr.drop_duplicates()
        res['ntod'] = len(tod_list)
        res['nproc'] = len(pr.data)
        res['ngood'] = len(pr.data[pr.data.liveDets >= 100])
        res['source_scans'] = source_scans
        res['mld'] = pr.data.liveDets.mean(skipna=False)

        #############
        # resp hist #
        #############
        res['resp_frac'] = op.join(p.o.cal.resp, "resp_frac.png")

        ############
        # pwv hist #
        ############
        res['hist_pwv'] = op.join(p.o.root, 'pwv_hist.png')

        ####################
        # detector numbers #
        ####################
        # get array data
        res['freq'] = "%s" % int(p.i.freq)
        array_data = moby2.scripting.get_array_data({
            'instrument': 'actpol',
            'array_name': p.i.ar,
            'season': p.i.season
        })
        # get number of tes detectors in the given freq
        dets_tes = array_data['det_uid'][(array_data['nom_freq']== p.i.freq) * (array_data['det_type'] == 'tes')]
        res['ndets_tes'] = len(dets_tes)
        # get the number of tes detectors with ff
        ff_file = cutParam.get_deep(("pathologyParams","calibration","flatfield"))
        ff_dict = moby2.util.MobyDict.from_file(ff_file)
        dets_ff = ff_dict.get("det_uid")
        res['ndets_tes_ff'] = len(set(dets_tes).intersection(set(dets_ff)))
        dets_stable = np.array(dets_ff)[np.array(ff_dict.get("stable"))]
        res['ndets_tes_ff_stable'] = len(set(dets_tes).intersection(set(dets_stable)))

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

        ############
        # view sel #
        ############
        res['view_sel'] = op.join(p.o.patho.root, 'view_sel.png')

        #############
        # view cuts #
        #############
        res['view_cuts'] = op.join(p.o.patho.root, 'view_cuts.png')

        #######################
        # cut threshold plots #
        #######################
        cuts_threshold = find_file(op.join(p.o.patho.root, "seasoncrit_hist*.png"))
        res['cuts_threshold'] = cuts_threshold

        ########################################
        # cut threshold plot (with thresholds) #
        ########################################
        res['hist_cuts'] = op.join(p.o.patho.root, 'hist_with_crits.png')

        ###################
        # triangular plot #
        ###################
        tri = op.join(p.o.patho.root, 'tri.png')
        if op.exists(tri):
            res['tri'] = tri

        ##################
        # killed_by_plot #
        ##################
        killed_by_plot = find_file(op.join(p.o.patho.root,"killedbyplot.png"))
        res['killed_by_plot'] = killed_by_plot

        ######################
        # live fraction plot #
        ######################
        res['live_fraction'] = find_file(op.join(p.o.root, "live_frac.png"))

        ###############
        # rms vs gain #
        ###############
        res['rms_vs_gain'] = find_file(op.join(p.o.root, "rms_gain.png"))

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
        res['planet_cal'] = find_file(op.join(p.o.cal.root, "peak_vs_loading*.png"))

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

# utility function to format changelog for better displaying
def format_changelog(changelog):
    lines = changelog.split('\n')
    output_lines = ""
    for i,l in enumerate(lines):
        if '- v' == l[:3]:
            output_lines += l + "\n"
        else:
            output_lines += "  " + l + "\n"
    return output_lines

def find_file(selector):
    results = glob.glob(selector)
    if len(results) > 0: return results[0]
    else: return ""
