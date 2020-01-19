"""This script aims to generate a report of the plots"""
TEMPLATE="""#+TITLE: ={tag}=
* Parameters
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
"""

import os, glob
import os.path as op
import moby2

def init(config):
    pass

def run(p):
    # load cut parameters
    par = moby2.util.MobyDict.from_file(p.i.cutParam)

    ################
    # cuts summary #
    ################
    live_par = par.get_deep(('pathologyParams', 'liveSelParams'))
    cuts_summary =  "|--------------+-------------------+-------|\n"
    cuts_summary += "| type         | crit              | apply |\n"
    cuts_summary += "|--------------+-------------------+-------|\n"
    for k, v in live_par.items():
        cuts_summary += "| %s | %s | %s |\n" % (k,
                                                v['absCrit'],
                                                v['apply'])
    cuts_summary += "|--------------+-------------------+-------|\n"

    #######################
    # cut threshold plots #
    #######################
    cuts_threshold = glob.glob(op.join(p.o.patho.root, "seasoncrit_hist*.png"))[0]
    ##################
    # killed_by_plot #
    ##################
    killed_by_plot = glob.glob(op.join(p.o.patho.root,"killedbyplot_*.png"))[0]
    ##################
    # flatfield plot #
    ##################
    flatfield_in = op.join(p.o.ff,"ff_%s_cal_input.png" % p.tag)
    flatfield_out = op.join(p.o.ff,"ff_%s_cal_output.png" % p.tag)

    #######################
    # generate org report #
    #######################
    report = TEMPLATE.format(
        tag=p.tag,
        cuts_summary=cuts_summary,
        cuts_threshold=cuts_threshold,
        killed_by_plot=killed_by_plot,
        flatfield_in=flatfield_in,
        flatfield_out=flatfield_out,
    )
    # save org report
    outfile = op.join(p.e.root, "report.org")
    print("Writing report: %s" % outfile)
    with open(outfile, "w") as f:
        f.write(report)
    # generate pdf report by converting org to pdf
    cmd="emacs --batch %s -f org-latex-export-to-pdf" % outfile
    print(cmd)
    os.system(cmd)
