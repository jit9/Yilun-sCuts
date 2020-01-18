"""This script aims to generate a report of the plots"""
TEMPLATE="""#+TITLE: Cuts report (={tag}=)

* Statistics
[[{cuts_threshold}]]

* Flatfield
[[{flatfield}]]
"""

import os, glob
import os.path as op

def init(config):
    pass

def run(p):
    cuts_threshold = glob.glob(op.join(p.o.patho.root, "seasoncrit_hist*.png"))[0]
    flatfield = op.join(p.o.ff,"ff_%s_cal.png" % p.tag)
    report = TEMPLATE.format(
        tag=p.tag,
        cuts_threshold=cuts_threshold,
        flatfield=flatfield
    )
    # generate org report
    outfile = op.join(p.e.root, "report.org")
    print("Writing report: %s" % outfile)
    with open(outfile, "w") as f:
        f.write(report)
    # generate pdf report
    cmd="emacs --batch %s -f org-latex-export-to-pdf" % outfile
    print(cmd)
    os.system(cmd)
