"""This module creates the killedbyplot that shows how many detectors
pass each cut criteria"""

from pathologyReport import pathoReport

def init(config):
    pass

def run(p):
    pr = pathoReport(filename=str(p.i.db))
    outfile = p.o.patho.root+"/killedbyplot.png"
    print("Saving plot: %s" % outfile)
    pr.killedbyplot(filename=outfile)
    outfile = p.o.patho.root+"/killedbyplot_400_800.png"
    print("Saving plot: %s" % outfile)
    pr.killedbyplot(filename=outfile, dets_lim=[400,800])
