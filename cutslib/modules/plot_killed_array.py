"""This module creates the killedbyplot that shows how many detectors
pass each cut criteria

Options:
    limit_{l,h}: lower and upper limit for the plot"""

from cutslib.pathologyReport import pathoReport

def init(config):
    global limit_l, limit_h
    limit_l = config.getint("limit_l", 400)
    limit_h = config.getint("limit_h", 800)

def run(p):
    global limit_l, limit_h
    pr = pathoReport(filename=str(p.i.db))
    outfile = p.o.patho.root+"/killedbyplot.png"
    print("Saving plot: %s" % outfile)
    pr.killedbyplot(filename=outfile)
    outfile = p.o.patho.root+"/killedbyplot_%d_%d.png" % (limit_l, limit_h)
    print("Saving plot: %s" % outfile)
    pr.killedbyplot(filename=outfile, dets_lim=[limit_l,limit_h])
