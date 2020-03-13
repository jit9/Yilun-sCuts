"""This module creates the killedbyplot that shows how many detectors
pass each cut criteria

Options:
    limit_{l,h}: lower and upper limit for the plot"""

from cutslib.pathologyReport import pathoReport

class Module:
    def __init__(self, config):
        self.limit_l = config.getint("limit_l", 400)
        self.limit_h = config.getint("limit_h", 800)
        self.type = config.get("type", "violin")

    def run(self, p):
        limit_l = self.limit_l
        limit_h = self.limit_h
        type = self.type

        pr = pathoReport(filename=str(p.i.db))
        outfile = p.o.patho.root+"/killedbyplot.png"
        print("Saving plot: %s" % outfile)
        pr.killedbyplot(filename=outfile, dets_lim=[limit_l,limit_h], type=type)
