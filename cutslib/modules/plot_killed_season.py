"""This module creates the killedbyplot as a time series 
that shows how many detectors pass each cut criteria as a function
of time."""

from pathologyReport import pathoReport

def init(config):
    global targets, add_pwv
    targets = config.get("targets", None)
    add_pwv = config.getboolean("add_pwv", False)

def run(p):
    global targets, add_pwv
    pr = pathoReport(filename=str(p.i.db))
    if add_pwv:
        pr.addPWV()

    # if targets are not specified, all are calculated
    if not targets:
        targets = ['corrLive', 'rmsLive', 'kurtLive', 'skewLive',
                   'normLive', 'darkRatioLive', 'MFELive',
                   'gainLive', 'DELive', 'jumpLive']
    
    for target in targets:
        outfile = p.o.patho.season.root+"/%s.png" % target
        print("Saving plot: %s" % outfile)
        try:
            pr.seasonplot(crit=target, filename=outfile)
        except KeyError:
            print("Key %s not found!" % target)
