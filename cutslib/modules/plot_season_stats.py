"""This module uses the SeasonStats class to produce various plots
for debugging purposes to be included in the report

"""

import os.path as op, numpy as np
import matplotlib.pyplot as plt

from cutslib import SeasonStats


class Module:
    def __init__(self, config):
        self.use_theta2 = config.getboolean('use_theta2', True)
        self.calibrate = config.getboolean('calibrate', False)
        self.tri = config.getboolean('tri',False)

    def run(self, p):
        # create SeasonStats object from the given tag
        ss = SeasonStats(tag=p.tag, use_theta2=self.use_theta2,
                         calibrate=self.calibrate)
        # produce resp hist
        ndet_wresp = np.sum(ss.stats['resp_sel']*ss.stats['tes_sel'][:,None],axis=0)
        plt.plot(np.sort(ndet_wresp),'k-')
        plt.xlabel('TOD')
        plt.ylabel('# of dets with resp')
        # plt.twinx()
        # plt.plot(np.sort(ndet_wresp)/np.sum(ss.tes_sel),'k-')
        # plt.ylabel("fraction to tes dets")
        plt.title("Dets with valid bias-step measurements")
        outfile = op.join(p.o.cal.resp, 'resp_frac.png')
        print(f"Saving plot: {outfile}")
        plt.savefig(outfile)
        plt.close()

        # crit histogram with thresholds
        outfile = op.join(p.o.patho.root, 'hist_with_crits.png')
        print(f"Saving plot: {outfile}")
        ss.hist()
        plt.savefig(outfile)
        plt.close()

        # pwv histogram
        outfile = op.join(p.o.root, 'pwv_hist.png')
        print(f"Saving plot: {outfile}")
        ss.hist_pwv()
        plt.savefig(outfile)
        plt.close()

        # sel plot
        outfile = op.join(p.o.patho.root, 'view_sel.png')
        print(f"Saving plot: {outfile}")
        ss.view_sel(ss.stats['sel'])
        plt.savefig(outfile)
        plt.close()

        # view cuts
        outfile = op.join(p.o.patho.root,'view_cuts.png')
        print(f"Saving plot: {outfile}")
        ss.view_cuts()
        plt.savefig(outfile)
        plt.close()

        # triangular plot
        # this is slow to produce so only do it when specified
        if self.tri:
            outfile = op.join(p.o.patho.root, 'tri.png')
            print(f"Saving plot: {outfile}")
            ss.tri(figsize=(30,30), filename=outfile)
            plt.close()

        print('Done')
