"""This module aims to anwser the question that how many dets are cut
versus pwv. The difference between this and the plot_live_fraction
module is that i aim to provide more information than the mean / median,
so it can be seen as an extended version of plot_live_fraction. It also
accounts for the missing tods which was not considerred in the
plot_live_fraction module.

Note that the results from this module depends on the results from the
collect_sel module.

"""
import os.path as op, numpy as np, pickle
import matplotlib.pyplot as plt

from cutslib.visual import set_plotstyle


####################
# utility function #
####################

def plot_smoothed(x, y, ax=None, window=10, mode='valid', **kwargs):
    """plot a moving average of a certain window"""
    kernel = np.ones(window)/window
    x_smt = np.convolve(x, kernel, mode=mode)
    y_smt = np.convolve(y, kernel, mode=mode)
    ax.plot(x_smt, y_smt, **kwargs)
    return ax


###############
# main module #
###############

class Module:
    def __init__(self, config):
        self.window = config.getint("window", 10)
        self.show_baseline = config.getboolean("show_baseline", True)

    def run(self, p):
        window = self.window
        show_baseline = self.show_baseline

        # first check whether we have the required file available
        sel_file = op.join(p.o.patho.root, "sel.pickle")
        if not op.exists(sel_file):
            raise ValueError("Missing required sel.pickle file! "\
                             "Have you run collect_sel module?")
        with open(sel_file, "rb") as f:
            data = pickle.load(f)
        # The problem is a bit tricky because some of the lists do not
        # have the same size. We will call the total number of TODs as
        # n_tot and the total number of tods with valid pathology
        # object as n_patho, these two numbers usually differ by
        # 10%. The lists that have n_tod entries include: pwv, tod,
        # tod_sel, while most of the other lists will have a dimension
        # with n_patho. The exceptions are calibration related sels
        # which are the same for all tods so they have shape of
        # n_dets. Let me get these numbers defined below
        n_tot = len(data['tod'])
        n_patho = data['gainLive'].shape[0]
        n_dets = len(data['ff_sel'])
        print(f"n_dets: {n_dets}; n_tot: {n_tot}; n_patho: {n_patho}")
        # first we try to find the num of dets with valid resp verses
        # pwv, for the tods without valid pathology object (n_patho),
        # we will fill them with 0. These tods can be found by tod_sel
        tod_sel = data['tod_sel']  # n_tot entries
        baseline = data['tes_sel']
        pwv = data['pwv']  # n_tot entries
        resp = np.zeros_like(pwv)
        resp[tod_sel] = np.sum(data['resp_sel']*baseline[None,:],
                               axis=1)
        # similarly for some other fields
        sel = np.zeros_like(pwv)
        sel[tod_sel] = np.sum(data['sel'], axis=1)
        presel = np.zeros_like(pwv)
        presel[tod_sel] = np.sum(data['presel'], axis=1)
        # we want to plot with pwv sorted in order
        ind = np.argsort(pwv)
        # start plotting
        set_plotstyle()  # default cutslib style
        fig, ax = plt.subplots(1,1,figsize=(10,8))
        # ax.plot(pwv[ind], resp[ind], label='with resp')
        # ax.plot(pwv[ind], sel[ind], label='uncut')
        # plt.plot(pwv[ind], presel[ind], label='presel')
        plot_smoothed(pwv[ind], resp[ind], ax=ax, label='with resp', window=window)
        plot_smoothed(pwv[ind], sel[ind], ax=ax, label='uncut', window=window)
        plot_smoothed(pwv[ind], presel[ind], ax=ax, label='presel', window=window)
        crit_keys = [k for k in list(data.keys()) if 'Live' in k]
        for k in crit_keys:
            if k in ['skewLive', 'kurtLive']:
                continue
            crit_num = np.zeros_like(pwv)
            crit_num[tod_sel] = np.sum(data[k]*data['resp_sel'],axis=1)
            plot_smoothed(pwv[ind], crit_num[ind], ax=ax, label=k, window=window)
        # plot reference line
        if show_baseline:
            ax.axhline(np.sum(baseline), label='TES', color='k', linestyle='-')
        ax.axhline(np.sum(data['ff_sel']*baseline), label='with ff',
                   color='k', linestyle='--')
        ax.axhline(np.sum(data['ff_stable']*baseline), color='k',
                   label='with ff stable', linestyle=':')
        ax.legend(loc='best',frameon=False,bbox_to_anchor=(1, -0.1),
                  ncol=3)
        ax.set_xlabel("pwv/sin(alt)/mm")
        ax.set_ylabel("num of dets")
        ax.set_xlim([0, max(pwv)])  # not show -1
        plt.title("num of dets of various cases (smoothed)")
        outfile = op.join(p.o.patho.root, "cuts_hist.png")
        plt.savefig(outfile, bbox_inches='tight')
        print("Written to:", outfile)
        plt.close()
