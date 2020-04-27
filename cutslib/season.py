"""Interactive with season stats."""

# general dependency
import numpy as np, pickle, copy, os.path as op
from matplotlib import pyplot as plt
from functools import reduce

# cutslib dependency
from .depot import Depot
from .util import update_if_not_exist, get_rundir, tag_to_afsv
from .pathologyReport import pathoReport


class SeasonStats:
    def __init__(self, tag, calibrate=False, use_theta2=False):
        """Show the season stats using the collected pickle
        file containing all pathological parameters

        Parameters
        ----------
        tag: tag associated with cuts run

        """
        # store metadata
        array, freq, season, ver = tag_to_afsv(tag)
        self.array, self.freq = array, freq
        self.season, self.ver = season, ver
        # store a selection for reference
        # define default plotting style, common style
        common_style = {
            'extend': 2,
            'type': 'percentile',
            'scale': 'log'
        }
        # field specific styles
        self.style = {
            'gainLive': {
                'name': 'Gain',
                'limits': [5, 95],
            },
            'corrLive': {
                'name': 'Corr',
                'scale': 'linear',
                'limits': [0.95, None],
                'type': 'abs',
            },
            'normLive': {
                'name': 'Norm',
                'limits': [5, 95],
            },
            'DELive': {
                'name': 'Drift Error',
                'limits': [5, 95],
            },
            'MFELive': {
                'name': "Mid-Frequency Error",
                'limits': [1, 99],
            },
            'rmsLive': {
                'name': "RMS",
                'limits': [5, 95],
            },
            'skewLive': {
                'name': "Skew",
                'scale': 'linear',
                'limits': [1, 99],
            },
            'kurtLive': {
                'name': 'Kurtosis',
                'scale': 'linear',
                'limits': [1, 99],
            },
            'jumpLive': {
                'name': 'Jump',
                'limits': [1, 99],
            }
        }
        # merge common styles into field specific ones
        # unless the field is specified
        for k, v in self.style.items():
            update_if_not_exist(v, common_style)
        # load pickle file
        pickle_file = Depot().get_deep(('Postprocess', tag,
                                        f'{tag}_results.pickle'))
        with open(pickle_file, "rb") as f:
            data = pickle.load(f)
        if calibrate:
            data['rmsLive'] *= data['resp'] * data['ff'][:,None]
            data['normLive'] *= data['resp'] * data['ff'][:,None]
            data['MFELive'] *= data['resp'] * data['ff'][:,None]
            data['DELive'] *= data['resp'] * data['ff'][:,None]
            data['jumpLive'] *= data['resp'] * data['ff'][:,None]
        if use_theta2:
            # treat corr as cos\theta ~ 1-\theta^2/2
            # -> \theta^2 = 2(1-\cos\theta)
            # if \theta is gaussian, theta^2 should be chi-squared
            data['corrLive'] = 2*(1-data['corrLive'])
            # also update style accordingly
            self.style['corrLive'].update({
                'name': 'theta^2',
                'scale': 'log',
                'limits': [1, 99],
                'type': 'percentile',
            })
        self.stats = data
        print(f"stats loaded with {len(data['name'])} tods")
        # save an empty sel for future restriction work
        self.sel = np.ones_like(data['sel'], dtype=bool)
        # load db file from the run that contains some cuts stats
        run_dir = get_rundir(tag)
        db = op.join(run_dir, f"{tag}.db")
        if op.exists(db):
            self.db = pathoReport(db)
            print("patho report loaded")

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]
        if item in self.stats:
            return self.stats[item]
        elif f"{item}Live" in self.stats:
            return self.stats[f"{item}Live"]
        else: raise ValueError(f"{item} does not exist!")

    # plotting utilities

    def count_matches(self, *conds):
        """multiple conditions"""
        conds += (self.stats['sel'],)
        sel = reduce(np.logical_and, conds)
        fig, axes = plt.subplots(2,2,figsize=(12,10))
        axes[0,0].plot(self.ctime, np.sum(sel, axis=0),'k.',alpha=0.3)
        axes[0,0].set_xlabel('ctime')
        ax = axes[0,0].twinx()
        idx = np.argsort(self.ctime)
        ax.plot(self.ctime[idx], self.pwv[idx], 'k-',alpha=0.1)
        ax.set_ylabel('pwv')
        ax.set_ylim(bottom=0)
        axes[0,1].plot(np.sum(sel, axis=1),'k.',alpha=0.3)
        axes[0,1].set_xlabel('det_uid')
        det_uid = np.arange(sel.shape[0])
        row, col = det_uid // 32, det_uid % 32
        axes[1,0].plot(row, np.sum(sel, axis=1),'k.',alpha=0.3)
        axes[1,0].set_xlabel('row')
        axes[1,1].plot(col, np.sum(sel, axis=1),'k.',alpha=0.3)
        axes[1,1].set_xlabel('col')
        for i in range(2):
            for j in range(2):
                axes[i,j].set_ylabel('n matches')
                axes[i,j].set_ylim(bottom=1)
        plt.tight_layout()
        self.sel = sel
        return self

    def hist(self, figsize=(20, 12), nbins=100, style={}, hist_opts={}, guideline=True):
        data = self.stats
        sel = data['sel'].astype(bool)
        sel = np.logical_and(self.sel, sel)
        # start to plot
        # get styles right -> make sure we don't overwrite defaults
        mystyle = copy.deepcopy(self.style)
        mystyle.update(style)
        style = mystyle  # a better name
        fields = list(style.keys())
        fig, axes = plt.subplots(3, 3, figsize=figsize)
        for i in range(len(fields)):
            ax = axes[i//3, i%3]
            f = fields[i]
            d = data[f][sel]
            # find axis limit
            lo, hi = self._find_limits(d, style[f])
            # get bins right
            if style[f]['scale'] == 'log':
                bins = np.logspace(np.log10(lo), np.log10(hi), nbins)
            else: bins = np.linspace(lo, hi, nbins)
            # actually plot it
            ax.hist(d, bins=bins, **hist_opts)
            # get axis scale right
            if style[f]['scale'] == 'log':
                ax.set_xscale('log')
            ax.set_title(style[f]['name'])
            ax.get_yaxis().set_visible(False)


    def tri(self, figsize=(20, 20), nbins=100, style={}, hist_opts={}, hist2d_opts={},
            filename=None, density=False):
        data = self.stats
        sel = data['sel'].astype(bool)
        # start to plot
        # get styles right -> make sure we don't overwrite defaults
        mystyle = copy.deepcopy(self.style)
        mystyle.update(style)
        style = mystyle  # a better name
        # loop over pairs of crit
        fields = list(style.keys())
        fig, axes = plt.subplots(len(fields), len(fields), figsize=figsize)
        for i in range(len(fields)):
            for j in range(len(fields)):
                if j == i:  # plot hist
                    f = fields[j]
                    d = data[f][sel]
                    # find axis limit
                    lo, hi = self._find_limits(d, style[f])
                    # get bins right
                    if style[f]['scale'] == 'log':
                        bins = np.logspace(np.log10(lo), np.log10(hi), nbins)
                    else:
                        bins = np.linspace(lo, hi, nbins)
                    # actually plot it
                    opts = {'density': density}
                    opts.update(hist_opts)
                    axes[i,j].hist(d, bins=bins, **opts)
                    # get axis scale right
                    if style[f]['scale'] == 'log':
                        axes[i,j].set_xscale('log')
                    axes[i,j].set_title(style[f]['name'])
                    axes[i,j].get_yaxis().set_visible(False)
                elif j < i:  # scatter plot
                    f1, f2 = fields[j], fields[i]  # x, y
                    d1, d2 = data[f1][sel], data[f2][sel]
                    # find axes limits
                    lo1, hi1 = self._find_limits(d1, style[f1])
                    lo2, hi2 = self._find_limits(d2, style[f2])
                    if style[f1]['scale'] == 'log':
                        bins1 = np.logspace(np.log10(lo1), np.log10(hi1), nbins)
                    else:
                        bins1 = np.linspace(lo1, hi1, nbins)
                    if style[f2]['scale'] == 'log':
                        bins2 = np.logspace(np.log10(lo2), np.log10(hi2), nbins)
                    else:
                        bins2 = np.linspace(lo2, hi2, nbins)
                    # set up plot opts
                    opts = {'cmap': plt.cm.RdYlBu_r, 'density': density}
                    opts.update(hist2d_opts)
                    # actually plot it
                    axes[i,j].hist2d(d1, d2, bins=[bins1, bins2], **opts)
                    # get axis scale right
                    if style[f1]['scale'] == 'log':
                        axes[i,j].set_xscale('log')
                    if style[f2]['scale'] == 'log':
                        axes[i,j].set_yscale('log')
                    # get axis range right
                    axes[i,j].set_xlim([lo1, hi1])
                    axes[i,j].set_ylim([lo2, hi2])
                    # get ticks right, which is hard
                    if j != 0:  # only show yticks on first col
                        axes[i,j].get_yaxis().set_visible(False)
                    else:  # show label on first col
                        axes[i,j].set_ylabel(style[f2]['name'])
                    if i != len(fields)-1:  # only show xticks on last row
                        axes[i,j].get_xaxis().set_visible(False)
                    else:  # show label on last row
                        axes[i,j].set_xlabel(style[f1]['name'])
                else:
                    axes[i,j].axis('off')
        plt.tight_layout()
        if filename:
            plt.savefig(filename)

    def reset(self):
        """reset sel"""
        self.sel = np.ones_like(self.stats['sel'], dtype=bool)
        return self

    def _find_limits(self, d, style):
        """find limits of plotting, d is data vector,
        style is a dict with fields like scale, limits, etc"""

        # find axis limit
        lo, hi = style['limits']
        if style['type'] == 'percentile':
            if not lo: lo = 0
            if not hi: hi = 100
            lo, hi = np.percentile(d, [lo, hi])
            # extend the range if extend is specified
            extend = style['extend']
            if style['scale'] == 'linear':
                ctr = (lo+hi)/2
                lo, hi = ctr-(ctr-lo)*extend, ctr+(hi-ctr)*extend
            else:
                lo, hi = lo/extend, hi*extend
        elif style['type'] == 'abs':
            if not lo: lo = np.min(d)
            if not hi: hi = np.max(d)
        return lo, hi
