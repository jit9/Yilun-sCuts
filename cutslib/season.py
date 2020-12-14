"""Interactive with season stats."""

# general dependency
import numpy as np, pickle, copy, os.path as op, pandas as pd
import h5py, moby2
from matplotlib import pyplot as plt
from functools import reduce
from tqdm import tqdm
from scipy.signal import savgol_filter

# cutslib dependency
from .depot import Depot, SharedDepot
from .util import update_if_not_exist, get_rundir, tag_to_afsv, deep_merge
from .util import get_cutParam
from .pathologyReport import pathoReport, pathoList
from .visual import array_plots
from .catalog import Catalog


class SeasonStats:
    def __init__(self, tag=None, depot=None, calibrate=False, abscal='201026',
                 use_theta2=False, sort=False, verbose=False, planet=True, rundb=True):
        """Show the season stats using the collected pickle file containing
        all pathological parameters. The most important attribute is called
        style which contains everything a plotting function needs to know such
        as a reasonable plotting range. It will also contain the cut criteria.
        The style attr has a structure like the following

        self.style = {'field_name': {
            'name': 'field name',
            'limits': [5, 95],    # plotting limits
            'type': 'percentile', # unit of limits, can be abs or percentile
            'extend': 3,          # if we want to extend the limits
            'scale': 'log'        # x-scale for plotting
            'crit': [lo, hi]      # cut thresholds
          },
          ...
        }
        Changing the criteria can be done by changing the crit field. The actual pickle
        parameters are stored inside self.ss.

        Parameters
        ----------
        tag: tag associated with cuts run
        depot: depot to load pickle file
        calibrate: whether to calibrate the pickle parameters
        abscal: if not None, it will attempt to load the abscal with the given tag
        use_theta2: whether one wants to use theta^2 in place of corrLive which is
          no thing but cos\theta ~ 1-\theta^2/2
        sort: whether to sort the loaded pickle file by ctime
        verbose: verbosity flag
        planet: if True, it will attempt to load planet calibration data if it exists
        rundb: if True, it will attempt to load the runtime output file that contains
          some useful statistics of the cuts run
        """
        # store metadata
        self.tag = tag
        array, freq, season, ver = tag_to_afsv(tag)
        self.array, self.freq = array, freq
        self.season, self.ver = season, ver
        # store a selection for reference
        # define default plotting style, common style
        common_style = {
            'extend': 3,
            'type': 'percentile',
            'scale': 'log'
        }
        # field specific styles
        # unfortunately i have now start to use this dict for all sort of things
        # and style is no longer an adequete name. I will leave it like this
        # for the time-being, just to note that this is not purely style but
        # also defining crits that lots of functions depend on
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
        pickle_file = Depot(depot).get_deep(('Postprocess', tag,
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
        if verbose: print(f"stats loaded with {len(data['name'])} tods in ss.stats")
        # save an empty sel for future restriction work
        self.select = np.ones_like(data['sel'], dtype=bool)
        # load db file from the run that contains some cuts stats
        if rundb:
            run_dir = get_rundir(tag)
            db = op.join(run_dir, f"{tag}.db")
            if op.exists(db):
                self.db = pathoReport(db)
                self.db.addPWV()
                if verbose: print("patho report loaded in ss.db")
                # get catalog merged in with patho db
                cat = Catalog().narrow_down(list(self.stats['name']))
                self.db.data = self.db.data.merge(cat.data,left_on='todName',right_on='tod_name')
                if verbose: print("catalog merged in ss.db")
                # get patholist
                self.pl = pathoList(db)
                if verbose: print("pathoList loaded in ss.pl")
        # load cutparam
        self.cutParam = moby2.util.MobyDict.from_file(get_cutParam(tag))
        # populate cuts crit applied
        live_par = self.cutParam.get_deep(('pathologyParams', 'liveSelParams'))
        for k, v in live_par.items():
            f = f"{k}Live"
            if f in self.style and v['apply']:
                self.style[f]['crit'] = v['absCrit']
        # make sure theta2 crit is propogated
        if use_theta2 and 'crit' in self.style['corrLive']:
            crit = self.style['corrLive']['crit']
            crit[0] = 2*(1-crit[0])
            crit[1] = 2*(1-crit[1])
            self.style['corrLive']['crit'] = crit
        if verbose: print("cuts thresholds loaded in ss.styles[*]['crit']")
        # try to load planet calibration dataframe
        if planet:
            planet_file = Depot(depot).get_deep(('Postprocess', tag, 'calibration',
                                                 f'{tag}.csv'))
            if op.exists(planet_file):
                self.planet = pd.read_csv(planet_file)
                if verbose: print("planet measurements loaded ss.planet:", planet_file)
        # load absolute calibration if that's what we want
        if abscal:
            abscal_file = SharedDepot().get_deep(('TODAbsCal',f'abscal_{abscal}.h5'))
            if verbose: print("Loading abscal:", abscal_file)
            with h5py.File(abscal_file, "r") as f:
                abscal_data = f['abscal'][:]
                # choose only one freq
                bmask = abscal_data['band_id'].astype(str) == f"f{self.freq:03d}"
                tod_id = abscal_data['tod_id'].astype(str)[bmask]
                cal = abscal_data['cal'][bmask]
                del abscal_data, bmask
            # next we want to get the abscal for the tods we have
            # this method below is okay but too slow
            # idx = np.array([np.where(tod_id == n)[0] for n in self.name])
            # here is a slightly faster approach
            inc_idx = np.where(np.isin(tod_id, self.name))[0]
            mat_idx = [np.where(self.name==tod_id[i])[0][0] for i in inc_idx]
            idx = inc_idx[np.argsort(mat_idx)]
            self.abscal = cal[idx]

        # sort values if that's what we want
        if sort: self.sort_values()

    def get_subset(self, todlist, verbose=True):
        match = np.isin(self.name, todlist)
        if verbose: print(f"match tods: {np.sum(match)}")
        new_ss = {}
        for k, v in self.ss.items():
            if v.shape[-1] == len(self.name):
                new_ss[k] = v[...,match]
        return new_ss

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]
        if item in self.stats:
            return self.stats[item]
        elif f"{item}Live" in self.stats:
            return self.stats[f"{item}Live"]
        else: raise ValueError(f"{item} does not exist!")

    def print_thresholds(self):
        for f in self.style:
            if 'crit' in self.style[f]:
                print(f"{f}: {self.style[f]['crit']}")

    def sel2hdf(self, filename, sel=None):
        """Save a given sel or the internal sel into a hdf file
        to be later digested by cuts pipeline"""
        if sel is None: sel = self.select
        assert sel.shape == self.stats['sel'].shape
        f = h5py.File(filename, "w")
        for i in tqdm(range(sel.shape[-1])):
            det_cuts = sel[:,i]
            # skip if nothing is cut
            if np.sum(det_cuts) == 0: continue
            obs = self.stats['name'][i]
            f[obs] = det_cuts

    def hdf2sel(self, filename):
        """Load hdf file into a sel"""
        f = h5py.File(filename, "r")
        sel = np.zeros_like(self.select)
        for i in range(sel.shape[-1]):
            obs = self.stats['name'][i]
            if obs in f:
                sel[:,i] = f[obs][:]
        return sel

    def sort_values(self):
        """Sort stats and db in the same order by ctime, since i don't have
        multiple observations at the same time now that we are using single
        frequency and array stats.

        """
        # depending on the meaning of the fields (axes), sort them accordingly. Here i am
        # avoiding have a general loop that checks the shape of things before sorting, because
        # i don't want to have the bug of sometimes having 1024 TODs get mistaken as 1024 dets, etc.
        t_fields = ['name', 'scan_freq', 'ctime', 'alt', 'pwv', 'tod_sel', 'tes_sel']
        td_fields = ['sel', 'psel', 'resp', 'resp_sel', 'cal','gainLive', 'gainLive_sel', \
                     'corrLive', 'corrLive_sel', 'normLive', 'normLive_sel', 'rmsLive', 'rmsLive_sel', \
                     'kurtLive', 'kurtLive_sel', 'skewLive', 'skewLive_sel', 'MFELive', 'MFELive_sel', \
                     'DELive', 'DELive_sel', 'jumpLive', 'jumpLive_sel']
        # sorted index
        sorted_idx = np.argsort(self.stats['ctime'])
        for f in t_fields:
            self.stats[f] = self.stats[f][sorted_idx]
        for f in td_fields:
            self.stats[f] = self.stats[f][:,sorted_idx]
        print("ss.stats sorted by ctime")
        # sort patho db in self.db if that's loaded
        if hasattr(self,'db'):
            self.db.data = self.db.data.sort_values(by=['ctime_x'])
            print("ss.db sorted by ctime")

    def update_critsel(self, verbose=True):
        """redo the sel based on new crit defined in the class attribute style, this
        goes through the list of pathologies and get sel based on the abs crit specified
        in the style attribute of this object, recalculate the total sel and repopulate
        that into the dictionary"""
        fields = [f for f in self.style if 'crit' in self.style[f]]
        for f in fields:
            # get crit corresponding to the field
            crit = self.style[f]['crit']
            # get pathology values
            values = self.stats[f]
            # apply the cuts
            sel = np.ones_like(values, dtype=bool)
            lo, hi = crit
            if crit[0] > crit[1]: lo, hi = crit[1], crit[0]
            if lo: sel *= values > lo
            if hi: sel *= values < hi
            # repopulate the sel back to the stats
            self.stats[f"{f}_sel"] = sel
            if verbose: print(f"-> {f}_sel updated: {np.sum(sel)} dets passed")
        # update total sel
        self.stats['sel'] = reduce(np.logical_and, [self.stats[f"{f}_sel"] for f in fields])
        if verbose: print("-> sel updated")
        return self

    def update_style(self, style={}):
        """Update the internal style in place, one use case is to change the crit and
        then call update_critsel to regenerate det cuts"""
        self.style = deep_merge(self.style, style)
        return self

    # plotting utilities
    def find_matches(self, *conds, alone=False):
        """Plot the detector matching multiple conditions against
        various quantities such as alt, loading, row, col for debugging purpose.
        It can be called as the following
        ss.find_matches(cond1, cond2, cond3, ...)

        Parameters
        ----------
        conds: list of different conditions
        alone: if True, it will bypass the original sel
        """
        if not alone: conds += (self.stats['sel'],)
        sel = reduce(np.logical_and, conds)
        fig, axes = plt.subplots(2,3,figsize=(15,10))
        axes[0,0].plot(self.ctime, np.sum(sel, axis=0),'k.',alpha=0.3)
        axes[0,0].set_xlabel('ctime')
        ax = axes[0,0].twinx()
        idx = np.argsort(self.ctime)
        ax.plot(self.ctime[idx], self.pwv[idx], 'k-',alpha=0.1)
        ax.set_ylabel('pwv')
        ax.set_ylim(bottom=0)
        axes[0,2].plot(self.alt, np.sum(sel, axis=0),'k.',alpha=0.3)
        axes[0,2].set_xlabel('alt')
        axes[0,1].plot(self.pwv, np.sum(sel, axis=0), 'k.', alpha=0.3)
        axes[0,1].set_xlim([0,3])
        axes[0,1].set_xlabel('loading (mm)')
        det_uid = np.arange(sel.shape[0])
        row, col = det_uid // 32, det_uid % 32
        axes[1,0].plot(np.sum(sel, axis=1),'k.',alpha=0.3)
        axes[1,0].set_xlabel('det_uid')
        axes[1,1].plot(row, np.sum(sel, axis=1),'k.',alpha=0.3)
        axes[1,1].set_xlabel('row')
        axes[1,2].plot(col, np.sum(sel, axis=1),'k.',alpha=0.3)
        axes[1,2].set_xlabel('col')
        for i in range(2):
            for j in range(3):
                axes[i,j].set_ylabel('n matches')
                axes[i,j].set_ylim(bottom=1)
        plt.tight_layout()
        self.select = sel
        return self

    def view_sel(self, sel=None):
        """View the selected mask as a function of various things for debugging"""
        if sel is None: sel = self.stats['sel']
        fig, axes = plt.subplots(3,1,figsize=(20, 10), sharex=True)
        idx = np.argsort(self.ctime)
        axes[0].plot(self.pwv[idx], 'k.', markersize=1)
        axes[0].set_ylim([0,3.5])
        axes[0].xaxis.set_visible(False)
        axes[0].set_ylabel('PWV / sin(alt) (cm)')
        axes[1].imshow(sel[:,idx], aspect="auto", cmap='Greys', origin='lower')
        axes[1].set_ylabel('Dets')
        axes[2].plot(np.sum(sel[:,idx], axis=0), 'k-', lw=0.2)
        axes[2].set_xlabel('TOD (ordered by ctime)')
        axes[2].set_ylabel('# of Live Dets')
        fig.subplots_adjust(hspace=0)

    def view_cuts(self, window=None, gain=True, mfe=True, jump=True, kurt=False, skew=False, **kwargs):
        """Plot how individual crit cut changes with pwv.

        """
        exclude = ['psel']
        if not kurt: exclude += ['kurtLive_sel']
        if not skew: exclude += ['skewLive_sel']
        if not jump: exclude += ['jumpLive_sel']
        if not gain: exclude += ['gainLive_sel']
        if not  mfe: exclude += ['MFELive_sel']
        fields = [k for k in self.stats if 'sel' in k and k not in exclude]
        if not window: window = int(self.select.shape[-1]/10)
        if window % 2 ==0: window += 1
        plt.figure(figsize=(10,8))
        for f in fields:
            if len(self.stats[f].shape) == 1: continue  # only 2d sel
            idx = np.argsort(self.pwv)
            # corrected by resp sel
            remain = np.sum(self.stats[f]*self.tes_sel[:,None], axis=0)[idx]
            # smooth remain
            remain = savgol_filter(remain, window_length=window, polyorder=2)
            plt.plot(self.pwv[idx], remain, '-', markersize=1, label=f, **kwargs)
            plt.xlim([0,4])
            plt.xticks(np.linspace(0, 4, 17))
            plt.xlabel('Loading (mm)')
            plt.ylabel('# of remaining dets')
        plt.legend()

    def hist_pwv(self, figsize=(8,6), **kwargs):
        """View histogram of PWV in the season"""
        fig, ax = plt.subplots(1,1,figsize=figsize)
        pwv = self.pwv[self.pwv>0]
        opts = {'rwidth': 0.75, 'bins': np.linspace(0,4)}
        opts.update(kwargs)
        ax.hist(pwv, **opts);
        ax.set_xlim([0,3])
        ax.set_xlabel('PWV/sin(alt) (mm)')
        ax.set_ylabel('# of TODs')

    def planet_peaks(self, c=None, ylim=[0,1e-10]):
        """View the planet peak measurements as function of optical loading"""
        assert hasattr(self, 'planet'), "Planet measurements not available!"
        df = self.planet
        plt.figure(figsize=(8,6))
        if c is not None:
            if isinstance(c, str): c=df[c]
            plt.scatter(df.loading, df.peak_mean, s=100, c=c,
                        edgecolor='None', cmap='copper')
        else: plt.plot(df.loading, df.peak_mean, '.')
        plt.xlabel('Loading (PWV / sin(alt))')
        plt.ylabel('Peak measurements')
        plt.ylim(ylim)
        clb = plt.colorbar()
        if isinstance(c, str): clb.set_label(c)

    def array_plots(self, field, dets=None, **kwargs):
        """Quick array plots for a given field"""
        if isinstance(field, str): field = getattr(self, field)
        if dets is None: dets = np.arange(self.stats['sel'].shape[0])
        assert dets.shape == field.shape, "dets and field mismatch!"
        array_plots(field, det=dets, array=self.array,season=self.season,fr=self.freq, **kwargs)

    def view_hist(self, field, sel=None, nbins=100, hist_opts={}, **kwargs):
        """View the histogram of a particular criteria field"""
        if sel is None: sel = self.select
        # get style
        style = copy.deepcopy(self.style[field])
        style.update(kwargs)
        fig, ax = plt.subplots(1,1,figsize=(10,8))
        d = self.stats[field][sel]
        lo, hi = self._find_limits(d, style)
        if lo == 0: lo += 0.01
        # get bins right
        if style['scale'] == 'log':
            bins = np.logspace(np.log10(lo), np.log10(hi), nbins)
        else: bins = np.linspace(lo, hi, nbins)
        # actually plot it
        ax.hist(d, bins=bins, **hist_opts)
        # show threshold if we have thresholds loaded
        if 'crit' in style:
            crit = style['crit']
            ax.axvline(crit[0],color='r',ls='--')
            ax.axvline(crit[1],color='r',ls='--')
        # get axis scale right
        if style['scale'] == 'log':
            ax.set_xscale('log')
        ax.set_title(style['name'])
        ax.get_yaxis().set_visible(False)

    def hist(self, sel=None, figsize=(20, 12), nbins=100, style={}, hist_opts={}, show_crit=True,
             show_sel=False, axes=None, show_all=False):
        """View histograms of all crit fields"""
        data = self.stats
        if sel is None:
            sel = data['sel'].astype(bool)
        # start to plot
        # get styles right -> make sure we don't overwrite defaults
        mystyle = copy.deepcopy(self.style)
        mystyle = deep_merge(mystyle, style)
        style = mystyle  # a better name
        fields = list(style.keys())
        if axes is None:
            fig, axes = plt.subplots(3, 3, figsize=figsize)
        for i in range(len(fields)):
            ax = axes[i//3, i%3]
            f = fields[i]
            if not show_all: d = data[f][sel]
            else: d = data[f]
            # find axis limit
            lo, hi = self._find_limits(d, style[f])
            # get bins right
            if style[f]['scale'] == 'log':
                bins = np.logspace(np.log10(lo), np.log10(hi), nbins)
            else: bins = np.linspace(lo, hi, nbins)
            # actually plot it
            ax.hist(d, bins=bins, **hist_opts)
            if show_crit and 'crit' in style[f]:
                crit = style[f]['crit']
                ax.axvline(crit[0],color='r',ls='--')
                ax.axvline(crit[1],color='r',ls='--')
            # get axis scale right
            if style[f]['scale'] == 'log':
                ax.set_xscale('log')
            ax.set_title(style[f]['name'])
            ax.get_yaxis().set_visible(False)
        if show_sel:
            axes = self.hist(self.select, axes=axes, hist_opts={'color':'r'})
        return axes

    def tri(self, figsize=(20, 20), nbins=100, style={}, hist_opts={}, hist2d_opts={},
            filename=None, density=False):
        data = self.stats
        sel = data['sel'].astype(bool)
        # start to plot
        # get styles right -> make sure we don't overwrite defaults
        mystyle = copy.deepcopy(self.style)
        mystyle = deep_merge(mystyle, style)
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
        self.select = np.ones_like(self.stats['sel'], dtype=bool)
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

    def report_tods(self, tods):
        """report the cuts for tods"""
        match = np.isin(self.db.data.todName, tods)
        columns = ['todName', 'length', 'liveDets', 'corrLive', 'normLive',
                   'DELive', 'MFELive', 'rmsLive', 'skewLive', 'kurtLive',
                   'corrLive_m', 'corrLive_s', 'normLive_m', 'normLive_s', 'DELive_m',
                   'DELive_s', 'MFELive_m', 'MFELive_s', 'rmsLive_m', 'rmsLive_s',
                   'skewLive_m', 'skewLive_s', 'kurtLive_m', 'kurtLive_s',
                   'obs_detail', 'date', 'hour_utc', 'alt', 'loading', 'pwv_source']
        return self.db.data[columns][match]

    def plot_stats(self, field, dets=None, nrand=10, crange=None,
                   hour=True, ylim=None, highlights=None, abscal=True, ylabel='',
                   title='', op=np.abs, dot_alpha=1):
        """Plot stats as a function of time
        Parameters
        ----------
        dets: detector list of interests
        nrand: randomly sample a few detectors to look at
        crange: ctime range of interests i.e. [1555000000, 1556000000]
        hour: whether to use hour as xaxis for the plot
        highlights: list of tods to highlight specifically

        """
        if dets is None:
            dets = np.where(self.ff_sel * self.tes_sel)[0]
            if nrand: dets = np.random.choice(dets, nrand, replace=False)
        # plot calibration with uncut dets
        if isinstance(field, str):
            field = np.ma.array(op(getattr(self, field)))
            field[~self.sel] = np.ma.masked
        else:
            field = field.copy()
        # apply absolute calibration per tod to account for pwv effects
        if abscal and hasattr(self, 'abscal'): field *= self.abscal[None,:]
        # xaxis
        if hour:
            tfunc = lambda x: (x - self.ctime.min())/3600
            ctime = tfunc(self.ctime)
            xlabel = f"hours since {self.ctime.min()} [h]"
        else:
            ctime = self.ctime
            xlabel = "ctime [s]"
        # unless required to show highlighted only, always plot all
        lines = plt.plot(ctime, field[dets].T, '.', markersize=1, alpha=dot_alpha)
        # highlight a list of tods if necessarily
        if highlights is not None:
            match = np.isin(self.name, highlights)
            print(f"match tod: {np.sum(match)}")
            plt.plot(ctime[match], field[np.ix_(dets, match)].T, 'x')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        # legend
        if len(dets)<=20:
            plt.legend(iter(lines), dets, bbox_to_anchor=(1.1,1),
                       loc="upper left", ncol=int(np.ceil(len(dets)/10)))
        if ylim is not None:
            plt.ylim(ylim)
        # also plot pwv
        pwv_ax = plt.gca().twinx()
        idx = np.argsort(ctime)
        pwv_ax.plot(ctime[idx], self.pwv[idx], 'k-', alpha=0.2, label='pwv/sin(alt)')
        pwv_ax.set_ylabel('pwv / sin(alt) [mm]')
        pwv_ax.set_ylim([0, 6])
        plt.legend()
        # ctime range
        if crange is not None:
            cstart, cend = crange
            cstart = max(self.ctime.min(), int(cstart))
            cend   = min(self.ctime.max(), int(cend))
            if hour: cstart, cend = tfunc(cstart), tfunc(cend)
            plt.xlim([cstart, cend])
        plt.title(title)

    def plot_cal(self, dets=None, nrand=10, crange=None, hour=True,
                 ylim=None, highlights=None, abscal=True, dot_alpha=1):
        return self.plot_stats('cal', dets=dets, nrand=nrand, dot_alpha=dot_alpha,
                               crange=crange, hour=hour, ylim=ylim, highlights=highlights,
                               abscal=abscal, title='calibration=ff*biasstep*abscal',
                               ylabel='calibration [uK/DAC]')
