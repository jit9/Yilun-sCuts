"""Utility functions related to glitch analysis"""
import numpy as np
from scipy import stats

import moby2

class TODSnippet:
    def __init__(self, tod=None, det_uid=None, tslice=None):
        """Create a TOD snippet with a subset of dets and a slice of samples

        Parameters
        ----------
        tod: base TOD object
        det_uid: list of dets of interests
        tslice: a slice object to get a subset of samples

        """
        if tod is None: return
        self.det_uid = det_uid
        self.data = tod.data[det_uid,tslice]
        self.tslice = tslice
        self.info = SnippetInfo.from_todinfo(tod.info)
    def demean(self):
        """Remove the mean of the snippet"""
        self._mean = self.data.mean(axis=1)
        self.data -= self._mean[:,None]
        return self
    def deslope(self):
        self._slope = (self.data[:,-1] - self.data[:,0]) / self.data.shape[-1]
        self.data -= self._slope[:,None] * np.arange(self.data.shape[-1])
        return self
    def __repr__(self):
        return f"TODSnippet(ndet={len(self.det_uid)},tslice={self.tslice})"
    def plot(self, demean=False, deslope=False, debuffer=0, **kwargs):
        return plot_snippet(self, demean, deslope, debuffer, **kwargs)
    def peaks_radius(self):
        """Peaks radius from the maximum peak"""
        return peaks_radius_from_max(self)
    def peaks(self):
        """Get peak amplitudes"""
        return np.max(np.abs(self.data), axis=1)
    def plot_array(self, values=None):
        """Plot det_uid on the array"""
        from cutslib.visual import array_plots
        if values is None: values = np.ones_like(self.det_uid)
        array_plots(values, det=self.det_uid, season=self.info.season, array=self.info.array)
    def plot_row_col(self, values=None):
        from matplotlib import pyplot as plt
        dets = self.det_uid
        ad = self.info.array_data
        row, col = ad['row'][dets], ad['col'][dets]
        if values is None: values = np.ones_like(self.det_uid)
        plt.scatter(row, col, c=values, s=20)
        plt.xlabel('row')
        plt.ylabel('col')
        return plt.gca()
    def plot_peaks_radius(self, bins=10, method='linear'):
        r, p = self.peaks_radius()
        return plot_peaks_radius_binned(r, p, bins, method)
    def print_info(self, fields=None):
        dets = self.det_uid
        if fields is None: fields = list(self.info.array_data.keys())
        print('array_info:')
        for k in fields:
            print(f'  {k}: {self.info.array_data[k][dets]}')
    def snr_template(self, template):
        """Return the snr for each det in the snippet after correlating
        with a known template

        Parameters
        ----------
        template: 1d np.ndarray

        """
        return np.apply_along_axis(lambda x: template_match(x, template), 1, self.data)
    def max_snr_template(self, template):
        """Return the maximum snr in the snippet after correlating
        with a known template"""
        return max(self.snr_template(template))

class SnippetInfo:
    def __init__(self):
        """Subset of TODInfo to keep in each TODSnippet"""
        pass
    @classmethod
    def from_todinfo(cls, todinfo):
        # fields to copy over
        fields = ['filename', 'basename', 'name', 'tod_id', 'instrument',
                  'ctime', 'day_str', 'time_str', 'season', 'array', 'array_data',
                  'sample_index', 'det_uid']
        self = cls()
        for f in fields:
            self.__dict__[f] = getattr(todinfo, f)
        return self

def get_glitch_snippets(tod, dets, cv):
    """Generate TODSnippet from a given CutsVector

    Parameters
    ----------
    tod: TOD object
    dets: list of dets of interests
    rng: time index range of the snippet of form [i_l, i_h]
    rm: whether to remove mean in each snippets

    Returns
    -------
    [TODSnippet] (list of TODSnippet)

    """
    snippets = []
    tslices = cv2slices(cv)
    for s in tslices:
        snippets.append(TODSnippet(tod, dets, s).demean())
    return snippets

def affected_snippets_from_cv(tod, cuts, cv, dets):
    """Get snippets from a given cut vector while only maintaining
    those dets that are affected in each range in the cv.

    Parameters
    ----------
    tod: base TOD object
    cuts: TODCuts object containing the base cuts to extract affected dets from
    cv: a CutsVector object that specifies ranges of interests
    dets: a narrow-down list of dets to look at

    """
    dets_events = dets_affected_in_cv(cuts, cv, dets)
    snippets = []
    for d, s in zip(dets_events, cv2slices(cv)):
        snippets.append(TODSnippet(tod, d, s))
    return snippets


def peaks_radius_from_max(snippet, ref='array'):
    """Get the peaks and corresponding radius from the highest peak

    Parameters
    ----------
    snippet: TODSnippet object
    ref: specify what coord to use, can be either array or sky

    Returns
    -------
    radius, peaks

    """
    # first make sure snippet has a mean of zero
    data = snippet.data - snippet.data.mean(axis=1)[:,None]
    # find peak height
    peaks = np.max(np.abs(data), axis=1)
    # find highest peak as the reference point
    imax = np.argmax(peaks)
    x = snippet.info.array_data[f'{ref}_x'][snippet.det_uid]
    y = snippet.info.array_data[f'{ref}_y'][snippet.det_uid]
    x_c = x[imax]
    y_c = y[imax]
    r = np.sqrt((x-x_c)**2+(y-y_c)**2)
    return r, peaks

def glitch_det_count(cuts, dets=None):
    """Count number of dets affected as a time series

    Parameters
    ----------
    cuts: TODCuts object

    """
    if dets is None:
        return np.sum([c.get_mask() for c in cuts.cuts], axis=0)
    else:
        return np.sum([c.get_mask() for i, c in enumerate(cuts.cuts)
                       if cuts.det_uid[i] in dets], axis=0)

def pcuts2mask(cuts):
    """Convert partial cuts to a 2d boolean mask"""
    mask = np.stack([c.get_mask() for c in cuts.cuts], axis=0)
    return mask

def is_cut(cv, t):
    """check if a specific time is cut in a det (provided CutVector)

    Parameters
    ----------
    cv: CutsVector
    t: time index

    Returns
    -------
    True if t is cut else False

    """
    for c in cv:
        if c[0] <= t <= c[1]:
            return True
    return False

def dets_affected_at_t(cuts, t):
    """Find dets affected by the given cuts at a specific time t

    Parameters
    ----------
    cuts: TODCuts object
    t: time index

    Returns
    -------
    [det_uid]

    """
    return [d for (d, cv) in zip(cuts.det_uid, cuts.cuts) if is_cut(cv, t)]

def pixels_affected_at_t(cuts, t, pr):
    """Find pixels affected by the given cuts at a specific time t

    Parameters
    ----------
    cuts: TODCuts object
    t: time index
    pr: PixelReader object

    Returns
    -------
    [pixel_id]

    """
    dets = dets_affected_at_t(cuts, t)
    pixels = np.unique(pr.dets2pixels(dets))
    return pixels

def fill_cv(data, cv, fill_value=0, inplace=True):
    """Fill an array-like data with a CutsVector, by default
    it acts on the last axis.

    """
    ss = cv2slices(cv)
    if not inplace: data = data.copy()
    for s in ss:
        data[...,s] = fill_value
    return data

def dets_affected_in_cv(cuts, cv, dets):
    """Find detectors affected in each range in a CutsVector

    Parameters
    ----------
    cuts (TODCuts): base TODCuts object to gather the dets affected info
    cv (CutsVector): specify the ranges of samples to find dets affected
    dets (boolean array): an narrowed-down list of dets to look at

    """
    get_dets = lambda x: np.where((np.sum(x,axis=1)>0)*dets)[0]
    return slices_map(get_dets, pcuts2mask(cuts), cv2slices(cv))

class PixelReader:
    def __init__(self, season='2016', array='AR3', mask=None):
        """Utility class to find information about each pixel (feedhorn), migrated
        from todloop.util.pixels.PixelReader to manage internally"""
        self._array_info = {
            'season': season,
            'array_name': array
        }
        self._array_pos = None
        self._freqs = None
        self.ad = moby2.scripting.get_array_data(self._array_info)
        self._pixel_dict = self.generate_pixel_dict()
        self._mask = mask
        self.calibrate_array(season=self._array_info['season'])
        self.get_adjacent_detectors = self.adjacent_detector_generator()

    def generate_pixel_dict(self):
        """Generate pixel dictionary that tells which detectors correspond
        to which pixel and frequencies """
        self._array_pos = np.vstack([self.ad['array_x'], self.ad['array_y']]).T
        self._freqs = np.sort(np.unique(self.ad['nom_freq']))[1:]  # gather freqs (exclude 0)
        pixel_dict = {}  # initialize empty pixel_dict

        for det_id in self.ad['det_uid']:
            if np.all(self._array_pos[det_id, :] == [0, 0]):  # not physical
                continue
            if self.ad['det_type'][det_id] != 'tes':  # remove non-tes
                continue
            dets = np.where(np.all(self._array_pos == self._array_pos[det_id, :], axis=1))[0]
            # make a dictionary of frequencies: f1: lower freq, f2: higher freq
            pol_dict = {
                'f1': [i for i in dets if self.ad['nom_freq'][i] == self._freqs[0] and
                       self.ad['det_type'][det_id] == 'tes'],
                'f2': [i for i in dets if self.ad['nom_freq'][i] == self._freqs[1] and
                       self.ad['det_type'][det_id] == 'tes']
            }
            pixel_id = dets[0]  # index pixel by the smallest det_uid
            pixel_dict[str(pixel_id)] = pol_dict

        return pixel_dict

    def calibrate_array(self, season):
        """Calibrate the array_data based on season since different
        seasons have different array_data units"""
        if season == '2017':
            self.ad['array_x'] /= 10000.0
            self.ad['array_y'] /= 10000.0

    def adjacent_detector_generator(self):
        """Generate a get_adjacent_pixels function
        Return a function to get adjacent detectors

        return: [int] function(int det)
        """
        # Find the adjacent detectors
        # Generate an empty list to store adjacent detector lists
        adj_dets = [None] * len(self.ad['array_x'])
        ar = self._array_pos

        for i in range(len(self.ad)):
            dis = np.dot((ar - ar[i, :])**2, [[1], [1]])
            _mask = ((dis < 0.6) & (dis > 0)).T & ((ar[:, 0] != 0) | (ar[:, 1] != 0))
            mask = _mask.flatten()
            indexes = np.where(mask is True)
            adj_dets[i] = list(list(indexes)[0])  # Normalize the np array output to list

        # Generate a function to access the data to make sure above procedures run once only
        def get_adjacent_detectors(detector):
            return adj_dets[detector]

        return get_adjacent_detectors

    def get_pixels(self):
        return [int(key) for key in self._pixel_dict]

    def get_f1(self, pixel):
        if self._mask is not None:
            return [det for det in self._pixel_dict[str(pixel)]['f1'] if self._mask[det] == 1]
        else:
            return self._pixel_dict[str(pixel)]['f1']

    def get_f2(self, pixel):
        if self._mask is not None:
            return [det for det in self._pixel_dict[str(pixel)]['f2'] if self._mask[det] == 1]
        else:
            return self._pixel_dict[str(pixel)]['f2']

    def get_dets(self, pixel):
        if self._mask is not None:
            return [det for det in self._pixel_dict[str(pixel)]['f1'] if self._mask[det] == 1]
        else:
            return self._pixel_dict[str(pixel)]['f1']

    def get_adjacent_pixels(self, pixel):
        all_adj_det = self.get_adjacent_detectors(pixel)
        return [int(det) for det in all_adj_det if str(det) in self._pixel_dict]


    def get_pixels_within_radius(self, pixel, radius):
        ar = self._array_pos
        dist = np.sqrt(np.sum((ar - ar[pixel, :])**2, axis=1))
        return [det for det in np.arange(1056)[dist < radius] if str(det) in self._pixel_dict]

    def plot(self, pixels=None):
        plt.plot(self.ad['array_x'], self.ad['array_y'], 'r.')
        if pixels:
            plt.plot(self.ad['array_x'][pixels], self.ad['array_y'][pixels], 'b.')

    def get_x_y(self, pixel):
        return self.ad['array_x'][pixel], self.ad['array_y'][pixel]

    def get_x(self, pixel):
        return self.ad['array_x'][pixel]

    def get_y(self, pixel):
        return self.ad['array_y'][pixel]

    def get_row_col(self, pixel):
        """Return row and col of pixels
        :param:
            pixel:  int or [int]
        :return:
            row, col"""
        return self.ad['row'][pixel], self.ad['col'][pixel]

    def get_row(self, pixel):
        """Return row of pixel(s)
        :param: int or [int]
        :return: int or [int]"""
        return self.ad['row'][pixel]

    def get_col(self, pixel):
        """Return col of pixel(s)
        :param: int or [int]
        :return: int or [int]"""
        return self.ad['row'][pixel]

    def get_row_col_array(self):
        """Return row and col of all pixels (array)
        :param:
            pixel:  int or [int]
        :return:
            [row], [col]"""
        return self.ad['row'], self.ad['col']

    def get_x_y_array(self):
        """Get the xy of the entire array for plotting"""
        return self.ad['array_x'], self.ad['array_y']

    def dets2pixels(self, dets):
        """Convert a list of det_uid to pixel_id

        Parameters
        ----------
        dets: [det_uid]
        pr: PixelReader

        Returns
        -------
        [pixel_id] with the same shape as det_uid which means
          there will be duplicated entries

        """
        pixels = []
        for det in dets:
            pixel = [p for p in self.get_pixels() if det in self.get_dets(p)]
            assert len(pixel) == 1, f"Det {det} shows up in no / multiple pixels!"
            pixels.append(pixel[0])
        return pixels

    @classmethod
    def for_tod(cls, tod):
        return cls(season=tod.info.season, array=tod.info.array)

class CutsVector(moby2.tod.cuts.CutsVector):
    """Wrapper around moby2 version"""
    def __invert__(self):
        return self.get_complement()
    def __add__(self, other):
        if isinstance(other, CutsVector) or isinstance(other, moby2.tod.cuts.CutsVector):
            return CutsVector.from_mask(self.get_mask() + other.get_mask())
        else: return NotImplemented
    def __radd__(self, other):
        return self.__add__(other)
    def __mul__(self, other):
        if isinstance(other, CutsVector) or isinstance(other, moby2.tod.cuts.CutsVector):
            return CutsVector.from_mask(self.get_mask() * other.get_mask())
        else: return NotImplemented
    def __rmul__(self, other):
        return self.__mul__(other)

class CutsMatrix:
    """Similar to CutsVector but higher dimentional"""
    def __init__(self, cvs=None):
        """Wrapper for a list of CutsVector

        Parameters
        ----------
        cvs: list of CutsVector

        """
        self.cvs = cvs
    def __add__(self, other):
        if isinstance(other, CutsVector):
            return CutsMatrix([cv + other for cv in self.cvs])
        elif isinstance(other, CutsMatrix):
            assert self.shape == other.shape, "Shape mismatch!"
            return CutsMatrix([cv1 + cv2 for cv1, cv2 in zip(self.cvs, other.cvs)])
        else: return NotImplemented
    def __iadd__(self, other):
        if isinstance(other, CutsVector):
            self.cvs = [cv + other for cv in self.cvs]
        elif isinstance(other, CutsMatrix):
            assert self.shape == other.shape, "Shape mismatch!"
            self.cvs = [cv1 + cv2 for cv1, cv2 in zip(self.cvs, other.cvs)]
        else: return NotImplemented
    def __radd__(self, other):
        return self.__add__(self, other)
    def __mul__(self, other):
        if isinstance(other, CutsVector):
            return CutsMatrix([cv * other for cv in self.cvs])
        elif isinstance(other, CutsMatrix):
            assert self.shape == other.shape, "Shape mismatch!"
            return CutsMatrix([cv1 * cv2 for cv1, cv2 in zip(self.cvs, other.cvs)])
        else: return NotImplemented
    def __iadd__(self, other):
        if isinstance(other, CutsVector):
            self.cvs = [cv * other for cv in self.cvs]
        elif isinstance(other, CutsMatrix):
            assert self.shape == other.shape, "Shape mismatch!"
            self.cvs = [cv1 * cv2 for cv1, cv2 in zip(self.cvs, other.cvs)]
        else: return NotImplemented
    def __rmul__(self, other):
        return self.__mul__(self, other)
    def __invert__(self):
        return CutsMatrix([~cv for cv in self.cvs])
    def get_mask(self):
        return np.stack([c.get_mask() for c in self.cvs], axis=0)
    @classmethod
    def from_mask(cls, mask):
        return cls([CutsVector.from_mask(m) for m in mask])
    @property
    def shape(self):
        if self.cvs is None: return (0,)
        return (len(self.cvs), self.cvs[0].nsamps)
    def __repr__(self):
        return f"CutsMatrix(shape={self.shape})"

#############
# utilities #
#############

def bin_data(x, y, bins=10, method='linear', err_method='std'):
    """Calculate binned statistics

    Parameters
    ----------
    x, y: series to bin
    bins: number of bins
    method: how to bin

    Returns
    -------
    (binned x (centers), binned y, std in each bin)

    """
    from scipy import stats
    if method == 'linear': x_ = x
    elif method == 'log': x_ = np.log(x)
    elif method == 'p2': x_ = x**0.5
    else: raise ValueError("Unsupported binning method")
    res = stats.binned_statistic(x_, y, bins=bins, statistic='mean')
    res_std = stats.binned_statistic(x_, y, bins=bins, statistic=err_method)
    bc_ = (res.bin_edges[1:] + res.bin_edges[:-1])/2
    if method == 'linear': bc = bc_
    elif method == 'log': bc = np.exp(bc_)
    elif method == 'p2': bc = bc_**2
    else: raise ValueError("Unsupported binning method")
    return bc, res.statistic, res_std.statistic

def cv2slices(cv):
    """Convert CutsVector to a list ot time slices"""
    return [slice(v[0],v[1],None) for v in cv]

def slices_map(func, data, slices):
    """Apply a function on the data inside given slices, by default
    the slice is applied on the last axis

    """
    res = []
    for s in slices:
        res.append(func(data[...,s]))
    return res

def template_match(data, template):
    """1D template search, assume data and template are both 1d np array"""
    # demean
    data = data - np.mean(data)
    # detrend
    slope = (data[-1] - data[0]) / len(data)
    data -= slope * np.arange(len(data))
    # cross-correlate with the template
    corr = np.correlate(data, template, mode='valid')
    # find noise-level after correlation
    nl = 0.741 * stats.iqr(data)
    # find snr of the maxima
    return np.max(np.abs(corr)) / nl

#########################
# visualization related #
#########################

def view_glitches(tod, dets, cv, ncol=5, figsize=(20,20), ymax=1e-11, ymin=None):
    nrow = int(np.ceil(len(cv)/ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=figsize)
    snippets = get_glitch_snippets(tod, dets, cv)
    for i, sub_data in enumerate(snippets):
        r, c = i//ncol, i%ncol
        axes[r, c].plot(sub_data.T, 'k-', alpha=0.1)
        if ymin is None: ymin = -ymax
        axes[r, c].set_ylim([-ymax, ymax])
    return axes

def plot_snippet(snippet, demean=False, deslope=False, debuffer=0, **kwargs):
    from matplotlib import pyplot as plt
    opts = {
        'color': 'k',
        'ls': '-',
        'alpha': 0.3,
    }
    opts.update(kwargs)
    data = snippet.data
    if demean: data = data - np.mean(data,axis=1)[:,None]
    if deslope:
        slope = (data[:,-1] - data[:,0]) / data.shape[-1]
        data = data - slope[:,None] * np.arange(data.shape[-1])
    if debuffer > 0:
        # make sure we don't go over the limit
        debuffer = min(debuffer, data.shape[-1]//2-1)
        data = data[:,debuffer:-debuffer]
    plt.plot(data.T, **opts)
    plt.xlabel('samps')
    return plt.gca()

def plot_peaks_radius_binned(r, p, bins=10, method='linear'):
    from matplotlib import pyplot as plt
    bc, val, err = bin_data(r, p, bins, method)
    plt.errorbar(bc, val, yerr=err)
    plt.xlabel('distances from max peak')
    plt.ylabel('peak heights')
    return plt.gca()
