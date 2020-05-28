"""Utility functions related to glitch analysis"""
import numpy as np

class TODSnippet:
    def __init__(self, tod, det_uid, tslice):
        """Create a TOD snippet with a subset of dets and a slice of samples

        Parameters
        ----------
        tod: base TOD object
        det_uid: list of dets of interests
        tslice: a slice object to get a subset of samples

        """
        self.det_uid = det_uid
        self.data = tod.data[det_uid,tslice]
        self.array_data = tod.info.array_data
        self.tslice = tslice
    def demean(self):
        """Remove the mean of the snippet"""
        self._mean = self.data.mean(axis=1)
        self.data -= self._mean[:,None]
        return self
    def __repr__(self):
        return f"TODSnippet(ndet={len(self.det_uid)},tslice={self.tslice})"

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
    assert tod is not None, "tod is needed to get metadata!"
    # first make sure snippet has a mean of zero
    snippet = snippet - snippet.mean(axis=1)[:,None]  # not in-place
    # find peak height
    peaks = np.max(np.abs(snippet), axis=1)
    # find highest peak as the reference point
    imax = np.argmax(peaks)
    x = tod.info.array_data[f'{ref}_x'][dets]
    y = tod.info.array_data[f'{ref}_y'][dets]
    x_c = x[imax]
    y_c = y[imax]
    r = np.sqrt((x-x_c)**2+(y-y_c)**2)
    return r, peaks

def glitch_det_count(cuts):
    """Count number of dets affected as a time series

    Parameters
    ----------
    cuts: TODCuts object

    """
    return np.sum([c.get_mask() for c in cuts.cuts], axis=0)

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

class PixelReader:
    def __init__(self, season='2016', array='AR3', mask=None):
        """Utility class to find information about each pixel (feedhorn), migrated
        from todloop.util.pixels.PixelReader to manage internally"""
        import moby2
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

#############
# utilities #
#############

def bin_data(x, y, bins=10, method='linear'):
    """calculate binned statistics"""
    if method == 'linear': x_ = x
    elif method == 'log': x_ = np.log(x)
    elif method == 'p2': x_ = x**0.5
    else: raise ValueError("Unsupported binning method")
    res = stats.binned_statistic(x_, y, bins=bins, statistic='mean')
    res_std = stats.binned_statistic(x_, y, bins=bins, statistic='std')
    bc_ = (res.bin_edges[1:] + res.bin_edges[:-1])/2
    if method == 'linear': bc = bc_
    elif method == 'log': bc = np.exp(bc_)
    elif method == 'p2': bc = bc_**2
    else: raise ValueError("Unsupported binning method")
    return bc, res.statistic, res_std.statistic

def cv2slices(cv):
    """Convert CutsVector to a list ot time slices"""
    return [slice(v[0],v[1],None) for v in cv]

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
