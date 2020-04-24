"""Immitate things in SO"""
from functools import reduce


class Axis:
    def __init__(self, count):
        self.count = count

class PathologyManager:
    """Only to look like one"""
    def __init__(self):
        self._fields = {}
        self._axes = {
            'dets': Axis(0)
        }
        # self._meta = {}  # not axis related metadata
        self.crits = {}
    def __getitem__(self, index):
        if index in self._axes:
            return self._axes[index]
        elif index in self._fields:
            return self._fields[index]
        else: raise ValueError(f"No {index} found!")
    def __getattr___(self, index):
        return self.__getitem__(index)
    def add(self, field, values):
        self._fields[field] = values
        return self
    @classmethod
    def for_tod(cls, tod):
        pman = cls()
        pman._axes['dets'].count = tod.data.shape[0]
        return pman
    def add_crit(self, field, limits=[5, 95], method='rel'):
        """add a criteria to find cuts.

        Parameters
        ----------
        field: field name to apply crit on
        limits: lengh 2 list with lower and higher limits, if no limits
          is needed put None accordingly, i.e. [2, None]
        method: can be 'rel' or 'abs'. When 'rel', limits refer to percentile
          and when 'abs', limits refer to absolute values
        """
        # TODO: support metadata
        assert len(limits) == 2
        lo, hi = limits
        self.crits.update({field: {
            'lo': lo, 'hi': hi,
            'method': method
        }})
        return self
    def clear_crit(self, fields=None):
        """remove crit for a given list of fields, defaults to all"""
        if not field:
            fields = list(self._fields.keys())
        for f in fields:
            if f in self.crits:
                del self.crits[f]
        return self
    def apply_crit(self):
        """apply crit to get a cut"""
        crits = [crit for crit in self.crit if crit in self._fields]
        ms = []  # masks
        for f in crits:
            v = self._fields[f]
            lo, hi = self.crits[f]['lo'], self.crits[f]['hi']
            method = self.crits[f]['method']
            m = np.ones(self.dets.count, dtype=bool)
            if method = 'rel':
                if lo: lo = np.percentile(v, lo)
                if hi: hi = np.percentile(v, hi)
            if lo: m *= (v >= lo)  # open
            if hi: m *= (v < hi)   # close
            ms.append(m)
        self.masks = ms  # store this as it includes lots of info
        return reduce(np.logical_and, ms)
