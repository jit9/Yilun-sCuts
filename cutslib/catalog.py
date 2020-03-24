import fitsio, os, numpy as np
import pandas as pd
import moby2
from cutslib import SharedDepot

class Catalog():
    """Pandas based class for ACTPol catalog"""
    def __init__(self, load_data=True, filename=None):
        # keep track of changes
        self.abuses = {'data': 'not loaded'}
        self.acqs = []
        self.data = None
        if load_data:
            self.load_data(filename)

    def load_data(self, filename=None):
        if not filename:
            filename = moby2.user_cfg.get_deep(('obs_catalog','filename'))
        self.filename = filename
        npcat = fitsio.read(filename)
        npcat = npcat.byteswap().newbyteorder()
        self.data = pd.DataFrame.from_records(npcat)
        self.data.index = pd.to_datetime(self.data.date)
        self.abuses.update({'data':'loaded'})
        return self

    def load_acqs(self, season=None, array=None, tag='171110'):
        """Load bias-step from a specific season/array/tag
        Args:
            season: scode format, numeric format will be auto-converted
            array: pa format, ar format will be aut-converted
            tag: postfix of the bias_step files
        """
        if (not season) or (not array):
            raise ValueError("Season and array needs to be specified!")
        # in case season and array are not in the expected format
        array = array.lower().replace('ar','pa')
        season = str(season)
        if season=='2016':
            season = 's16'
        elif season=='2017':
            season = 's17'
        elif season=='2018':
            season = 's18'
        elif season=='2019':
            season = 's19'
        else:
            raise ValueError("Unknown season: %s" % season)
        # load bias step db filepath
        sd = SharedDepot()
        fpath = sd.get_deep(('BiasStepTimes',
                             f'intervals_{season}_{array}_{tag}.txt'))
        df = pd.read_fwf(fpath)
        # fix column name (remove # in the first column name)
        columns = list(df.columns)
        columns[0] = columns[0].replace('# ', '')
        df.columns = columns
        # if biasstep_tag is not in the catalog, add it
        if 'bs_tag' not in self.data.columns:
            self.data['bs_tag'] = ['']*len(self.data)
        # fill biasstep_tag
        for i, r in self.data.iterrows():
            # search for bias_step
            ctime = r['ctime']
            res = df[np.logical_and(df.ctime0 <= ctime, df.ctime1 > ctime)]
            if len(res) == 0:
                continue
            elif len(res) == 1:
                self.data.loc[i,'bs_tag'] = res['biasstep_tag'].values[0]
            else:
                raise ValueError("Unexpected happens!")

        return self

    def select(self, query={}):
        """Select a subset of catalog
        Args:
            query (dict): a dict containing the conditions
              each key should match a column name in catalog
        """
        # TODO: this is only equality based now, i should support
        # range based as well
        sel = np.ones(len(self.data))
        for k in query:
            if k in self.data.columns:
                if k == 'array':
                    target = query[k].replace("pa","ar")
                else:
                    target = query[k]
                sel *= self.data[k] == target
        sel = sel.astype(bool)
        self.data = self.data[sel]
        self.abuses.update(query)
        return self

    def __repr__(self):
        return f"Catalog(n_entry={len(self.data)})"

    def __str__(self):
        repr = "Catalog:\n" \
               f"=> n_entry={len(self.data)}\n" \
               f'=> abuses={str(self.abuses).replace(" ","")}'
        return repr

    def __len__(self):
        return len(self.data)
