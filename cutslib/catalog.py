import fitsio, glob, os, numpy as np
import pandas as pd
import moby2


class Catalog():
    """Pandas based class for ACTPol catalog"""
    def __init__(self):
        # keep track of changes
        self.abuses = {}
        self.acqs = []
        self.data = None

    def load_data(self, filename=None):
        if not filename:
            filename = moby2.user_cfg.get_deep(('obs_catalog','filename'))
        self.filename = filename
        npcat = fitsio.read(filename)
        npcat = npcat.byteswap().newbyteorder()
        self.data = pd.DataFrame.from_records(npcat)
        self.data.index = pd.to_datetime(self.data.date)
        return self

    def load_acqs(self, array, bs_dir):
        """Load bias-step from a specific array
        Args:
            array: pa4,5,6 for example
            bs_dir: dir to bias-step files
        """
        # use .cal as an example
        files = glob.glob(os.path.join(bs_dir, array) + '/*/*.cal')
        # get basenames as ctimes (str)
        acqs = np.array([int(os.path.basename(f).split('.cal')[0]) \
                         for f in files])
        acqs.sort()
        ctime = self.data.ctime
        i_prev_bs = acqs.searchsorted(ctime)-1
        last_bs = acqs[i_prev_bs]
        self.data['prev_bs_ctime'] = last_bs
        self.query['acqs'] = 'loaded'
        # FIXME: need to ignore bs prior to the last IV
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
                sel *= self.data[k] == query[k]
        sel = sel.astype(bool)
        self.data = self.data[sel]
        self.abuses.update(query)
        return self

    def __repr__(self):
        return f"Catalog(n_entry={len(self.data)})"

    def __str__(self):
        repr = "Catalog:\n" \
               f"=> n_entry={len(self.data)}\n" \
               f"=> select={self.query}"
        return repr

    def __len__(self):
        return len(self.data)
