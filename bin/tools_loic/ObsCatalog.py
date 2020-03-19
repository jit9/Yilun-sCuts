import fitsio
import pandas as pd
import moby2


class Catalog():
    """Pandas based class for ACTPol catalog"""
    def __init__(self, filename=filename):
        npcat = fitsio.read(filename)
        npcat = npcat.byteswap().newbyteorder()
        self.data = pd.DataFrame.from_records(npcat)
        self.data.index = pd.to_datetime(self.data.date)
