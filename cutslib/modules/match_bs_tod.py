"""This module matches the bias-step to TODs"""

class Module:
    def __init__(self, config):
        self.filename = config.get("obs_catalog", None)
        self.tag = config.get("tag")
        self.tag_out = config.get("tag_out", None)
        self.bs_ver = config.get("bs_ver", "171110")

    def run(self, p):
        filename = self.filename
        tag = self.tag
        tag_out = self.tag_out
        bs_ver = self.bs_ver
        # load obs catalog
        import os, shutil
        from cutslib import Catalog
        from cutslib.util import to_pa, to_scode
        cat = Catalog(filename=filename)
        # down-select data
        season, array = p.i.season, p.i.ar
        cat.select({'season': season, 'array': array})
        # load bias-steps
        cat.load_acqs(season=season,array=array,version=bs_ver)
        # start copying files
        for i,r in cat.data.iterrows():
            fin = os.path.join(p.depot, 'biasstep', p.i.season, tag) + \
                  "/" + r['bs_tag']+".cal"
            # unless otherwise specified, use tag as tag_out
            if not tag_out:
                tag_out = tag
            array = to_pa(array)
            scode = to_scode(season)
            fout = os.path.join(p.depot, 'Calibration',
                                f"{array}_{scode}_bs_{tag_out}",
                                r['tod_name'][:5], r['tod_name']+".cal")
            # check whether input file exists
            if not os.path.isfile(fin):
                print("Warning: %s not found" % fin)
                continue
            # check whether output dir exists
            dout = os.path.dirname(fout)
            if not os.path.exists(dout):
                print("Creating %s" % dout)
                os.makedirs(dout)
            # start copying
            print("Writing: %s" % os.path.basename(fout))
            shutil.copy(fin, fout)
        print("Done!")
