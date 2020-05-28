"""This module generates a gain correction required by enki,
at this moment i will fill the list with default values"""

from cutslib import Catalog
import json

class Module:
    def __init__(self, config):
        self.sel = json.loads(config.get('sel'))
        self.outfile = config.get("outfile")
        self.fill_value = config.getfloat("fill_value", None)

    def run(self, p):
        cat = Catalog()
        cat.select(self.sel)
        with open(self.outfile, "w") as f:
            f.write(f"#                  tod_id band_id          cal\n")
            for i, r in cat.data.iterrows():
                todname = r.tod_name
                if 'ar4' in todname:
                    freqs = ['f150', 'f220']
                elif 'ar5' in todname:
                    freqs = ['f090', 'f150']
                else:
                    freqs = ['f090', 'f150']
                # get calibration
                if self.fill_value:
                    cal = self.fill_value
                for freq in freqs:
                    f.write(f"{todname} {freq}  {cal}\n")
        print(f"Written: {self.outfile}")
