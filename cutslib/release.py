"""Utility function for working with release file"""
import yaml

from .depot import Depot
from .util import tod_to_tag, tag_to_afsv, get_tod_fcode, nowarn

def ver(tag):
    splits = tag.split("_")
    if len(splits) == 6: return splits[-2]
    else: return splits[-1]

class Release:
    def __init__(self, release_file=None, tag='20201005', depot=None):
        """Load release file

        Parameters
        ----------
        release_file: path to release file
        tag: if no release file is provided, it will be loaded based on release tag
        depot: depot to load release from if non-default one is to be used

        """
        self.tag = tag
        # load release file
        if release_file is None:
            depot = Depot(depot)
            release_file = depot.get_deep((f'release_{tag}', 'release.txt'))
        with nowarn():
            with open(release_file, "r") as f:
                self.release = yaml.load(f.read())
    def tags2filedb(self):
        """Convert tags into strings recognized by filedb"""
        text = "{yilun} = " + f"'{self.release['depot']}'\n"
        for dset, tags in self.release['tags'].items():
            pa, fcode, scode, _ = dset.split('_')
            ptag = ver(tags['tag_partial'])
            ctag = ver(tags['tag_out'])
            gtag = ver(tags['tag_cal'])
            text += f"if season == {scode[1:]} and pa == {pa[2:]} and freq == '{fcode}': " + \
                    f"ctag, ptag, gtag = '{ctag}', '{ptag}', '{gtag}'\n"
        text += 'cut_quality = "{yilun}/TODCuts/pa{pa}_{freq}_s{season}_c11_{ctag}/{t5}/{id}.cuts"\n'
        text += 'cut_basic   = "{yilun}/TODCuts/pa{pa}_{freq}_s{season}_c11_{ptag}_partial/{t5}/{id}.cuts"\n'
        text += 'gain        = "{yilun}/Calibration/pa{pa}_{freq}_s{season}_c11_{gtag}/{t5}/{id}.cal"\n'
        return text

    def tags(self, key):
        return self.release['tags'][key]

    def tod_tags(self, tod):
        key = tod_to_tag(tod)
        return self.tags(key)

    def tod_pcuts(self, tod):
        tag = self.tod_tags(tod)['tag_partial']
        depot = self.release['depot']
        tod, _ = get_tod_fcode(tod)
        return f"{depot}/TODCuts/{tag}/{tod[:5]}/{tod}.cuts"

    def tod_cuts(self, tod):
        tag = self.tod_tags(tod)['tag_out']
        depot = self.release['depot']
        tod, _ = get_tod_fcode(tod)
        return f"{depot}/TODCuts/{tag}/{tod[:5]}/{tod}.cuts"

    def tod_cal(self, tod):
        tag = self.tod_tags(tod)['tag_cal']
        depot = self.release['depot']
        tod, _ = get_tod_fcode(tod)
        return f"{depot}/Calibration/{tag}/{tod[:5]}/{tod}.cal"
