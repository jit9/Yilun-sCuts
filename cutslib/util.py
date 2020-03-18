"""This script contains common utility function for all other scripts"""

import os, pickle
import moby2

def tag_to_afsv(tag, ar=True):
    """Parse array(a), freq(f), season(s) and version(v) from tag
    An example tag is pa4_f150_s17_c11_v0

    Args:
        tag (str): i.e. pa4_f150_s17_c11_v0
        ar (bool): whether season should be prefixed with ar or pa
    Returns:
        array, freq, season, version
    """
    array = tag.split('_')[0]
    freq = int(tag.split('_')[1][1:])  # from f150 -> int(150)
    season = '20'+tag.split('_')[2][1:]  # from s17 -> 2017
    version = tag.split('_')[-1]
    if array[0] == 'a':  # start with ar
        if not ar:       # but don't want ar
            array = 'pa'+array[-1]  # from ar4 -> pa4
    else:                # start with pa
        if ar:           # but want ar
            array = 'ar'+array[-1]
    return (array, freq, season, version)

def mkdir(dir):
    if not os.path.exists(dir):
        print("Creating directory: %s" % dir)
        os.makedirs(dir)
    return dir

def parse_tag(cutparam):
    params = parse_param(cutparam)
    return params.get('tag_out')

def parse_depot(cutparam):
    params = parse_param(cutparam)
    return params.get('depot')

def parse_param(cutparam):
    return moby2.util.MobyDict.from_file(cutparam)

def pickle_load(filename):
    """Load pickle file in a py2/3 compatible way"""
    with open(filename, "rb") as f:
        try:
            data = pickle.load(f)
        except UnicodeDecodeError:
            f.seek(0)  # fix 'cannot find MARK' bug
            data = pickle.load(f, encoding='latin1')
    return data
