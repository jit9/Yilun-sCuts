"""Detectors related recipes"""

def get_ff_unstable(cpar):
    """get the unstable detectors from the flatfield
    Args:
        cpar: path to the cutparams file
    Example:
        cuts dets get_ff_unstable cutparams_v3.par
    """
    from moby2.util import MobyDict
    import numpy as np, sys
    cPar = cpar.replace("cutp","cutP")
    ff_file = MobyDict.from_file(cPar).get_deep(('pathologyParams',
                                                 'calibration',
                                                 'flatfield'))
    ff = MobyDict.from_file(ff_file)
    idx = np.where(np.array(ff['stable']) == False)[0]
    det_uid = list(np.array(ff['det_uid'])[idx])
    print("det_uid=%s" % det_uid)

def union(*dicts):
    """merge det_uid in different files into one dict file by union
    Args:
        dicts: dict files
    Example:
        cuts dets union exclude*.dict"""
    from moby2.util import MobyDict
    all_dets = set()
    for l in dicts:
        all_dets = all_dets.union(set(MobyDict.from_file(l)['det_uid']))
    all_dets = list(all_dets)
    all_dets.sort()
    print("det_uid = %s" % all_dets)

def intersection(*dicts):
    """merge det_uid in different files into one dict file by intersection
    Args:
        dicts: dict files
    Example:
        cuts dets intersection exclude*.dict"""
    from moby2.util import MobyDict
    for i, l in enumerate(dicts):
        if i == 0:
            all_dets = set(MobyDict.from_file(l)['det_uid'])
        else:
            all_dets = all_dets.intersection(set(MobyDict.from_file(l)['det_uid']))
    all_dets = list(all_dets)
    all_dets.sort()
    print("det_uid = %s" % list(all_dets))
