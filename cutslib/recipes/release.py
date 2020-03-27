"""Relevant recipes for generating releases that can be
digested by the mapping crew
"""

def tags(*tags):
    """get release for input tags"""
    import os.path as op, glob
    import moby2
    lookup = {}
    # for each tag, find the latest cuts version
    for tag in tags:
        # find latest version
        cutparams = glob.glob(op.join(tag, 'cutparams_v*'))
        ver_latest = max([int(c.split('.par')[0][-1]) for c in cutparams])
        cpar_name = "cutparams_v{}.par".format("ver_latest")
        cpar_path = op.join(tag, cpar_name)
        # load cutparam
        cpar = moby2.util.MobyDict.from_file(cpar_path)
        lookup[tag] = {
            'tag_out': cpar.get('tag_out'),
            'tag_cal': cpar.get('tag_cal'),
            'tag_partial': cpar.get('tag_partial'),
            'tag_planet': cpar.get('tag_planet'),
            'tag_source': cpar.get('tag_source'),
        }
    import json
    print(json.dumps(lookup, indent=2, sort_keys=True))
