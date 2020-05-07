"""Relevant recipes for generating releases that can be
digested by the mapping crew
"""

###############
# interactive #
###############

def tags(*tags):
    """get release for input tags"""
    import os.path as op, glob, datetime, os, shutil
    import moby2
    from cutslib import Depot
    release = {}
    # add meta information
    # date
    today = datetime.date.today()
    release['version'] = prompt("Version", today.strftime("%Y%m%d"))
    release['date'] = today.strftime("%b %d, %Y")
    release['tags'] = {}
    release['author'] = prompt("Author", "Yilun Guan")
    release['depot'] = os.environ.get("CUTS_DEPOT")
    # for each tag, find the latest cuts version and relevant tags
    reports = []
    for tag in tags:
        # find latest version
        cutparams = glob.glob(op.join(tag, 'cutparams_v*'))
        ver_latest = max([int(c.split('.par')[0][-1]) for c in cutparams])
        cpar_name = "cutparams_v{}.par".format(ver_latest)
        cpar_path = op.join(tag, cpar_name)
        # load cutparam
        cpar = moby2.util.MobyDict.from_file(cpar_path)
        release['tags'][tag] = {
            'params': op.abspath(cpar_path),
            'tag_out': cpar.get('tag_out'),
            'tag_cal': cpar.get('tag_cal'),
            'tag_partial': cpar.get('tag_partial'),
            'tag_planet': cpar.get('tag_planet'),
            'tag_source': cpar.get('tag_source'),
            'tag_cmb': cpar.get('tag_cmb'),
        }
        # find report
        tag_out = cpar.get('tag_out')
        report = Depot().get_deep(('Postprocess', tag_out,'report', f'{tag_out}.pdf'))
        reports.append(report)
    import json
    print(json.dumps(release, indent=2))
    # release dir
    release_dir = op.join(release['depot'], f"release_{release['version']}")
    # create it if doesn't exists
    outdir = prompt("Write to", release_dir)
    if not op.exists(outdir): os.makedirs(outdir)
    # write release file
    outfile = op.join(outdir, 'release.txt')
    with open(outfile, "w") as f:
        f.write(json.dumps(release, indent=2))
    print("Release written: %s" % outfile)
    # copy reports to the release dir
    for report in reports:
        if not op.exists(report): continue
        outfile = op.join(outdir, op.basename(report))
        shutil.copyfile(report, outfile)
        print("Report written: %s" % outfile)


#####################
# utility functions #
#####################

def prompt(prompt=None, default=None):
    """A wrapper for input to allow default"""
    if not default:
        res = input(prompt)
        return res
    else:
        res = input("%s: (default: %s) " % (prompt, default))
        if res == "":
            return default
        else:
            return res
