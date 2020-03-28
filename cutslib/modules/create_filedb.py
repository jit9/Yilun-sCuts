"""This module creates the filedb.txt file required by enki for
map-making.  It makes use of the filedb template stored in the
template folder and populate the tags based on season and array.

"""

import os.path as op, json
from jinja2 import Template
import cutslib
from cutslib.util import tag_to_afsv

class Module:
    def __init__(self, config):
        self.cut_release = config.get("cut_release", None)
        self.outfile = config.get("outfile", "filedb.txt")
        self.template = config.get("template", "default")

    def run(self, p):
        template = self.template
        cr = self.cut_release
        outfile = self.outfile

        # load template
        if template == 'default':
            # for default template we use the one in cutslib
            template = op.abspath(op.join(op.dirname(cutslib.__file__),
                                          '..', 'templates',
                                          'filedb_template.txt'))
        with open(template, "r") as f:
            tpl = Template(f.read())

        # load release
        cr_file = p.depot.get_deep((f'release_{cr}.txt',))
        with open(cr_file, "r") as f:
            rl = json.loads(f.read())

        # load tags
        tags = list(rl['tags'].keys())
        lines = []
        for tag in tags:
            tag_out = rl['tags'][tag]['tag_out']
            tag_partial = rl['tags'][tag]['tag_partial']
            # parse afsv
            ar, freq, season, cver = tag_to_afsv(tag_out)
            _, _, _, pver = tag_to_afsv(tag_partial)
            # create line
            line = f'if pa == {ar[-1]} and freq == "f{freq:03d}": '\
                   f'ctag, ptag = "{cver}", "{pver}"'
            lines.append(line)
        res = tpl.render(cut_tags="\n".join(lines))

        # writing output file
        # check output dir exists
        outdir = op.abspath(op.dirname(outfile))
        if not op.exists(outdir):
            os.makedirs(outdir)
            print("Creating: %s" % outdir)
        with open(outfile, "w") as f:
            # write down path
            f.write(res)
        print("Writing to: %s" % outfile)
