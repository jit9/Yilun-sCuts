"""This module initialize the cutparams with a default template

Options
-------
season: season in scode format
arrays: list of arrays seperated by space
cut_software: version tag of software, default to c11

This module will generate freqs for all freqs for the arrays listed
in the specified season, if arrays are not provided, it will use
the default list.

"""
import os, os.path as op
from jinja2 import Template

import cutslib
from cutslib.util import pas_per_season, freqs_per_pa
from cutslib.util import to_ar, to_pa, to_scode, to_season


class Module:
    def __init__(self, config):
        self.season = to_scode(config.get("season", None))
        arrays = config.get("arrays", None)
        if not arrays:
            self.arrays = pas_per_season(self.season)
        else:
            self.arrays = arrays.split()
        self.cut_software = config.get("cut_software", "c11")
        self.source_scan = config.get("source_scan")
        self.dets_live = config.get("dets_live")
        self.dets_exclude = config.get("dets_exclude")
        self.dets_dark = config.get("dets_dark")
        self.flatfield = config.get("flatfield")
        self.pointing = config.get("pointing")
        self.tod_offsets = config.get("tod_offsets")
        self.biasstep = config.get("biasstep")
        self.force = config.getboolean("force", False)

    def run(self, p):
        season = self.season
        softver = self.cut_software
        # loop over arrays and freqs
        for ar in self.arrays:
            # make sure our format is correct
            pa = to_pa(ar)
            freqs = freqs_per_pa(pa)
            for freq in freqs:
                # a dictionary to be used to autofill variables
                # in the configuration
                autofill = {
                    'season': to_season(season),
                    'ar': to_ar(pa),
                    'pa': pa,
                    'scode': to_scode(season),
                    'shared_depot': p.shared_depot.root,
                    'freq': freq,
                }
                outdir = f"{pa}_{freq}_{season}_{softver}"
                if not op.exists(outdir):
                    print(f"Creating: {outdir}")
                    os.makedirs(outdir)
                # load cutparam template and render it
                template = op.abspath(op.join(op.dirname(cutslib.__file__),
                                              '..', 'templates',
                                              'cutparams_v0.par'))
                with open(template, "r") as f:
                    tpl = Template(f.read())
                # get source scan by making a copy of the template
                # and fill up relevant terms
                source_scan = self.source_scan.format(**autofill)
                # generate slurmjob name like p5s180 which stands
                # for pa5 s18 f090
                jobname = pa[0]+pa[2]+season+freq[1]
                # render template
                res = tpl.render(
                    depot=p.depot.root,
                    tag_out=f"{pa}_{freq}_{season}_{softver}_v0",
                    source_scan=source_scan,
                    jobname=jobname
                )
                # generate cutparams in each folder
                outfile = op.join(outdir, 'cutparams_v0.par')
                if op.exists(outfile):
                    if not self.force:
                        ans = input(f"{outfile} exists, overwrite? y/n ")
                        to_write = ans == 'y'
                    else:
                        to_write = True
                    if to_write:
                        with open(outfile, "w") as f:
                            f.write(res)
                        print(f"Written to: {outfile}")

                # load cutParams template and render it
                template = op.abspath(op.join(op.dirname(cutslib.__file__),
                                              '..', 'templates',
                                              'cutParams_v0.par'))
                with open(template, "r") as f:
                    tpl = Template(f.read())
                pointing = self.pointing.format(**autofill)
                tod_offsets = self.tod_offsets.format(**autofill)
                dets_live = self.dets_live.format(**autofill)
                dets_dark = self.dets_dark.format(**autofill)
                dets_exclude = self.dets_exclude.format(**autofill)
                flatfield = self.flatfield.format(**autofill)
                biasstep = self.biasstep.format(**autofill)
                res = tpl.render(
                    depot=p.depot.root,
                    pointing=pointing,
                    tod_offsets=tod_offsets,
                    dets_live=dets_live,
                    dets_dark=dets_dark,
                    dets_exclude=dets_exclude,
                    flatfield=flatfield,
                    biasstep_tag=biasstep,
                )
                # generate cutparams in each folder
                outfile = op.join(outdir, 'cutParams_v0.par')
                if op.exists(outfile):
                    if not self.force:
                        ans = input(f"{outfile} exists, overwrite? y/n ")
                        to_write = ans == 'y'
                    else:
                        to_write = True
                    if to_write:
                        with open(outfile, "w") as f:
                            f.write(res)
                        print(f"Written to: {outfile}")
