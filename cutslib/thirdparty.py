"""This module contains overwrite of existing parts that
extends the original functionaility of cuts pipeline"""

import moby2
from todloop import Routine

from .pathologies_tools import reportPathologies


class reportPathologiesMod(reportPathologies):
    def __init__(self, params):
        """Extends the original class to support
        loading dets cut from a specified place"""
        reportPathologies.__init__(self, params)

    def get_det_cuts(self, pa, pathop):
        """Load cuts from a given tag instead of generating
        a new one based on pathology object"""
        # get cuts from depot
        depot = moby2.util.Depot(self.params.get('depot'))
        tag_cuts = self.params.get('thirdparty_tag_cuts')
        c_obj = depot.read_object(moby2.TODCuts,
                                  tag=tag_cuts,
                                  tod=pa.tod)
        return c_obj

class PathologyReportMod(Routine):
    def __init__(self, **params):
        """This routine aims to generate a Pathology Report object
        as is done for the moby2 script"""
        Routine.__init__(self)
        self.cutparam = params.get("cutparam", None)
        self.report = None

    def initialize(self):
        self.report = reportPathologiesMod(self.cutparam)
        # update the depot file name with rank if needed
        # if mpi, label depot file with the rank to avoid racing condition
        self.report.depot_file += '.%d' % self.get_rank()

    def execute(self, store):
        # get obs name
        obs = self.get_name()
        # append results
        self.report.appendResult(obs)
