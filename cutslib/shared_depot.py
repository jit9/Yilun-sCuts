"""Module that wraps around shared_depot"""

import os, glob
import moby2

class SharedDepot:
    """Thin wrapper of the act_shared_depot"""
    def __init__(self, path=None):
        # if path not otherwise specified, use default in .moby2
        if not path:
            self.root = moby2.user_cfg.get_deep(('depots',
                                                 'actpol_shared',
                                                 'path'))
        else:
            self.root = path

    def get_deep(self, uri):
        return os.path.join(self.root, *uri)

    def list(self, uri=()):
        """List contents for interactive use
        Args:
            uri (tuple): path (with wildcard) of interests.
        Examples:
            sd.list('*')
            sd.list(('BiasStepTimes','*'))
        """
        path = os.path.join(self.root,*uri)
        return [os.path.basename(f) for f in glob.glob(path)]

    def __repr__(self):
        return f"SharedDepot(root={self.root})"
