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


class Depot(SharedDepot):
    """Thin wrapper of the depot. Comparing with the Depot object in
    moby2.util.Depot, this is mainly to provide a navigational helper to
    avoid repeating os.path.join procedures

    """
    def __init__(self, path=None):
        # if path not otherwise specified, use default in environment
        if not path:
            self.root = os.environ.get("CUTS_DEPOT")
        else:
            self.root = path
        # compatibility with moby2 depot
        self.moby2_depot = moby2.util.Depot(self.root)

    def __repr__(self):
        return f"Depot(root={self.root})"

    def read_object(self, *args, **kwargs):
        return self.moby2_depot.read_object(*args, **kwargs)

    def write_object(self, *args, **kwargs):
        try:
            res = self.moby2_depot.write_object(*args, **kwargs)
        except FileExistsError:  # unlikely but possible race condition
            kwargs['make_dirs'] = False
            res = self.moby2_depot.write_object(*args, **kwargs)
        return res

    def get_full_path(self, *args, **kwargs):
        return self.moby2_depot.get_full_path(*args, **kwargs)
