import numpy as np
import inspect
from .errors import PreselectionError

# Decorator for the lazy one: convert preselectors to partial functions
# without the need to specify the corrmat parameter and give them reasonable
# operation rules. For example:
#   logical or:    selector1 + selector2
#   logical and:   selector1 * selector2
#   fallback rule: selector1 >> selector2
# self is always evaluated before other.
def preselector(fun):
    class presel_func:
        def __init__(self, *args, **kwargs):
            self.args = args
            self.kwargs = kwargs
        def __call__(self, cc):
            return fun(cc, *self.args, **self.kwargs)
        def __mul__(self, other):
            def fun_prod(cc):
                return self(cc)*other(cc)
            return preselector(fun_prod)()
        def __add__(self, other):
            def fun_prod(cc):
                try: presel = self(cc)
                except PreselectionError:
                    return other(cc)
                try: opresel = other(cc)
                except PreselectionError:
                    return presel
                return presel+opresel
            return preselector(fun_prod)()
        def __rshift__(self, other):
            def fun_prod(cc):
                try: return self(cc)
                except PreselectionError:
                    return other(cc)
            return preselector(fun_prod)()
    # get help from someone
    presel_func.__doc__ = fun.__doc__
    # get arguments from someone after removing the partial args (first 1)
    args = [v for k, v in inspect.signature(fun).parameters.items()][1:]
    presel_func.__signature__ = inspect.Signature(parameters=args)
    return presel_func


@preselector
def by_median(corrmat, min_corr=0.6, min_sel=10, min_frac=None, dets=None):
    """Preselect by median. The number of detectors needs to be larger
    than the maximum of min_sel and ndet*min_frac

    Parameters
    ----------
    min_corr: minimum correlation for preselected detectors
    min_sel: minimum number of preselected detectors before failing
    min_frac: minimum fraction of preselected detectors before failing
    dets: if we want to narrow down to look at a subset of dets only
    """
    if dets is None: dets = np.ones(corrmat.shape[0], dtype=bool)
    # select dets with median above `min_corr`
    sel = np.nanmedian(np.abs(corrmat)[np.ix_(dets,dets)], axis=1) > min_corr
    # check whether we have enough dets
    if min_frac is None: min_frac=min_sel/corrmat.shape[0]
    min_dets = max(min_sel, round(corrmat.shape[0]*min_frac))
    if np.sum(sel) < min_dets:
        raise PreselectionError(f"sel={np.sum(sel)}, min_dets={min_dets}")
    # massage index to have the shape of input det numbers
    presel = np.zeros(corrmat.shape[0], dtype=bool)
    presel[dets] = sel
    return presel

@preselector
def by_mask(corrmat, mask):
    """Preselect by mask. A simple wrapper for boolean mask

    Parameters
    ----------
    mask: boolean mask

    """
    return mask
