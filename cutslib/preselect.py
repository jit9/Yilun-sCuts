import numpy as np

class PreselectionError(RuntimeError): pass

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
    return presel_func


@preselector
def by_median(corrmat, min_corr=0.6, min_sel=10, min_frac=0.1):
    # select dets with median above `min_corr`
    sel = np.nanmedian(np.abs(corrmat), axis=1) > min_corr
    # check whether we have enough dets
    min_dets = max(min_sel, int(corrmat.shape[0]*min_frac))
    if np.sum(sel) < min_dets:
        raise PreselectionError
    return sel
