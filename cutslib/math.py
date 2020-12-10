import numpy as np

def nzmean(x, mask=None, **kwargs):
    """mean of non-zero values, using masked array"""
    if mask is None: mask = (x==0)
    return np.mean(np.ma.array(x, mask=mask), **kwargs)
