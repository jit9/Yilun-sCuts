"""This script contains common utility function for all other scripts"""

import os, pickle, numpy as np
import moby2

def tag_to_afsv(tag, ar=True):
    """Parse array(a), freq(f), season(s) and version(v) from tag
    An example tag is pa4_f150_s17_c11_v0

    Args:
        tag (str): i.e. pa4_f150_s17_c11_v0
        ar (bool): whether season should be prefixed with ar or pa
    Returns:
        array, freq, season, version
    """
    array = tag.split('_')[0]
    freq = int(tag.split('_')[1][1:])  # from f150 -> int(150)
    season = '20'+tag.split('_')[2][1:]  # from s17 -> 2017
    version = tag.split('_')[-1]
    # if tag has postfix like _partial, use the one before that
    if version[0] != 'v':
        version = tag.split('_')[-2]
    if array[0] == 'a':  # start with ar
        if not ar:       # but don't want ar
            array = 'pa'+array[-1]  # from ar4 -> pa4
    else:                # start with pa
        if ar:           # but want ar
            array = 'ar'+array[-1]
    return (array, freq, season, version)


def mkdir(dir, comm=None, rank=0):
    if comm: return mkdir(dir, None, comm.Get_rank())
    if rank == 0:
        if not os.path.exists(dir):
            print("Creating directory: %s" % dir)
            os.makedirs(dir)
        return dir

def parse_tag(cutparam):
    params = parse_param(cutparam)
    return params.get('tag_out')

def parse_depot(cutparam):
    params = parse_param(cutparam)
    return params.get('depot')

def parse_param(cutparam):
    return moby2.util.MobyDict.from_file(cutparam)

def pickle_load(filename):
    """Load pickle file in a py2/3 compatible way"""
    with open(filename, "rb") as f:
        try:
            data = pickle.load(f)
        except UnicodeDecodeError:
            f.seek(0)  # fix 'cannot find MARK' bug
            data = pickle.load(f, encoding='latin1')
    return data

def to_scode(season):
    season = str(season)
    if season=='2016':
        return 's16'
    elif season=='2017':
        return 's17'
    elif season=='2018':
        return 's18'
    elif season=='2019':
        return 's19'
    else:
        return season

def to_season(scode):
    return scode.replace('s','20')

def to_pa(array):
    return array.lower().replace('ar','pa')

def to_ar(array):
    return array.lower().replace('pa','ar')

def pas_per_season(season):
    if season == 's13':
        return ['pa1']
    elif season == 's14':
        return ['pa1', 'pa2']
    elif season == 's15':
        return ['pa1', 'pa2', 'pa3']
    elif season == 's16':
        return ['pa2', 'pa3', 'pa4']
    elif season in ('s17', 's18', 's19'):
        return ['pa4', 'pa5', 'pa6']
    else:
        raise ValueError('No arrays found for season {}'.
                         format(season))

def freqs_per_pa(pa):
    if pa in ('pa1', 'pa2'):
        return ['f150']
    elif pa in ('pa3', 'pa5', 'pa6'):
        return ['f090', 'f150']
    elif pa == 'pa4':
        return ['f150', 'f220']
    else:
        raise ValueError('No freqs found for season {}'.
                         format(pa))

# functions copied from pixell.utils to reduce dependency
def cumsum(a, endpoint=False):
	"""As numpy.cumsum for a 1d array a, but starts from 0. If endpoint is True, the result
	will have one more element than the input, and the last element will be the sum of the
	array. Otherwise (the default), it will have the same length as the array, and the last
	element will be the sum of the first n-1 elements."""
	res = np.concatenate([[0],np.cumsum(a)])
	return res if endpoint else res[:-1]

def moveaxis(a, o, n):
	if o < 0: o = o+a.ndim
	if n < 0: n = n+a.ndim
	if n <= o: return np.rollaxis(a, o, n)
	else: return np.rollaxis(a, o, n+1)

def allgather(a, comm):
	"""Convenience wrapper for Allgather that returns the result
	rather than needing an output argument."""
	a   = np.asarray(a)
	res = np.zeros((comm.size,)+a.shape,dtype=a.dtype)
	if np.issubdtype(a.dtype, np.string_):
		comm.Allgather(a.view(dtype=np.uint8), res.view(dtype=np.uint8))
	else:
		comm.Allgather(a, res)
	return res

def allgatherv(a, comm, axis=0):
	"""Perform an mpi allgatherv along the specified axis of the array
	a, returning an array with the individual process arrays concatenated
	along that dimension. For example gatherv([[1,2]],comm) on one task
	and gatherv([[3,4],[5,6]],comm) on another task results in
	[[1,2],[3,4],[5,6]] for both tasks."""
	a  = np.asarray(a)
	fa = moveaxis(a, axis, 0)
	# mpi4py doesn't handle all types. But why not just do this
	# for everything?
	must_fix = np.issubdtype(a.dtype, np.str_) or a.dtype == bool
	if must_fix:
		fa = fa.view(dtype=np.uint8)
	ra = fa.reshape(fa.shape[0],-1) if fa.size > 0 else fa.reshape(0,np.product(fa.shape[1:],dtype=int))
	N  = ra.shape[1]
	n  = allgather([len(ra)],comm)
	o  = cumsum(n)
	rb = np.zeros((np.sum(n),N),dtype=ra.dtype)
	comm.Allgatherv(ra, (rb, (n*N,o*N)))
	fb = rb.reshape((rb.shape[0],)+fa.shape[1:])
	# Restore original data type
	if must_fix:
		fb = fb.view(dtype=a.dtype)
	return moveaxis(fb, 0, axis)

def decode_array_if_necessary(arr):
    """Given an arbitrary numpy array arr, decode it if it is of type S
    and we're in a version of python that doesn't like that

    """
    try:
        np.array(["a"],"S")[0] in "a"
        return arr
    except TypeError:
        if arr.dtype.type is np.bytes_:
            return np.char.decode(arr)
        else:
            return arr

def encode_array_if_necessary(arr):
    """Given an arbitrary numpy array arr, encode it it if it is of type U
    and we're in a version of python that doesn't like that
    """
    try:
        np.array(["a"],"S")[0] in "a"
        return arr
    except TypeError:
        if arr.dtype.type is np.str_:
            return np.char.encode(arr)
        else:
            return arr
