import moby2
import database
from database import Database as DB
import datetime
import numpy as np

def get_previous_BS(ctime, array):
    """Find the ctime of the previous BS"""
    
    if not isinstance(ctime, int):
        ctime = np.asarray(ctime)
        year = datetime.date.fromtimestamp(np.median(ctime)).year
    else:
        year = datetime.date.fromtimestamp(ctime).year
    db = DB(year)
    db.load_acqs()
    if isinstance(array, str):
        db.acqs = db.acqs[db.acqs.array==array]
        time_ivs = db.acqs.ctime[db.acqs.suffix=='iv'].values
        time_ivs.sort()
        time_bs = db.acqs.ctime[db.acqs.suffix=='bc1_step'].values
        time_bs.sort()
        previous_iv = time_ivs.searchsorted(ctime) - 1
        previous_bs = time_bs.searchsorted(ctime) - 1
        last_bs = time_bs[previous_bs]
        if np.any(time_bs[previous_bs] < time_ivs[previous_iv]):
            print "No BS available between last IV and TOD for at least one TOD. Return 0 for these TODs."
            last_bs[time_bs[previous_bs] < time_ivs[previous_iv]] = 0
    else:
        arrays = np.unique(array)
        last_bs = np.zeros_like(ctime)
        for ar in arrays:
            ar_idx = array == ar
            sel = db.acqs.array == ar
            acqs = db.acqs[sel]
            time_ivs = acqs.ctime[acqs.suffix=='iv'].values
            time_ivs.sort()
            time_bs = acqs.ctime[acqs.suffix=='bc1_step'].values
            time_bs.sort()
            previous_iv = time_ivs.searchsorted(ctime[ar_idx]) - 1
            previous_bs = time_bs.searchsorted(ctime[ar_idx]) - 1
            last_bs[ar_idx] = time_bs[previous_bs]
            if np.any(time_bs[previous_bs] < time_ivs[previous_iv]):
                print "No BS available between last IV and TOD for at least one TOD. Return 0 for these TODs."
                last_bs[ar_idx][time_bs[previous_bs] < time_ivs[previous_iv]] = 0
            print ar, (last_bs == 0).sum()
                    

    return last_bs

def get_next_BS(ctime, array):
    """Find the ctime of the next BS"""
    
    if not isinstance(ctime, int):
        ctime = np.asarray(ctime)
        year = datetime.date.fromtimestamp(ctime[0]).year
    else:
        year = datetime.date.fromtimestamp(ctime).year
    db = DB(year)

    db.load_acqs()
    if isinstance(array, str):
        db.acqs = db.acqs[db.acqs.array==array]
        time_ivs = db.acqs.ctime[db.acqs.suffix=='iv'].values
        time_ivs.sort()
        time_bs = db.acqs.ctime[db.acqs.suffix=='bc1_step'].values
        time_bs.sort()
        next_iv = time_ivs.searchsorted(ctime)
        next_bs = time_bs.searchsorted(ctime)
        following_bs = time_bs[next_bs]
        if np.any(time_bs[next_bs] > time_ivs[next_iv]):
            print "No BS available between last IV and TOD for at least one TOD. Return 0 for these TODs."
            next_bs[time_bs[next_bs] > time_ivs[next_iv]] = 0
    else:
        arrays = np.unique(array)
        following_bs = np.zeros_like(ctime)
        for ar in arrays:
            ar_idx = array == ar
            sel = db.acqs.array == ar
            acqs = db.acqs[sel]
            time_ivs = acqs.ctime[acqs.suffix=='iv'].values
            time_ivs.sort()
            time_bs = acqs.ctime[acqs.suffix=='bc1_step'].values
            time_bs.sort()
            next_iv = time_ivs.searchsorted(ctime[ar_idx])
            next_iv[next_iv==time_bs.size] = 0
            next_bs = time_bs.searchsorted(ctime[ar_idx])
            next_bs[next_bs==time_bs.size] = 0
            following_bs[ar_idx] = time_bs[next_bs]
            if np.any(next_bs == time_bs.size):
                print "There are TODs after the last BS. Return 0 for these TODs."
                following_bs[ar_idx][next_bs==time_bs.size] = 0.
            # if np.any(time_bs[next_bs] > time_ivs[next_iv]):
            #     print "No BS available between next IV and TOD for at least one TOD. Return 0 for these TODs."
            #     following_bs[ar_idx][time_bs[next_bs] > time_ivs[next_iv]] = 0
            print ar, (following_bs == 0).sum()

    return following_bs

def get_previous_IV(ctime, array):
    """Find the ctime of the previous IV"""
    
    if not isinstance(ctime, int):
        ctime = np.asarray(ctime)
        year = datetime.date.fromtimestamp(ctime[0]).year
    else:
        year = datetime.date.fromtimestamp(ctime).year
    db = DB(year)

    db.load_acqs()
    db.acqs = db.acqs[db.acqs.array==array]
    time_ivs = db.acqs.ctime[db.acqs.suffix=='iv'].values
    time_ivs.sort()
    previous_iv = time_ivs.searchsorted(ctime) - 1
    return time_ivs[previous_iv]


def get_next_IV(ctime, array):
    """Find the ctime of the next IV"""
    
    if not isinstance(ctime, int):
        ctime = np.asarray(ctime)
        year = datetime.date.fromtimestamp(ctime[0]).year
    else:
        year = datetime.date.fromtimestamp(ctime).year
    db = DB(year)

    db.load_acqs()
    db.acqs = db.acqs[db.acqs.array==array]
    time_ivs = db.acqs.ctime[db.acqs.suffix=='iv'].values
    time_ivs.sort()
    previous_iv = time_ivs.searchsorted(ctime)
    return time_ivs[previous_iv]

def get_between(ctime_start, ctime_stop, array, suffix):
    """Find all acqs between two ctimes"""
     
    if not isinstance(ctime, int):
        ctime = np.asarray(ctime)
        year = datetime.date.fromtimestamp(ctime[0]).year
    else:
        year = datetime.date.fromtimestamp(ctime).year
    db = DB(year)

    db.load_acqs()
    sel_array = db.acqs.array == array
    db.acqs = db.axqs[sel_array]
    sel_time = np.logical_and(db.acqs.ctime_start > ctime_start, db.acqs.ctime_start<ctime_stop)
    sel_suffix = db.acqs.suffix == suffix
    
    ctimes = db.acqs.ctime_start[np.logical_and(sel_suffix, sel_ctime)]
    return ctimes
