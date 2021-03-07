"""convenient functions to load tods"""
import os.path as op
import numpy as np, h5py
import yaml
import moby2
from moby2.util import MobyDict
from . import util
from .depot import Depot, SharedDepot
from .pathologies import Pathologies, get_pathologies
from .pathologies_tools import get_pwv

glitchp = {'nSig': 10., 'tGlitch' : 0.007, 'minSeparation': 30, \
           'maxGlitch': 50000, 'highPassFc': 6.0, 'buffer': 200 }

def load_tod(todname, tag=None, planet=None, partial=None, depot=None,
             abscal='201026', verbose=True, release='20201005', rd=True, fs=False,
             fcode=None, autoloads=['cuts','planet','partial','patho','cal','abscal','pwv'],
             **kwargs):
    if ':' in todname: todname, fcode = todname.split(':')
    opts = {'filename': todname, 'repair_pointing':True, 'read_data': rd, 'fix_sign': fs}
    opts.update(kwargs)
    tod = moby2.scripting.get_tod(opts)
    if verbose: print("tod loaded")
    # load metadata with release tag
    if len(autoloads)>0 and release:
        depot = Depot(path=depot)
        release_file = depot.get_deep((f'release_{release}', 'release.txt'))
        with util.nowarn():
            with open(release_file, "r") as f:
                rl = yaml.load(f.read())
        if verbose: print(f"Loaded tags from {release_file}")
        pa = tod.info.array.replace('ar','pa')
        scode = tod.info.season.replace('20','s')
        if tag: fcode = tag.split('_')[1]
        if not fcode: raise ValueError("Missing fcode")
        key = f"{pa}_{fcode}_{scode}_c11"
        tags = rl['tags'][key]
        # try to load iv
        if 'iv' in autoloads:
            cal = moby2.scripting.get_calibration(
                {'type':'iv', 'source':'data'}, tod=tod)
            tod.cal = cal
        # try to load cuts
        if 'cuts' in autoloads:
            try:
                cuts = moby2.scripting.get_cuts({
                    'depot': depot.root,
                    'tag': tags['tag_out'],
                }, tod=tod)
                tod.cuts = cuts
                if verbose: print(f"-> cuts loaded in tod.cuts tagged: {tags['tag_out']}")
            except:
                print("Warning: cuts not loaded successfully")
        # load planet cuts
        if 'planet' in autoloads:
            try:
                cuts = moby2.scripting.get_cuts({
                    'depot': depot.root,
                    'tag': tags['tag_planet']
                }, tod=tod)
                tod.planet = cuts
                if verbose: print(f"-> planet cuts loaded in tod.planet tagged: {tags['tag_planet']}")
            except:
                print("Warning: planet cuts not loaded successfully")
        # load partial
        if 'partial' in autoloads:
            try:
                cuts = moby2.scripting.get_cuts({
                    'depot': depot.root,
                    'tag': tags['tag_partial']
                }, tod=tod)
                tod.partial = cuts
                if verbose: print(f"-> partial cuts loaded in tod.partial tagged: {tags['tag_partial']}")
            except:
                print("Warning: partial cuts not loaded successfully")
        # load calibrations
        if 'cal' in autoloads:
            try:
                cal = depot.read_object(moby2.Calibration, tod=tod, tag=tags['tag_cal'])
                tod.cal = cal
                if verbose: print(f"-> cal loaded in tod.cal tagged: {tags['tag_cal']}")
            except:
                print("Warning: calibrations not loaded successfully")
        # load pathologies
        if 'patho' in autoloads:
            try:
                patho = get_pathologies({
                    'depot': depot.root,
                    'tag': tags['tag_out']
                }, tod=tod)
                tod.patho = patho
                if verbose: print(f"-> patho loaded in tod.patho tagged: {tags['tag_out']}")
            except:
                print("Warning: patho not loaded successfully")
        # get pwv
        if 'pwv' in autoloads:
            try:
                pwv = get_pwv([tod.ctime[0]])
                tod.pwv = pwv[0]
                if verbose: print("-> pwv loaded in tod.pwv")
            except:
                print("Warning: pwv not loaded successfully")
        # get abscal
        if 'abscal' in autoloads:
            try:
                abscal_file = SharedDepot().get_deep(('TODAbsCal',f'abscal_{abscal}.h5'))
                with h5py.File(abscal_file, "r") as f:
                    abscal_data = f['abscal'][:]
                    bmask = abscal_data['band_id'].astype(str) == fcode
                    todmask = abscal_data['tod_id'].astype(str) == tod.info.name
                    cal = abscal_data['cal'][bmask*todmask][0]
                    del abscal_data, bmask
                    tod.abscal = cal
                if verbose: print(f"-> abscal loaded in tod.abscal tagged: {abscal}")
            except:
                print("Warning: abscal not loaded successfully")
    elif tag:
        # try to load cuts
        try:
            cuts = moby2.scripting.get_cuts({
                'depot': depot.root,
                'tag': tag,
            }, tod=tod)
            tod.cuts = cuts
            print("-> cuts loaded in tod.cuts")
        except:
            print("Warning: cuts not loaded successfully")
        # load planet cuts
        if not planet: planet = tag+'_planet'
        try:
            cuts = moby2.scripting.get_cuts({
                'depot': depot.root,
                'tag': planet
            }, tod=tod)
            tod.planet = cuts
            print("-> planet cuts loaded in tod.planet")
        except:
            print("Warning: planet cuts not loaded successfully")
        # load planet cuts
        if not partial: partial = tag+'_partial'
        try:
            cuts = moby2.scripting.get_cuts({
                'depot': depot.root,
                'tag': partial
            }, tod=tod)
            tod.partial = cuts
            print("-> partial cuts loaded in tod.partial")
        except:
            print("Warning: partial cuts not loaded successfully")
        # load calibrations
        try:
            cal = depot.read_object(moby2.Calibration, tod=tod, tag=tag)
            tod.cal = cal
            print("-> cal loaded in tod.cal")
        except:
            print("Warning: calibrations not loaded successfully")
        # load pathologies
        try:
            patho = get_pathologies({
                'depot': depot.root,
                'tag': tag
            }, tod=tod)
            tod.patho = patho
            print("-> patho loaded in tod.patho")
        except:
            print("Warning: patho not loaded successfully")
        # get pwv
        try:
            pwv = get_pwv([tod.ctime[0]])
            tod.pwv = pwv[0]
            print("-> pwv loaded in tod.pwv")
        except:
            print("Warning: pwv not loaded successfully")
    return tod

def get_tes(tod, mask=True):
    """return tes dets mask"""
    if mask: return tod.info.array_data['det_type'] == 'tes'
    else: return np.where(tod.info.array_data['det_type'] == 'tes')[0]

def quick_transform(tod, steps=[], safe=False, glitchp=glitchp, verbose=True,
                    opts={}):
    """Short-hand notation for some common processings. Steps specified
    will be executed in order. Here are some possible steps:
      detrend: detrend tod for each det
      demean: remove the mean of the tod for each det
      demean_nospike: remove the mean with glitches masked
      cal: calibrate tod with tod.cal
      abscal: abs calibration using tod.abscal
      fill_cuts: fill cuts using tod.cuts
      ff_mce: find and fix mce
      ff_glitch: find and fix glitches
      get_iv: get iv calibration value
    """
    for step in steps:
        if verbose: print(f"step: {step}")
        if step == 'detrend':
            moby2.tod.detrend_tod(tod)
        elif step == 'demean':
            tod.data -= np.mean(tod.data, axis=1)[:,None]
        elif step == 'demean_nospike':
            if step in opts:
                # allow specifying externally
                cuts = opts[step].get('cuts',None)
            if cuts is None:
                if not hasattr(tod, 'cuts'):
                    print(f"tod.cuts missing, skipping step: {step}")
                    if not safe: continue
                cuts = tod.cuts
            demean_nospike(tod, cuts)
        elif step == 'cal':
            if not hasattr(tod, 'cal'):
                print(f"tod.cal missing, skipping step: {step}")
                if not safe: continue
            tod.data[tod.cal.det_uid,:] *= tod.cal.cal[:,None]
        elif step == 'abscal':
            if not hasattr(tod, 'abscal'):
                print("tod.abscal missing, skipping step: {step}")
            tod.data *= tod.abscal
        elif step == 'fill_cuts':
            if step in opts:
                # allow specifying externally
                cuts = opts[step].get('cuts',None)
            if cuts is None:
                if not hasattr(tod, 'cuts'):
                    print(f"tod.cuts missing, skipping step: {step}")
                    if not safe: continue
                cuts = tod.cuts
            moby2.tod.cuts.fill_cuts(tod, cuts)
        elif step == 'filter_gain':  # remove filter gain
            moby2.tod.remove_filter_gain(tod)
        elif step == 'ff_mce':  # find and fill mce cuts
            mce_cuts = moby2.tod.get_mce_cuts(tod)
            tod.mce_cuts = mce_cuts
            moby2.tod.fill_cuts(tod, mce_cuts, no_noise=True)
        elif step == 'ff_glitch':  # find and fill glitch cuts
            pcuts = moby2.tod.get_glitch_cuts(tod=tod, params=glitchp)
            moby2.tod.fill_cuts(tod, pcuts, no_noise=True)
            tod.pcuts = pcuts
        elif step == 'get_iv':
            cal = moby2.scripting.get_calibration({'type':'iv', 'source':'data'}, tod=tod)
            tod.cal = cal
        else:
            raise NotImplementedError

def demean_nospike(tod, cuts):
    """Remove the mean of tod without using the glitches"""
    mask = np.stack([c.get_mask() for c in cuts.cuts], axis=0)
    assert tod.nsamps >= cuts.nsamps
    data_masked = np.ma.masked_array(tod.data[:,:cuts.nsamps], mask)
    del mask
    mean = np.mean(data_masked, axis=1)
    del data_masked
    tod.data -= mean[:,None]
    del mean
    return tod

# IO related

def get_partial(tod, params=None, depot=None, tag=None, force=False, write=True):
    """get precomputed partial cuts if it exists, otherwise compute it on
    the run"""
    if not depot: depot = Depot()
    exists = op.exists(depot.get_full_path(moby2.TODCuts, tag=tag, tod=tod))
    # if we want to skip creating partial cuts, load from depot
    if exists and not force:
        cuts = depot.read_object(moby2.TODCuts, tag=tag, tod=tod)
    else:
        # Generate and save new glitch cuts
        # note calbol may not be implemented...
        cuts = moby2.tod.get_glitch_cuts(tod=tod, params=params)
        # find mce cuts
        mce_cuts = moby2.tod.get_mce_cuts(tod)
        # merge it with the partial cuts
        cuts.merge_tod_cuts(mce_cuts)
        # write to depot, not needed here
        if write: depot.write_object(cuts, tag=tag, tod=tod)
    return cuts

def get_ff(tod, flatfield):
    """return flatfield readings"""
    ff_ = moby2.detectors.RelCal.from_dict(flatfield)
    ff_sel, ff = ff_.get_property('cal', tod.det_uid, default=0)
    _, stable  = ff_.get_property('stable', tod.det_uid, default=False)
    return ff, ff_sel, stable

def get_resp(tod, param):
    """return responsivity measurement"""
    resp = moby2.scripting.products.get_calibration(param, tod.info).cal
    resp_sel = resp != 0
    return resp, resp_sel

def get_dets(params):
    live = MobyDict.from_file(params['live'])['det_uid']
    dark = MobyDict.from_file(params['dark'])['det_uid']
    excl = MobyDict.from_file(params['exclude'])['det_uid']
    return live, dark, excl

