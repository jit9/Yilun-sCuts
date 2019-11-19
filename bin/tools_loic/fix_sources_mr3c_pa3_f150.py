import h5py
import os
import joblib
import multiprocessing
import numpy as np
import moby2
from moby2.scripting import products
import sys

args = sys.argv

# obs = args[2]
# source_list = '/home/lmaurin/actpol_data_shared/BrightSources/mr3c_merged_gt2mK.txt'

params = moby2.util.MobyDict.from_file(args[1])
params_fix = moby2.util.MobyDict.from_file(args[2])
numcores = params_fix['ncores']
from moby2.util.database import TODList
tl = TODList.from_file(params['source_scans'])


# def new_partial_cuts(obs, params, params_fix):
#     new_tag_source = params_fix["tag_source"]
#     new_tag_partial = params_fix["tag_partial"]
    
#     IV_wait = params.get("IV_wait",0)
#     rebias_wait = params.get("rebias_wait",0)
#     end = params.get("sample_end",None)
    
    
#     depot = moby2.util.Depot(params.get('depot'))
#     cutParams = params.get('cutParams')
#     cpar = moby2.util.MobyDict.from_file(cutParams)
#     cpar['source_cuts'] = params_fix['source_cuts']
#     pathop = cpar['pathologyParams']
#     glitchp = cpar['glitchParams']
    
#     # Fill cuts parameter                                                                        
#     no_noise = not(cpar.get("fillWithNoise",True))
    
#     user_config = moby2.util.get_user_config()
#     moby2.pointing.set_bulletin_A(params=user_config.get('bulletin_A_settings'))
    
    
#     # Find IV/rebias gap time                                                                
#     ct = int(obs.split("/")[-1].split(".")[0])
#     ctimes = (ct-IV_wait,ct)
#     from moby2.instruments import actpol as inst
#     try:
#         config_file = params.get('manifest_conf', None)
#         db = inst.TODDatabase(config_file=config_file)
#     except:
#         db = None
#     if (db is None) or (db.db is None):
#         offset_IV = 0
#         offset_rebias = 0
#     else:
#         recs = db.select_acqs(suffix='iv', ctime=ctimes)
#         if len(recs) > 0:
#             offset_IV = (IV_wait - ct + recs[-1].ctime)*400
#         else:
#             offset_IV = 0
#         ctimes = (ct-rebias_wait,ct)
#         recs = db.select_acqs(suffix='bc1_step', ctime=ctimes)
#         if len(recs) > 0:
#             offset_rebias = (rebias_wait - ct + recs[-1].ctime)*400
#         else:
#             offset_rebias = 0
#     offset = max( offset_IV, offset_rebias, params.get("offset",0)*400 )
    
    
#     params_loadtod = {
#         'filename':obs,
#         'fix_sign':True,
#         'start':offset,
#         'end':end,
#         'repair_pointing':True
#     }
#     params_loadtod.update(params.get("params_loadtod",{}))
#     tod = moby2.scripting.get_tod(params_loadtod)
    
#     # CUT MCE FRAME ERROR                                                                    
#     mce_cuts = moby2.tod.get_mce_cuts(tod)
#     moby2.tod.fill_cuts(tod, mce_cuts, no_noise = no_noise)
    
#     # CUT SOURCES
#     if params_fix.has_key('sigurds_cuts'):
#         f = h5py.File(params_fix.get('sigurds_cuts'), 'r+')
#         grp = f[tod.info.basename]
#         grp.attrs['_moby2_class_name'] ='tod_flags'
#         grp.attrs['_moby2_class_version'] = 1
#         flags_sources = moby2.tod.TODFlags.from_hdf(grp)
#         pos_cuts_sources = flags_sources.get_cuts('cut')
#     else:
#         tod.fplane = products.get_focal_plane(cpar['pointing'], tod.info)
#         mask_params = cpar.get_deep(('source_cuts','mask_params'),{})
#         shift_params = cpar.get_deep(('source_cuts', 'mask_shift_generator'))
#         if shift_params is not None:
#             mask_params['offset'] = products.get_pointing_offset(
#                 shift_params, tod=tod, source_offset=True)
#         if mask_params['offset'] is None:
#             mask_params['offset'] = (0.,0.)
#         if max(mask_params['offset']) > 20./60:
#             mask_params['map_size'] = max(mask_params['offset']) + 10./60
	
#         matched_sources = moby2.ephem.get_sources_in_tod(
#             tod=tod, source_list=cpar.get_deep(('source_cuts','source_list')), pointing_shift=mask_params['offset'])
	
#         pos_cuts_sources = moby2.TODCuts.for_tod(tod, assign=False)
#         print "Found %i sources for this TOD" %len(matched_sources)
#         for source in matched_sources:
#             source_cut = moby2.tod.get_source_cuts(
#                 tod, source[1], source[2], **mask_params)
#             pos_cuts_sources.merge_tod_cuts(source_cut)
#     depot.write_object(pos_cuts_sources,
#                        tag=new_tag_source,
#                        force=True, tod=tod, make_dirs=True)
#     moby2.tod.fill_cuts(tod, pos_cuts_sources, no_noise = no_noise)
    
#     # CUT PLANETS
#     if os.path.exists(
#             depot.get_full_path(
#                 moby2.tod.TODCuts, tag=params['tag_planet'], tod=tod)):
#         pos_cuts_planets = depot.read_object(moby2.TODCuts,
#                                              tag=params['tag_planet'], tod=tod)
#         moby2.tod.fill_cuts(tod, pos_cuts_planets, no_noise = no_noise)
    
#     # GLITCH CUTS
#     cuts_partial = moby2.tod.get_glitch_cuts(tod=tod, params=glitchp)
#     cuts_partial.merge_tod_cuts(mce_cuts)        
#     depot.write_object(cuts_partial, tag=new_tag_partial, tod=tod,
#                        make_dirs = True,
#                        force=True)
#     return pos_cuts_sources, cuts_partial

def replace_pcuts(obs, params, params_fix):
    old_cuts_tag = params['tag_out']
    old_partial_tag = params['tag_partial']
    new_cuts_tag = params_fix['tag_out']
    new_partial_tag = params_fix['tag_partial']

    depot_path = params.get( 'depot', '/data/actpol/depot' )
    depot = moby2.util.Depot(depot_path)

    cuts = moby2.scripting.get_cuts({'depot':depot_path,'tag':old_cuts_tag}, tod=obs)
    pcuts = moby2.scripting.get_cuts({'depot':depot_path,'tag':old_partial_tag}, tod=obs)
    pcuts_new = moby2.scripting.get_cuts({'depot':depot_path,'tag':new_partial_tag}, tod=obs)
    cuts_new = cuts.copy()
    ld = cuts.get_uncut()
    nsamps = np.min([cuts.nsamps, pcuts.nsamps, pcuts_new.nsamps])
    for d in ld:
        cuts_mask = cuts.cuts[d].get_mask()
        pcuts_mask = pcuts.cuts[d].get_mask()
        pcuts_new_mask = pcuts_new.cuts[d].get_mask()
        cuts_new.cuts[d] = moby2.tod.cuts.CutsVector.from_mask(
            np.logical_or( cuts_mask[:nsamps] != pcuts_mask[:nsamps], pcuts_new_mask[:nsamps] ) )
    depot.write_object(cuts_new, tag=new_cuts_tag, tod=obs, make_dirs=True, force=True)

def fix_mr3(obs, params, params_fix):
    print "Fix mr3 for %s"%obs
    depot = moby2.util.Depot(params.get('depot'))
    old_cuts = os.path.exists(
        depot.get_full_path(
            moby2.tod.TODCuts, tag=params['tag_out'], tod=obs))
    new_cuts = os.path.exists(
        depot.get_full_path(
            moby2.tod.TODCuts, tag=params_fix['tag_out'], tod=obs))
    print "Old cuts = %r, New cuts = %r" %(old_cuts, new_cuts)
    if not new_cuts and old_cuts:
        try:
        # new_partial_cuts(obs, params, params_fix)
            replace_pcuts(obs, params, params_fix)
        except:
            print 'Failed for %s' %obs
    else:
        print 'Already exists for %s' %obs

joblib.Parallel(n_jobs=numcores)(
    joblib.delayed(fix_mr3)(obs, params, params_fix) for obs in tl)
