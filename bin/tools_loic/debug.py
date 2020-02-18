import os
import joblib
import multiprocessing
import numpy as np
import moby2
from moby2.scripting import products
import sys

args = sys.argv

obs = '1443426331.1443459471.ar2'

params = moby2.util.MobyDict.from_file(args[1])
params_fix = moby2.util.MobyDict.from_file(args[2])
numcores = params_fix['ncores']
from moby2.util.database import TODList
tl = TODList.from_file(params['source_scans'])

new_tag_source = params_fix["tag_source"]
new_tag_partial = params_fix["tag_partial"]

IV_wait = params.get("IV_wait",0)
rebias_wait = params.get("rebias_wait",0)
end = params.get("sample_end",None)


depot = moby2.util.Depot(params.get('depot'))
cutParams = params.get('cutParams')
cpar = moby2.util.MobyDict.from_file(cutParams)
cpar['source_cuts'] = params_fix['source_cuts']
pathop = cpar['pathologyParams']
glitchp = cpar['glitchParams']

# Fill cuts parameter                                            \
    
no_noise = not(cpar.get("fillWithNoise",True))

user_config = moby2.util.get_user_config()
moby2.pointing.set_bulletin_A(params=user_config.get('bulletin_A_\
settings'))

# Find IV/rebias gap time                                        \
    
ct = int(obs.split("/")[-1].split(".")[0])
ctimes = (ct-IV_wait,ct)
from moby2.instruments import actpol as inst
try:
    config_file = params.get('manifest_conf', None)
    db = inst.TODDatabase(config_file=config_file)
except:
    db = None
if (db is None) or (db.db is None):
    offset_IV = 0
    offset_rebias = 0
else:
    recs = db.select_acqs(suffix='iv', ctime=ctimes)
    if len(recs) > 0:
        offset_IV = (IV_wait - ct + recs[-1].ctime)*400
    else:
        offset_IV = 0
    ctimes = (ct-rebias_wait,ct)
    recs = db.select_acqs(suffix='bc1_step', ctime=ctimes)
    if len(recs) > 0:
        offset_rebias = (rebias_wait - ct + recs[-1].ctime)*400
    else:
        offset_rebias = 0
offset = max( offset_IV, offset_rebias, params.get("offset",0)*400 )
 

 
params_loadtod = {
    'filename':obs,
    'fix_sign':True,
    'start':offset,
    'end':end,
    'repair_pointing':True
}
params_loadtod.update(params.get("params_loadtod",{}))
tod = moby2.scripting.get_tod(params_loadtod)

# CUT MCE FRAME ERROR                                            \
    
mce_cuts = moby2.tod.get_mce_cuts(tod)
moby2.tod.fill_cuts(tod, mce_cuts, no_noise = no_noise)

# CUT SOURCES                                                     
tod.fplane = products.get_focal_plane(cpar['pointing'], tod.info)
#    pointing_shift = (0,0)                                           
mask_params = cpar.get_deep(('source_cuts','mask_params'),{})
shift_params = cpar.get_deep(('source_cuts', 'mask_shift_generato\
r'))
if shift_params is not None:
    mask_params['offset'] = products.get_pointing_offset(
        shift_params, tod=tod, source_offset=True)
    if mask_params['offset'] is None:
        mask_params['offset'] = (0.,0.)
    if max(mask_params['offset']) > 20./60:
        mask_params['map_size'] = max(mask_params['offset']) + 10./60

matched_sources = [('J2000', 1.4596188534428578, 0.38414696836395196, 'source')]

pos_cuts_sources = moby2.TODCuts.for_tod(tod, assign=False)
print(tod.data.shape)
print("Found %i sources for this TOD" %len(matched_sources))
for source in matched_sources:
    print(source)
    print(mask_params)
    source_cut = moby2.tod.get_source_cuts(
        tod, source[1], source[2], **mask_params)
    pos_cuts_sources.merge_tod_cuts(source_cut)

wand = moby2.pointing.ArrayWand.for_tod_source_coords(tod, ref_coord=(source[1],source[2]),scan_coords=True)




