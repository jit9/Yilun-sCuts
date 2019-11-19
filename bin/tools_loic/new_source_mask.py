import numpy as np
import moby2
from moby2.scripting import products
moby2.pointing.set_bulletin_A()           

obs = '1475836450.1475851025.ar2' # Tau A
#obs = '1477209039.1477222347.ar2'
cuts = moby2.scripting.get_cuts({'depot':'/data/actpol/depot','tag':'mr3c_pa2_s16'}, tod=obs)
ld = cuts.get_uncut()
pointing = {
    'source': 'fp_file',
    'filename': "/data/actpol/actpol_data_shared/RelativeOffsets/template_ar2_s16_170131.txt"}                        
pointing_shift =  {
    'source':'file',
    'filename':'/data/actpol/actpol_data_shared/TODOffsets/tod_offsets_2016_170131.txt',
    'columns': [0,3,4],
    'rescale_degrees': 1./60}
offset = products.get_pointing_offset(
        pointing_shift, tod=obs, source_offset=True)
offset_ra = np.deg2rad(offset[1])
offset_dec = np.deg2rad(offset[0])


tod = moby2.scripting.get_tod({'filename':obs, 'read_data':True, 'fix_sign':True,'repair_pointing':True})
moby2.tod.detrend_tod(tod)
moby2.tod.remove_median(tod)
fc_l = 1.
df_l = 0.5
filt = moby2.tod.filters.sine2highPass( tod=tod, fc=fc_l, df=df_l)
moby2.tod.filter.apply_simple(tod.data,filt, detrend=True,dets=ld)

tod.fplane = products.get_focal_plane(pointing, tod.info)

f = open('/home/lmaurin/actpol_data_shared/BrightSources/sources.txt', 'r')
source_list = f.readlines()
source_list = [(s.strip('\n'), 'source') for s in source_list]
f.close()
matched_sources = moby2.ephem.get_sources_in_patch(
    tod=tod, source_list=source_list)


tod.alt -= np.deg2rad(offset[0])
tod.az += np.deg2rad(offset[1])
wand = moby2.pointing.ArrayWand.for_tod(tod, coords='ra_dec')
ra, dec = wand.get_coords(tod.fplane)
ra = np.rad2deg(ra)
dec = np.rad2deg(dec)

map_pix = 0.5 / 60
maskMap = moby2.mapping.fits_map.spaceMap.simpleMap(
    (ra.min(), ra.max()), (dec.min(), dec.max()), (map_pix,map_pix), dtype='float32', wtype='int32')
maskMap.x = maskMap.x
maskMap.y = maskMap.y
radius = 5./60
for source in matched_sources:
    r, d = np.rad2deg(source[1]), np.rad2deg(source[2])
    my_mask = ((maskMap.x-r)**2 + (maskMap.y-d)**2 < radius**2)
    maskMap.data[my_mask] = 1.
ny, nx = maskMap.data.shape


gridding = moby2.pointing.GridPixelization.forFitsMap(maskMap)
gridding.dx0 *= -1
proj = moby2.pointing.WandProjector(wand, gridding, tod.fplane)
# Warning this creates TOD-sized data...                                                                                                 
mask_vect = np.asarray( proj.deproject_from_map(maskMap.data), dtype=bool )




