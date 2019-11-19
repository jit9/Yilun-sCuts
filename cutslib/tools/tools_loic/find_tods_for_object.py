import numpy as np
import moby2
from moby2.instruments import actpol
from moby2.util.database import TODList
from moby2.ephem import actEphem
from moby2.scripting.products import get_filebase
db = actpol.TODDatabase(config_file='/data/manifest_conf/manifest_2014.conf')
fb = get_filebase()
from matplotlib import pyplot as plt

ids = db.select_tods()

eph = actEphem.ACTEphem()
Sun = eph._objects['Sun']
Moon = eph._objects['Moon']
    

def angdist( x, y ):
    az1, alt1 = x
    az2, alt2 = y
    return np.arccos( 
        np.cos(alt1)*np.cos(alt2)*np.cos(az1-az2)
        + np.sin(alt1)*np.sin(alt2) )



def find_tods_for_object( coords_object ):
    ids = db.select_tods(obs_detail='deep56')
    tl = TODList()
    ids_sel = []
    plt.ioff()
    plt.figure()
    for id in ids:
        RA = np.asarray([
            id.min_RA_south, id.max_RA_south, id.min_RA_north, id.max_RA_north
            ])
        RA[RA>180] -= 360
        coords_tod = (
            RA[0], RA[1], RA[2], RA[3],
            id.min_dec, id.max_dec)
        inside = inside_tod(coords_object, coords_tod)
        if inside: color = 'r'
        else: color = 'b'
        plt.plot(
            [RA[0], RA[1],
             RA[3],RA[2],
             RA[0]],
            [id.min_dec, id.min_dec,
             id.max_dec, id.max_dec,
             id.min_dec],
            color, alpha=0.2)
        if inside: 
            tl.append(id.basename)
            ids_sel.append(id)
    return tl, ids_sel

def inside_tod(coords_point, coords_tod ):
    RA_min_0, RA_max_0, RA_min_1, RA_max_1, DEC_0, DEC_1= coords_tod
    RA,DEC = coords_point
    DEC_min = min(DEC_0, DEC_1)
    DEC_max = max(DEC_0, DEC_1)
    if DEC > DEC_max or DEC < DEC_min:
        return False
    elif (DEC - DEC_0)*(RA_min_1-RA_min_0) - (RA-RA_min_0)*(DEC_1-DEC_0) > 0:
        return False
    elif (DEC - DEC_0)*(RA_max_1-RA_max_0) - (RA-RA_max_0)*(DEC_1-DEC_0) < 0:
        return False
    else:
        return True


rad2deg = 180 / np.pi
deg2rad = np.pi / 180
# plt.ioff()
# plt.figure()


f = open('distance_sun_moon.txt', 'w')
f.write('# tod   sun   moon\n')
f.close()

for id in ids:
#id = ids[5000]
    try:
        filename = fb.filename_from_name(id.basename, single=True)
        tod = moby2.scripting.products.get_tod(
            {'filename':filename, 'read_data':False, 'repair_pointing':True} )
        
        eph.set_ctime(tod.info.ctime)
        
        Sun.compute(tod.info.ctime)
        Sun_alt, Sun_az = eph.radec_to_altaz(Sun.ra*rad2deg,
                                             Sun.dec*rad2deg)
        
        Moon.compute(tod.info.ctime)
        Moon_alt, Moon_az = eph.radec_to_altaz(Moon.ra*rad2deg,
                                               Moon.dec*rad2deg)
        
        
        # plt.plot(tod.az*rad2deg, tod.alt*rad2deg, 'b.')
        # plt.plot(Sun_az*rad2deg, Sun_alt*rad2deg, 'y*')
        # plt.plot(Moon_az*rad2deg, Moon_alt*rad2deg, 'ko')
        
        dSun = 180 / np.pi * min(
            angdist([tod.az.min(), tod.alt.mean()], [Sun_az, Sun_alt]),
            angdist([tod.az.max(), tod.alt.mean()], [Sun_az, Sun_alt]) )
        dMoon = 180 / np.pi * min(
            angdist([tod.az.min(), tod.alt.mean()], [Moon_az, Moon_alt]),
            angdist([tod.az.max(), tod.alt.mean()], [Moon_az, Moon_alt]) )
        print "Sun distance: %.2f deg, Moon distance: %.2f deg" %(dSun,
                                                                  dMoon)
        dsm = angdist([Moon_az, Moon_alt], [Sun_az, Sun_alt]) * rad2deg
        print "distance sun-moond: %.2f deg" %dsm
        print ""
        f = open('distance_sun_moon.txt', 'a')
        f.write('%s   %.2f   %.2f\n' %(id.basename, dSun, dMoon))
        f.close()
    except:pass
