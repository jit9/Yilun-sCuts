"""This seems to be a useful script that returns the list of TODs
that scan across a given source"""

import numpy as np
from database import Database as DB
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, Longitude, Latitude, Angle
import sys
import moby2


params = moby2.util.MobyDict.from_file(sys.argv[1])

year = params.get('year')
array = params.get('array')
array = array.upper()
source = params.get('source')

db = DB(year)
db.load_tods()
db.tods = db.tods[db.tods.array==array]
db.tods = db.tods[db.tods.mean_alt<90.]
#db.tods = db.tods[np.logical_and(db.tods.scan_speed>0.,np.isfinite(db.tods.scan_speed))]
times = Time((db.tods.ctime_start+db.tods.ctime_end)/2., format='unix')
min_az = Longitude(db.tods.min_az, unit='deg')
max_az = Longitude(db.tods.max_az, unit='deg')
throw = Angle(db.tods.max_az - db.tods.min_az, unit='deg')
mean_alt = Latitude(db.tods.mean_alt, unit='deg')


act_site = EarthLocation(lon=-1.183116812024908 * u.rad,
                         lat=-0.40070141631911815 * u.rad,
                         height=5188 * u.m)

source_coords = SkyCoord(ra=source[0], dec=source[1], unit='deg',
                         location=act_site,
                         obstime=times)
altaz = source_coords.transform_to('altaz')

min_coord = SkyCoord(min_az, mean_alt, frame='altaz', location=act_site)
max_coord = SkyCoord(max_az, mean_alt, frame='altaz', location=act_site)

# d_min = altaz.separation(min_coord)
# d_max = altaz.separation(max_coord)
#sel = (d_min<throw+1*u.deg) * (d_max<throw+1*u.deg) * ( np.abs(altaz.alt-mean_alt) < 1*u.deg)

sel_alt = np.abs( mean_alt - altaz.alt ) < 1.5*u.deg
sel_az = np.logical_and(
    altaz.az - min_coord.az > -1.5*u.deg,
    max_coord.az - altaz.az > -1.5*u.deg)
sel = np.logical_and(sel_alt,sel_az)

tods = db.tods.name[sel]
tods.to_csv(
    'tods_ra%.2f_dec%.2f_%s_%s.txt' %(source[0],source[1],array,year),
    index=False)
