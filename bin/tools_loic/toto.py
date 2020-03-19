"""Get the number of sources in a given list of TODs. Likely obsolete now"""

import moby2
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, Longitude, Latitude, Angle
import astropy

moby2.pointing.set_bulletin_A()
obs = '1439386380.1439582265.ar2'
from moby2.util.database import TODList
todlist = TODList.from_file(
    '/home/lmaurin/TODLists/2015_ar3_season_nohwp.txt')

for obs in todlist:
    print(obs)
    try:
        tod = moby2.scripting.get_tod({'filename':obs, 'read_data':False,
                                       'repair_pointing':True})
        sources = moby2.ephem.get_sources_in_tod(
            tod,
            source_list='/home/lmaurin/actpol_data_shared/BrightSources/mr3c_merged_gt2mK.txt')
        f = open('Nsources_gt2mK_in_tod.txt', 'a')
        f.write('{} {}\n'.format(obs, len(sources)))
        f.close()
    except: pass
