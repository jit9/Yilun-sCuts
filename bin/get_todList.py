"""This script generates a list of TODs for preliminary cuts study.
It takes parameters from an external parameter file. See `todlist.par`
as an example.

Example configuration script:

######################################################################
# output format and season array specification
output = '{season}_{array}_uranus_nohwp.txt'
season = '2017'
array = 'ar4'

# particular observations of interests
obs_detail = ['uranus']
# obs_type = "scan"

# whether to include hwp
hwp_epoch = False

# whether to select a randomized subset
do_random = False
nsamples = 1000
######################################################################

This parameter file specifies a set of criteria used to narrow
down to a subset of tods.

Usage:
>>> python get_todList.py tod.par

"""

import moby2
import sys
import fitsio
import pandas as pd
import numpy as np
fb = moby2.scripting.get_filebase()

# load parameter file
params = moby2.util.MobyDict.from_file( sys.argv[1] )

# load catalog of tods and convert to pandas dataframe
filename = params.get('ObsCatalog','/home/lmaurin/actpol_data_shared/ObsCatalog/catalog.fits')
npcat = fitsio.read(filename)
npcat = npcat.byteswap().newbyteorder()
catalog = pd.DataFrame.from_records(npcat)
catalog.index = pd.to_datetime(catalog.date)

# select tods based on the specified criteria
sel = np.ones(catalog.shape[0], dtype=bool)
for k in params.keys():
    if k in catalog.columns:
        print '%s: ' %k,
        p = params[k]
        if type(p) == list:
            sel *= catalog[k].isin(p)
        elif type(p) == bool:
            if p:
                sel *= catalog[k] != 'none'
            else:
                sel *= catalog[k] == 'none'
        else:
            sel *= catalog[k] == p
        print 'left with %i TODs' %sel.sum()
    else:
        print '%s is not a valid parameter, skip it' %k

# prepare output
output_filename = params['output'].format(**params)
print 'Output %i tods to %s' %(sel.sum(),output_filename)

# whether to randomize the list
if not params["do_random"]:
    catalog.tod_name[sel].to_csv(output_filename, index=False)
else:
    catalog.tod_name[sel].sample(n=params["nsamples"]).to_csv(output_filename, index=False)

######################################################################
# end of get_todList.py

