import moby2, sys, os
import numpy as np
from moby2.scripting.pathologies_tools import pathoList
import database as DB
import pandas as pd



params = moby2.util.MobyDict.from_file(sys.argv[1])
array = params.get("array")
season = params.get("season")
cuts_db = params.get('cuts_db')


# Load the cuts database into a DataFrame
pl = pathoList( cuts_db )


pl.addEphem()
keys = ['todName', 'liveDets',
        'hour', 'hourAfterSunset', 'hourAfterSunrise',
        'length', 'liveFrac']
PL = pd.DataFrame.from_dict( {k:pl.data[k] for k in keys} )
del pl

# Load the TOD infos database
db = DB.Database(params.get('season'))
db.create_query('tods')
db.add_condition( ('array', '=', params.get('array')) )
db.make_query()

# Merge them, add time index
output = pd.merge(db.data, PL, left_on='name', right_on='todName', how='left')
output.index = pd.to_datetime(output.ctime_start, unit='s')
output.sort_index(inplace=True)
output['PWV'] = output.pwv
output.length[np.abs(output.length) > 1e3] = np.NaN

# Create flag for mapping
output['flag'] = np.zeros(len(output), dtype=int)
output.flag[~np.isnan(output.liveDets)] += 1

sel1 = np.ones(len(output), dtype=bool)
selParams = params.get('selParams', {})
for k in selParams.keys():
    if "lt" in selParams[k].keys():
        sel1 *= (output[k] < selParams[k]['lt']) & ~np.isnan(output[k])
    if "gt" in selParams[k].keys():
        sel1 *= (output[k] > selParams[k]['gt']) & ~np.isnan(output[k])
    print sel1.sum()

include_time = np.loadtxt(params.get('include_time'), dtype=int )
sel2 = np.zeros(len(output), dtype=bool)
for start, end in include_time:
    sel2 += (output.ctime_start > start) & (output.ctime_start < end)
print sel2.sum()

output.flag[sel1*sel2] += 1




# Some intersting statistics
output.flag.groupby(output.obs_detail).value_counts()
output.liveDets.groupby(output.obs_detail).median()

# output.flag.groupby(output.hwp_epoch).value_counts()
# output.liveDets.groupby(output.hwp_epoch).median()

output.length.groupby(output.obs_detail).sum() / 60.
# output.length.groupby(output.hwp_epoch).sum() / 60.

output['liveDetsFrac'] = output.liveDets / output.liveDets.max()
output.liveDetsFrac.groupby(output.obs_detail).median()
# output.liveDetsFrac.groupby(output.hwp_epoch).median()

output.liveFrac.groupby(output.obs_detail).median()
output['effectiveTime'] = output.length * output.liveDetsFrac * (1-output.liveFrac)
output.effectiveTime.groupby(output.obs_detail).sum() / 60.
# output.effectiveTime.groupby(output.hwp_epoch).sum() / 60.
