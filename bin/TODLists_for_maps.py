import fitsio
import moby2, sys, numpy, ephem, os
from moby2.scripting.pathologies_tools import pathoList, get_pwv
from moby2.util.database import TODList
from moby2.instruments import actpol
from moby2.util import ctime as ct
from moby2.scripting.products import get_filebase
import database as DB
import pandas as pd
np = numpy


params = moby2.util.MobyDict.from_file(sys.argv[1])
cuts_db = params.get('cuts_db', None)
obs_details = params.get("obs_details", ["deep6"])
array = params.get("array", "ar1")


if cuts_db is not None:
    pl = pathoList( cuts_db )
else:
    cpar = moby2.util.MobyDict.from_file(params.get('cutparams'))
    pl = pathoList( os.path.join(cpar.get("outdir"),cpar.get("report")+".db") )
#pl.addPWV2results()
pl.removeDuplicates()
pl.addEphem()
Ndata = pl.ndata

keys = ['todName', 'liveDets', 'hour', 'hourAfterSunset', 'hourAfterSunrise']
PL = pd.DataFrame.from_dict( {k:pl.data[k] for k in keys} )

# db = DB.Database(params.get('season'))
# db.create_query('tods')
# db.add_condition( ('array', '=', params.get('array')) )
# db.make_query()

# sel = db.data.obs_type != 'stare'
# output = pd.merge(db.data[sel], PL, left_on='name', right_on='todName', how='left')
# output.index = pd.to_datetime(output.ctime_start, unit='s')
# output.sort_index(inplace=True)
# output['PWV'] = output.pwv
# output.PWV[~np.isfinite(output.PWV)] = 0
# output['flag'] = np.zeros(len(output), dtype=int)
# output.flag[~np.isnan(output.liveDets)] += 1

filename = '/home/lmaurin/actpol_data_shared/ObsCatalog/catalog.fits'
npcat = fitsio.read(filename)
npcat = npcat.byteswap().newbyteorder()
catalog = pd.DataFrame.from_records(npcat)
catalog.index = pd.to_datetime(catalog.date)

sel = np.logical_and( catalog.obs_type != 'stare', catalog.season == str(params.get('season')) )
sel = np.logical_and( sel, catalog.array == params.get('array'))
output = pd.merge(catalog[sel], PL, left_on='tod_name', right_on='todName', how='left')
output.index = pd.to_datetime(output.ctime, unit='s')
output.sort_index(inplace=True)
output['PWV'] = output.pwv
output.PWV[~np.isfinite(output.PWV)] = 0
output['flag'] = np.zeros(len(output), dtype=int)
output.flag[~np.isnan(output.liveDets)] += 1

print("%i TODs" %(len(output)))
print("%i were processed" %((output.flag == 1).sum()))
sel1 = np.ones(len(output), dtype=bool)
selParams = params.get('selParams', {})
for k in list(selParams.keys()):
    if "lt" in list(selParams[k].keys()):
        sel1 *= (output[k] < selParams[k]['lt']) & ~np.isnan(output[k])
    if "gt" in list(selParams[k].keys()):
        sel1 *= (output[k] > selParams[k]['gt']) & ~np.isnan(output[k])
    print("%i passed the %s criteria (and anteriors)" %(sel1.sum(),k))

include_time = np.loadtxt(params.get('include_time'), dtype=int )
sel2 = np.zeros(len(output), dtype=bool)
for start, end in include_time:
    sel2 += np.logical_and(output.ctime > start, output.ctime < end)
print("%i were taken inside the observation times" %(sel2.sum()))

output.flag[sel1*sel2] += 1
print("%i are good for mapping" %((output.flag==2).sum()))

output['mean_az'] = output.az + output.az_throw / 2.
output['hour2'] = output.index.hour + output.index.minute/60.

outDir = params.get("outDir", ".")
if not(os.path.isdir(outDir)): os.makedirs(outDir)


for field in np.unique(output.obs_detail):
    suboutput = output[output.obs_detail==field]
    filename = 'selectedTODs_%s.txt' %field
    suboutput.to_csv( os.path.join(outDir,filename),
                      header=['# TOD', 'hour', 'altitude', 'azimuth',
                              'PWV', 'cut_status', 'field'],
                      columns=['tod_name', 'hour2', 'alt', 'mean_az',
                               'pwv', 'flag', 'obs_detail'],
                      index=False,
                      na_rep='-', float_format='%.2f', sep='\t')

filename = 'selectedTODs_allfields.txt'
output.to_csv( os.path.join(outDir,filename),
                  header=['# TOD', 'hour', 'altitude', 'azimuth',
                          'PWV', 'cut_status', 'field'],
                  columns=['tod_name', 'hour2', 'alt', 'mean_az',
                           'pwv', 'flag', 'obs_detail'],
                  index=False,
                  na_rep='-', float_format='%.2f', sep='\t')
