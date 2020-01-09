"""This script is adapted version of the TODList_for_maps.py script in
moby2. It takes the observation details and try to find the list of
tods that correspond to each of the field of observations"""

import fitsio
import moby2, sys, numpy as np, ephem, os
from moby2.util.database import TODList
from moby2.instruments import actpol
from moby2.util import ctime as ct
from moby2.scripting.products import get_filebase
import pandas as pd

from .pathologies_tools import pathoList, get_pwv


def init(config):
    global obs_details, cuts_db, include_time, outdir, selParams, obs_catalog
    obs_details = config.get('obs_details').split()
    cuts_db = config.get('cuts_db', None)
    include_time = config.get('include_time', None)
    outdir = config.get('outdir', ".")
    selParams = config.get('selParams', {
        "liveDets": {"gt": 150},
        "PWV": {"lt": 3},
    })
    obs_catalog = config.get('obs_catalog', None)

def run(p):
    global obs_details, cuts_db, include_time, outdir, selParams, obs_catalog

    if cuts_db is not None:
        pl = pathoList( cuts_db )
    else:
        cpar = moby2.util.MobyDict.from_file(p.i.cutparam)
        pl = pathoList( str(p.i.db) )
    #pl.addPWV2results()
    pl.removeDuplicates()
    pl.addEphem()
    Ndata = pl.ndata

    keys = ['todName', 'liveDets', 'hour', 'hourAfterSunset', 'hourAfterSunrise']
    PL = pd.DataFrame.from_dict( {k:pl.data[k] for k in keys} )

    filename = obs_catalog
    npcat = fitsio.read(filename)
    npcat = npcat.byteswap().newbyteorder()
    catalog = pd.DataFrame.from_records(npcat)
    catalog.index = pd.to_datetime(catalog.date)

    sel = np.logical_and( catalog.obs_type != 'stare', catalog.season == p.i.season)
    sel = np.logical_and( sel, catalog.array == p.i.ar)
    output = pd.merge(catalog[sel], PL, left_on='tod_name', right_on='todName', how='left')
    output.index = pd.to_datetime(output.ctime, unit='s')
    output.sort_index(inplace=True)
    output['PWV'] = output.pwv
    output.PWV[~np.isfinite(output.PWV)] = 0
    output['flag'] = np.zeros(len(output), dtype=int)
    output.flag[~np.isnan(output.liveDets)] += 1

    print "%i TODs" %(len(output))
    print "%i were processed" %((output.flag == 1).sum())
    sel1 = np.ones(len(output), dtype=bool)

    for k in selParams.keys():
        if "lt" in selParams[k].keys():
            sel1 *= (output[k] < selParams[k]['lt']) & ~np.isnan(output[k])
        if "gt" in selParams[k].keys():
            sel1 *= (output[k] > selParams[k]['gt']) & ~np.isnan(output[k])
        print "%i passed the %s criteria (and anteriors)" %(sel1.sum(),k)

    include_time = np.loadtxt(include_time, dtype=int )
    sel2 = np.zeros(len(output), dtype=bool)
    for start, end in include_time:
        sel2 += np.logical_and(output.ctime > start, output.ctime < end)
    print "%i were taken inside the observation times" %(sel2.sum())

    output.flag[sel1*sel2] += 1
    print "%i are good for mapping" %((output.flag==2).sum())

    output['mean_az'] = output.az + output.az_throw / 2.
    output['hour2'] = output.index.hour + output.index.minute/60.

    outDir = outdir + "/" + p.tag
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
