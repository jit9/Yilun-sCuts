"""This script is adapted version of the TODList_for_maps.py script in
moby2. It takes the observation details and try to find the list of
tods that correspond to each of the field of observations

Parameters:
------------
cuts_db (str): patho db file, if not interested in the default
source_scan (bool): whether we want to restrict to source scan list
selParams (str): dictionary of crits to select TODs from catalog

"""

import fitsio
import moby2, sys, numpy as np, ephem, os
from moby2.util.database import TODList
from moby2.instruments import actpol
from moby2.util import ctime as ct
from moby2.scripting.products import get_filebase
import pandas as pd
from cutslib.pathologies_tools import pathoList, get_pwv
from cutslib import Catalog


class Module:
    def __init__(self, config):
        self.cuts_db = config.get('cuts_db', None)
        self.include_time = config.get('include_time', None)
        self.outdir = config.get('outdir', ".")
        self.source_scan = config.getboolean('source_scan', False)
        self.selParams = config.get('selParams', {
            "liveDets": {"gt": 150},
            "PWV": {"lt": 3},
        })
        # load dict from string if selParams is given by str
        if type(self.selParams) == str:
            import json
            self.selParams = json.loads(self.selParams)

    def run(self, p):
        cuts_db = self.cuts_db
        include_time = self.include_time
        outdir = self.outdir
        selParams = self.selParams
        source_scan = self.source_scan
        cpar = moby2.util.MobyDict.from_file(p.i.cutparam)
        if cuts_db is not None:
            pl = pathoList( cuts_db )
        else:
            pl = pathoList( str(p.i.db) )
        # pl.addPWV2results()
        pl.removeDuplicates()
        pl.addEphem()
        Ndata = pl.ndata

        keys = ['todName', 'liveDets', 'hour', 'hourAfterSunset', 'hourAfterSunrise']
        PL = pd.DataFrame.from_dict( {k:pl.data[k] for k in keys} )
        # load catalog
        cat = Catalog()
        catalog = cat.data
        # base sel for season, array, etc
        sel = np.logical_and(catalog.obs_type != 'stare',
                             catalog.season == p.i.season)
        sel = np.logical_and( sel, catalog.array == p.i.ar)
        # if we want to restrict to source_scan only
        if source_scan:
            from moby2.util.database import TODList
            cpar = moby2.util.MobyDict.from_file(p.i.cutparam)
            source_list = TODList.from_file(cpar.get('source_scans'))
            sel = np.logical_and(sel, catalog.tod_name.isin(source_list))

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

        for k in list(selParams.keys()):
            if "lt" in list(selParams[k].keys()):
                sel1 *= (output[k] < selParams[k]['lt']) & ~np.isnan(output[k])
            if "gt" in list(selParams[k].keys()):
                sel1 *= (output[k] > selParams[k]['gt']) & ~np.isnan(output[k])
            print("%i passed the %s criteria (and anteriors)" %(sel1.sum(),k))

        include_time = np.loadtxt(include_time, dtype=int )
        sel2 = np.zeros(len(output), dtype=bool)
        for start, end in include_time:
            sel2 += np.logical_and(output.ctime > start, output.ctime < end)
        print("%i were taken inside the observation times" %(sel2.sum()))

        output.flag[sel1*sel2] += 1
        print("%i are good for mapping" %((output.flag==2).sum()))

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
