from database import Database as DB
import sys, os
import moby2
import numpy as np
import pandas as pd

# params = moby2.util.MobyDict.from_file(sys.argv[1])


# year = params.get('year')
# max_diff = params.get('max_diff', 3600.)
# depot_path = params.get('depot')
# structure = params.get('structure')
# tag_in = params.get('tag_in')
# tag_out = params.get('tag_out')

year = 2016
array = 'AR4'
max_diff = 3600.
depot_path = '/data/actpol/depot'
tag_cal_out = 'pa4_s16_BScal'
tag_height_out = 'pa4_s16_BSheight'


depot = moby2.util.Depot(depot_path)
db = DB(year)

db.load_acqs()
db.load_tods()

# arrays = np.unique(db.tods.array)
# for ar in arrays:
sel_acqs = db.acqs.array == array
sel_tods = db.tods.array == array


merged_table = pd.merge(db.acqs[sel_acqs], db.tods[sel_tods],
                        how='outer',
                        left_on='id', right_on='mce_data_acq_id')

ctime_bs = merged_table.ctime[merged_table.suffix=='bc1_step']
ctime_bs.sort()
ctime_tods = merged_table.ctime[merged_table.suffix=='dat']
idx = np.digitize(ctime_tods, ctime_bs) - 1    
sel = ~ merged_table.name.isnull()


for ct_tod, ct_bs in zip(ctime_tods,ctime_bs.values[idx]):
    name = merged_table.name[merged_table.ctime == ct_tod].values[0]
    t_diff = ct_tod - ct_bs
    if name is not None and name is not np.NaN and t_diff < max_diff:
        try:
            cal_in = moby2.util.MobyDict.from_file(
                os.path.join('/data/actpol/season4/bias_step_calibration_v2/mce4',
                             str(ct_bs)[:5],
                             str(ct_bs)+'.cal'))
            cal_out = moby2.Calibration(det_uid = cal_in['det_uid'])
            cal_out.set_property(['cal'],[cal_in['cal']])
            depot.write_object(cal_out, tag=tag_cal_out, tod=name)
        except:pass
        try:
            height_in = moby2.util.MobyDict.from_file(
                os.path.join('/data/actpol/season4/bias_step_calibration/mce4',
                             str(ct_bs)[:5],
                             str(ct_bs)+'.H'))
            height_out = moby2.Calibration(det_uid = height_in['det_uid'])
            height_out.set_property(['cal'],[height_in['cal']])
            depot.write_object(height_out, tag=tag_height_out, tod=name)
        except: pass
    
        
