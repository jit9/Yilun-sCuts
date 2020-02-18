import match_acqs
from database import Database as DB
import os
import shutil

year = 2017
path_input = '/data/actpol/depot/bias_step/2017/calibration_20190802'
path_output = '/data/actpol/depot/Calibration'
tag_output = 'calibration_20190802'

db = DB(year)
db.load_tods()
sel = (db.tods.ctime_start > 1483228800) & (db.tods.ctime_start < 1516900000)
db.tods = db.tods[sel]

db.tods['previous_BS_ctime'] = match_acqs.get_previous_BS(db.tods.ctime_start, db.tods.array)
#db.tods['next_BS_ctime'] = match_acqs.get_next_BS(db.tods.ctime_start, db.tods.array)


tot = 0
f = open('missing_BS.txt', 'w')
f.write('# Ctimes of missing BS\n')
f.close()
for ctime_tod, ctime_bs, array, todname in zip(db.tods.ctime_start, db.tods.previous_BS_ctime, db.tods.array, db.tods.name):
    filename_input = os.path.join(path_input, 'pa%s' %array[-1], str(ctime_bs)[:5], '%s.cal' %ctime_bs)
    filename_output = os.path.join(path_output, 'pa%s_s%s_bs_%s' %(array[-1],str(year)[-2:],tag_output), str(ctime_tod)[:5], '%s.cal' %todname)
    
    if os.path.exists(filename_input) and not os.path.exists(filename_output): 
        tot+=1
        f = open('missing_BS.txt', 'a')
        f.write('%s\n' %ctime_bs)
        f.close()
        if not os.path.exists(os.path.join(path_output, 'pa%s_s%s_bs_%s' %(array[-1],str(year)[-2:],tag_output), str(ctime_tod)[:5])):
            os.makedirs(os.path.join(path_output, 'pa%s_s%s_bs_%s' %(array[-1],str(year)[-2:],tag_output), str(ctime_tod)[:5]))
        shutil.copy(filename_input, filename_output)
    else: print(filename_output, "already exists")
print(tot, db.tods.array.size)

db.tods.to_csv('%s_bs_%s_for_tod.txt' %(year, tag_output), columns=['name', 'previous_BS_ctime'], index=False, header=['#name', 'BS_ctime'], sep='\t') 
