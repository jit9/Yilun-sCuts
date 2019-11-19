import moby2
import pandas as pd
import zipfile



obslist= ['1409892592.1409892632.ar2', # 758 live dets
          '1414970015.1414970622.ar2', # 384 live dets
          '1411609539.1411609591.ar2'] # 19 livr dets

for obs in obslist:
    tod = moby2.scripting.get_tod({'filename':obs, 'read_data':True, 'fix_sign':True,'repair_pointing':True})
    moby2.tod.remove_median(tod)
    cal = moby2.scripting.get_calibration({'type':'iv','source':'data'}, tod=tod)
    moby2.tod.apply_calibration(tod.data, cal.det_uid, cal.cal)
    
    data = pd.DataFrame(data=tod.data.T, index=tod.ctime, columns=['tesdata%04i' %i for i in xrange(tod.data.shape[0])])
    
    
    
    z = zipfile.ZipFile(tod.info.filename, 'r')
    files = z.namelist()
    z.close()
    channels = [f.strip('.slm') for f in files if not 'tesdata' in f and not 'Extras' in f]
    for channel_name in channels:
        try:
            channel = tod.get_hk(channel_name,fix_gaps=True)
            assert channel.size == data.shape[0]
            data[channel_name] = channel
        except:
            pass
            
            
    data.to_csv(tod.info.basename+'.csv', index_label='ctime')




