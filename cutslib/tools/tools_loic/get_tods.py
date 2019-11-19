import sys, os
import moby2
from moby2.util.database import TODList
moby2.pointing.set_bulletin_A()
# moby2.pointing.set_bulletin_A()



cparFile = "/home/lmaurin/cuts/2015/AR1/MR1_PA1_2015/cutParams_MR1_PA1_2015_nohwp.par"

source_scans = '/home/lmaurin/TODLists/2015_AR1_season.txt'#params['source_scans']
obsnames = TODList.from_file(source_scans)
for obs in obsnames:
    out = os.popen('python get_sources_in_tod.py %s %s' %(obs, cparFile))
