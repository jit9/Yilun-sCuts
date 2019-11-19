import moby2
import cPickle
import numpy as np
import sys


params = moby2.util.MobyDict.from_file(sys.argv[1])
cPar = moby2.util.MobyDict.from_file(params['cutParams'])


detector_lists = cPar['pathologyParams']['detectorLists']

input_live = moby2.util.MobyDict.from_file(detector_lists['live'])
input_live['det_uid'] = np.asarray(input_live['rows'])*32 + np.asarray(input_live['cols'])
input_exclude = moby2.util.MobyDict.from_file(detector_lists['exclude'])
input_exclude['det_uid'] = np.asarray(input_exclude['rows'])*32 + np.asarray(input_exclude['cols'])
input_dark = moby2.util.MobyDict.from_file(detector_lists['dark'])


season = cPickle.load(open(params['critfile'], 'r'))
liveSum = (season['sel']).sum(axis=1)
exclude = np.arange(season['sel'].shape[0])[liveSum==0]

set_live = set(input_live['det_uid'])
live_in = len(set_live)
set_exclude = set(exclude)
set_live.difference_update(set_exclude)
live_out = len(set_live)
print "%i detectors have been excluded" %(live_in-live_out)

def generate_mobydict(input_set):
    uid = np.array(list(input_set))
    col = uid%32
    row = uid/32
    output = moby2.util.MobyDict()
    output['det_uid'] = uid.tolist()
    output['rows'] = row.tolist()
    output['cols'] = col.tolist()
    return output

output_live = generate_mobydict(set_live)
output_exclude = generate_mobydict(set_exclude)

output_live.write_to_file('live_%s.dict' %params['tag_out'])
output_exclude.write_to_file('exclude_%s.dict' %params['tag_out'])

