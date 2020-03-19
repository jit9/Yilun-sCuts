"""Generate the list of sources associated with a given TOD."""

import moby2
import sys
moby2.pointing.set_bulletin_A()


obs, cparFile = sys.argv[1:]
cpar = moby2.util.MobyDict.from_file(cparFile)

tod = moby2.scripting.get_tod({'filename':obs,
                               'read_data':False})
tod.fplane = moby2.scripting.products.get_focal_plane(
    cpar['pointing'], tod.info)
matched_sources = moby2.ephem.get_sources_in_patch(
    tod=tod, source_list=None)

if len(matched_sources) > 0:
    f = open('2015_PA1_sources.txt','a')
    f.write('%s %s\n' %(obs, matched_sources))
    f.close()
