"""The hdf5 file for source cuts has a multiple write problem. This
is because there are a few missing attributes that need to be added
manually for each TOD. This means that running over each TOD requires
a writing permission. This is not allowed for hdf5 file, so this is
a problem. A solution to this is to first fill-up the required fields
beforehand to make sure the hdf5 file is correct, then only invoke
the read permission during parallel processing. This should solve
the multiple read problem. """

import h5py
import moby2
from moby2.util.database import TODList
from moby2.scripting import products
import sys

# parameters
freq = 'f090'
input_file = "/home/yguan/data/src_cuts_s17v5_{}.hdf".format(freq)

# load pickle file
params = moby2.util.MobyDict.from_file(sys.argv[1])

# setup logger
logger = moby2.util.log.logger
logger.init_from_params(params.get('moby_options', {}))

# load source scans
source_scans = params.get('source_scans')
obsnames = TODList.from_file(source_scans)

# prepare to calculate offset
IV_wait = params.get("IV_wait",0)
rebias_wait = params.get("rebias_wait",0)
end = params.get("sample_end",None)
if end is not None:
    end = -end


f = h5py.File(input_file, "r+")

for obs in list(f.keys()):
    # Find IV/rebias gap time
    ct = int(obs.split("/")[-1].split(".")[0])
    ctimes = (ct-IV_wait,ct)
    if products._get_instrument() == "actpol":
        from moby2.instruments import actpol as inst
    elif products._get_instrument() == "mbac":
        from moby2.instruments import mbac as inst
    try:
        config_file = params.get('manifest_conf', None)
        db = inst.TODDatabase(config_file=config_file)
    except:
        db = None
    if (db is None) or (db.db is None):
        logger.trace(0, "Database not accessible, IV/rebias offsets set to 0")
        offset_IV = 0
        offset_rebias = 0
    else:
        recs = db.select_acqs(suffix='iv', ctime=ctimes)
        if len(recs) > 0:
            offset_IV = (IV_wait - ct + recs[-1].ctime)*400
        else:
            offset_IV = 0
        logger.trace(0,"IV offset set to %d"%offset_IV)
        ctimes = (ct-rebias_wait,ct)
        recs = db.select_acqs(suffix='bc1_step', ctime=ctimes)
        if len(recs) > 0:
            offset_rebias = (rebias_wait - ct + recs[-1].ctime)*400
        else:
            offset_rebias = 0
        logger.trace(0,"Rebias offset set to %d"%offset_rebias)

    offset = max( offset_IV, offset_rebias, params.get("offset",0)*400 )
    logger.trace(0,"Total offset set to %d" %offset)

    grp = f[obs]    
    grp.attrs['_moby2_class_name'] = 'tod_flags'
    grp.attrs['_moby2_class_version'] = 1

    print("Pre-loading tod data...")
    params_loadtod = {
        'filename':obs,
        'read_data': False,
    }
    tod = moby2.scripting.get_tod(params_loadtod)
    if not end:
        end=0
    sample_count = tod.nsamps - offset + end
    if sample_count <= 0:
        print("Not enough data")
        continue
    else:
        grp.attrs['sample_count'] = sample_count

f.close()

