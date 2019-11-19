import joblib
import multiprocessing
import moby2
import sys

params = moby2.util.MobyDict.from_file(sys.argv[1])

numcores = params.get( 'ncores', 8 )
tl = moby2.util.database.TODList.from_file(params['source_list'])
input_cuts_tag = params['input_cuts']
input_partial_tag = params['input_partial']
output_cuts_tag = params['output_cuts']
depot_path = params.get( 'depot', '/data/actpol/depot' )
depot = moby2.util.Depot(depot_path)

for obs in tl:
def replace_pcuts(obs, depot, params):
    input_cuts_tag = params['input_cuts']
    input_partial_tag = params['input_partial']
    output_cuts_tag = params['output_cuts']
    depot_path = params.get( 'depot', '/data/actpol/depot' )
    cuts = moby2.scripting.get_cuts({'depot':depot_path,'tag':input_cuts_tag}, tod=obs)
    pcuts = moby2.scripting.get_cuts({'depot':depot_path,'tag':input_cuts_tag+'_partial'}, tod=obs)
    pcuts_new = moby2.scripting.get_cuts({'depot':depot_path,'tag':input_partial_tag}, tod=obs)
    cuts_new = cuts.copy()
    ld = cuts.get_uncut()
    for d in ld:
        cuts_mask = cuts.cuts[d].get_mask()
        pcuts_mask = pcuts.cuts[d].get_mask()
        pcuts_new_mask = pcuts_new.cuts[d].get_mask()
        cuts_new.cuts[d] = np.logical_or( cuts_mask != pcuts_mask, pcuts_new_mask )
    depot.write_object(cuts_new, tag=output_cuts_tag, tod=obs, make_dirs=True, force=True)
    return 0


joblib.Parallel(n_jobs=numcores)(
    joblib.delayed(replace_pcuts)(obs, depot, params) for obs in tl)
