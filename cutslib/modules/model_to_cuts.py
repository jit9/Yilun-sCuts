"""This module takes in a machine learning based model (from mlpipe)
and generate cuts+cal based on its prediction."""

import pickle, os.path as op, numpy as np
import moby2


class Module:
    def __init__(self, config):
        self.model_file = config.get("model_file")
        self.tag_out = config.get("tag_out", None)

    def run(self, p):
        model_file = self.model_file
        tag_out = self.tag_out
        # load pickle file
        pickle_file = p.o.pickle_file
        with open(pickle_file, "rb") as f:
            data = pickle.load(f)
        # get data shape
        ndet, ntod = data['sel'].shape
        # load model file
        with open(model_file, "rb") as f:
            model = pickle.load(f)
        # gather training data
        # loop over tod
        for i in np.arange(ntod):
            obs = data['name'][i]
            features = []
            print("TOD: %d/%d" % (i,len(data['name'])))
            # loop over keys
            for k in model.features:
                # (ndet, ntod)
                feat = data[k]
                if np.ndim(feat) == 2:
                    features.append(feat[:,i][:,None])
                # (ndet or ntod)
                elif np.ndim(feat) == 1:
                    if len(feat) == ndet:
                        features.append(feat[:,None])
                    elif len(feat) == ntod:
                        features.append(np.ones((ndet,1))*feat[i])
                    else:
                        raise ValueError("feat not understood")
                else:
                    raise ValueError("feat not understood")
            # get prediction
            features = np.hstack(features)
            pred = model.predict(features).astype(bool)
            # generate cuts based on prediction
            try:
                tod = moby2.scripting.get_tod({'filename': obs,
                                               'read_data': False})
            except IOError as e:
                print("Failed to read tod, skipping...")
                continue
            cuts = moby2.TODCuts.for_tod(tod, assign=False)
            cuts.set_always_cut(~pred)
            # TODO: merge in all other cuts
            # write out
            if tag_out:
                depot = moby2.util.Depot(p.depot)
                depot.write_object(cuts, tod=tod,
                                   force=True, tag=tag_out)
