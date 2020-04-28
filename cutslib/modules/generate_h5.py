"""Convert pickle file to h5 file"""
import os.path as op, h5py, pickle, numpy as np
from cutslib.util import pickle_load

class Module:
    def __init__(self, config):
        """Module to convert pickle file to h5 file
        Args:
            randomize: if tods should be chosen randomly
            limit: number of tods to export to h5
            pickle_file_label: if want to use a different pkl file for label
        """
        self.outdir = config.get("outdir",".")
        self.train_portion = config.getfloat("train_portion", 0.6)
        self.validate_portion = config.getfloat("validate_portion", 0.2)
        self.test_portion = config.getfloat("test_portion", 0.2)
        self.randomize = config.getboolean("randomize", True)
        self.limit = config.getint("limit", None)
        self.pickle_file_label = config.get("pickle_file_label", None)

    def run(self, p):
        outdir = self.outdir
        train_portion = self.train_portion
        validate_portion = self.validate_portion
        test_portion = self.test_portion
        randomize = self.randomize
        limit = self.limit
        pickle_file_label = self.pickle_file_label
        outfile = op.join(outdir, f"{p.tag}.h5")
        # specify keys in the data of interests
        # data per tod per det
        keys_ptd = ['psel', 'resp', 'resp_sel', 'cal', 'corrLive',
                    'DELive', 'MFELive', 'rmsLive', 'skewLive',
                    'kurtLive', 'normLive', 'gainLive', 'jumpLive']
        # data per det
        keys_pd = ['ff', 'stable']
        # data per tod
        keys_pt = ['alt', 'pwv']
        # load pickle file
        pickle_file = p.o.pickle_file
        data = pickle_load(p.o.pickle_file)
        # if a different pickle file is used for label
        # load it now
        if pickle_file_label:
            data_label = pickle_load(pickle_file_label)
        # get tod names
        obs = np.array(data['name'])
        # if a differet pickle file is used for label
        # find their intersection
        if pickle_file_label:
            obs_label = np.array(data_label['name'])
            obs = np.array(list(set(obs).intersection(set(obs_label))))
        # whether to randomize tod lists
        if randomize:
            np.random.shuffle(obs)
        # restrict number of tods to output
        if limit:
            limit = min(limit, len(obs))
            obs = obs[:limit]
        nobs = len(obs)
        # find their corresponding indices
        indices = np.array([data['name'].index(e) for e in obs])
        if pickle_file_label:
            indices_label = np.array([data_label['name'].index(e) for \
                                      e in obs])
        # get train,test,validate indices
        train_idx = int(nobs*train_portion)
        validate_idx = train_idx + int(nobs*validate_portion)
        train_indices = indices[:train_idx]
        validate_indices = indices[train_idx:validate_idx]
        test_indices = indices[validate_idx:]
        if pickle_file_label:
            train_indices_label = indices_label[:train_idx]
            validate_indices_label = indices_label[train_idx:validate_idx]
            test_indices_label = indices_label[validate_idx:]
        # get live dets
        dets = np.where(data['live'])[0]
        # create output h5 file if it doesn't exist
        if op.isfile(outfile):
            print("File %s exists, updating instead" % outfile)
            hf = h5py.File(outfile, 'a')
        else:
            # file doesn't exist
            print("Create file: %s" % outfile)
            hf = h5py.File(outfile, 'w')
        # fill up each group
        indices_list = [train_indices, validate_indices, test_indices]
        if pickle_file_label:
            indices_list_label = [train_indices_label,
                                  validate_indices_label,
                                  test_indices_label]
        groups = ['train', 'validate', 'test']
        for ig in range(len(groups)):
            idx = indices_list[ig]
            gname = groups[ig]
            if pickle_file_label:
                idx_label = indices_list_label[ig]
            # create or get existing group
            try:
                print("Create group %s" % gname)
                g = hf.create_group(gname)
            except ValueError:
                print("Group %s exists, update instead" % gname)
                g = hf[gname]
            # loop over indices: ii represents i of indices
            for ii in range(len(idx)):
                # get corresponding index in data['name]
                i = idx[ii]
                if pickle_file_label:
                    # get corresponding index in data_label['name']
                    i_label = idx_label[ii]
                    # make sure we are doing the right thing
                    try:
                        assert data['name'][i] == data_label['name'][i_label]
                    except AssertionError:
                        print(data['name'][i])
                        print(data_label['name'][i_label])
                        import sys;sys.exit(-1)
                print("TOD: %d/%d"%(ii+1,len(idx)))
                for det_uid in dets:
                    # name of dataset
                    name_ = "%s.%d" % (data['name'][i],det_uid)
                    # create or get dataset
                    try:
                        d = g.create_dataset(name_, dtype='float')
                    except RuntimeError:
                        d = g[name_]
                    # store all keys
                    for k in keys_ptd:
                        d.attrs[k] = data[k][det_uid, i]
                    for k in keys_pt:
                        d.attrs[k] = data[k][i]
                    for k in keys_pd:
                        d.attrs[k] = data[k][det_uid]
                    # store label
                    # depends on which pickle file we use
                    if pickle_file_label:
                        d.attrs['label'] = data_label['sel'][det_uid,
                                                             i_label]
                    else:
                        d.attrs['label'] = data['sel'][det_uid, i]
