"""Convert pickle file to h5 file"""
import os.path as op, h5py, pickle, numpy as np


class Module:
    def __init__(self, config):
        self.outfile = config.get("outfile","dataset.h5")
        self.train_portion = config.getfloat("train_portion", 0.6)
        self.validate_portion = config.getfloat("validate_portion", 0.2)
        self.test_portion = config.getfloat("test_portion", 0.2)
        self.randomize = config.getboolean("randomize", True)
        self.limit = config.getint("limit", None)

    def run(self, p):
        outfile = self.outfile
        train_portion = self.train_portion
        validate_portion = self.validate_portion
        test_portion = self.test_portion
        randomize = self.randomize
        limit = self.limit
        # specify keys in the data of interests
        # data per tod per det
        keys_ptd = ['psel', 'resp', 'respSel', 'cal', 'corrLive', 'DELive', 'MFELive', \
                   'rmsLive', 'skewLive', 'kurtLive', 'normLive', 'gainLive', 'jumpLive']
        # data per det
        keys_pd = ['ff', 'stable']
        # data per tod
        keys_pt = ['alt', 'pwv']
        # load pickle file
        pickle_file = p.o.pickle_file
        with open(p.o.pickle_file, "rb") as f:
            data = pickle.load(f)
        # get tod names
        obs = np.array(data['name'])
        if limit:
            obs = obs[:limit]
        nobs = len(obs)
        # whether to randomize tod lists
        indices = np.arange(nobs)
        if randomize:
            np.random.shuffle(indices)
        # get train,test,validate indices
        train_idx = int(nobs*train_portion)
        validate_idx = train_idx + int(nobs*validate_portion)
        train_indices = indices[:train_idx]
        validate_indices = indices[train_idx:validate_idx]
        test_indices = indices[validate_idx:]
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
        groups = ['train', 'validate', 'test']
        for idx, gname in zip(indices_list, groups):
            # create or get existing group
            try:
                print("Create group %s" % gname)
                g = hf.create_group(gname)
            except ValueError:
                print("Group %s exists, update instead" % gname)
                g = hf[gname]
            for c,i in enumerate(idx):
                print("TOD: %d/%d"%(c+1,len(idx)))
                for det_uid in dets:
                    # name of dataset
                    name_ = "%s.%d" % (obs[i],det_uid)
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
                    d.attrs['label'] = data['sel'][det_uid, i]
