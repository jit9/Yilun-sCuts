"""This script aims to produce a 2d histogram of rms and input gain
"""

import pickle, numpy as np, os.path as op
import matplotlib.pyplot as plt

def init(config):
    global xbins, ybins, xmin, xmax, ymin, ymax
    xbins = config.getint("xbins", 100)
    ybins = config.getint("ybins", 100)
    xmin = config.getfloat("xmin", 0)
    xmax = config.getfloat("xmax", 1e-15)
    ymin = config.getfloat("ymin", 0)
    ymax = config.getfloat("ymax", 4e-16)

def run(p):
    global xbins, ybins, xmin, xmax, ymin, ymax
    # load pickle file
    pickle_file = p.o.pickle_file
    with open(pickle_file, "r") as f:
        data = pickle.load(f)
    # calibrate rmsLive
    input_gain = data['resp'] * data['ff'][:,np.newaxis]
    data['rmsLive'] *= input_gain
    data['normLive'] *= input_gain
    data['MFELive'] *= input_gain
    data['DELive'] *= input_gain
    data['jumpLive'] *= input_gain
    # make 2d histogram
    # x: rms, y: input gain
    # get bin edges
    x_edges = np.linspace(xmin, xmax, xbins+1)
    y_edges = np.linspace(ymin, ymax, ybins+1)
    X, Y = np.meshgrid(x_edges, y_edges)
    sel = data['sel'].astype('bool')
    rms = data['rmsLive'][sel]
    gain = input_gain[sel]
    H, _, __ = np.histogram2d(rms, gain, bins=(x_edges,y_edges))
    # make individual hists for reference
    # hist for gain
    plt.hist(gain, bins=ybins)
    plt.xlabel("gain")
    outfile = op.join(p.o.root, "gain_hist.png")
    print("Saving plot: %s" % outfile)
    plt.savefig(outfile)
    plt.clf()
    # hist for rms
    plt.hist(rms, bins=xbins)
    plt.xlabel("rms")
    outfile = op.join(p.o.root, 'rms_hist.png')
    print("Saving plot: %s" % outfile)
    plt.savefig(outfile)
    plt.clf()
    # hist for bias-step
    plt.hist(data['resp'][sel], bins=xbins)
    plt.xlabel("bias-step")
    outfile = op.join(p.o.root, 'bs_hist.png')
    print("Saving plot: %s" % outfile)
    plt.savefig(outfile)
    plt.clf()
    # hist for flatfield
    plt.hist(data['ff'], bins=xbins)
    plt.xlabel("flatfield")
    outfile = op.join(p.o.root, 'ff_hist.png')
    print("Saving plot: %s" % outfile)
    plt.savefig(outfile)
    plt.clf()
    # start to produce the 2d histogram
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(X,Y,H, cmap="hot")
    plt.colorbar(mesh)
    plt.xlabel("rms")
    plt.ylabel("input gain")
    outfile = op.join(p.o.root, "rms_gain.png")
    print("Saving plot: %s" % outfile)
    plt.savefig(outfile)
