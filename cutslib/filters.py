import numpy as np

def gen_freqs(ndata, dt):
    """Generate frequency vector of specified length and corresponding
    to a specified total time period."""
    dn = 2      # if you like the central frequency to be negative, change dn to 1
    return 1/(ndata*dt) * np.hstack((np.arange(0, (ndata+dn)//2),
                                     np.arange(-(ndata+dn)//2+dn, 0)))

def highpass_sine2(f, fc=1.0, df=0.1):
    """highpass filter using sine square function"""
    filt = np.zeros(len(f))
    filt[np.abs(f) > fc + df/2] = 1
    sel = (np.abs(f) > fc - df/2)*(np.abs(f) < fc + df/2)
    filt[sel] = np.sin(np.pi/2/df*(np.abs(f[sel]) - fc + df/2))**2
    return filt

def lowpass_sine2(f, fc=1.0, df=0.1):
    """lowpass filter using sine square function"""
    filt = np.zeros(len(f))
    filt[np.abs(f) < fc - df/2] = 1.0
    sel = (np.abs(f) > fc - df/2)*(np.abs(f) < fc + df/2)
    filt[sel] = np.sin(np.pi/2*(1 - 1/df*(np.abs(f[sel]) - fc + df/2)))**2
    return filt

def gaussian(f, t_sigma=None, gain=1, fc=0):
    """gaussian filter"""
    sigma = 1 / (2*np.pi*t_sigma)
    return gain * np.exp(-0.5*(np.abs(f)-fc)**2/sigma**2)

