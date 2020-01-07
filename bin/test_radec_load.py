import moby2
import sys
import matplotlib.pyplot as plt
import numpy as np

tod_name = sys.argv[-1]
tod = moby2.scripting.get_tod({'filename': tod_name, 'repair_pointing': True, 'read_data': False})
error = np.sum(np.logical_or(tod.alt > np.pi/2, tod.alt < 0)) > 0
if error:
    print("Alt falls outside range of 0 and pi/2")

if 0:
    plt.plot(tod.az, tod.alt)
    plt.xlabel("Az")
    plt.ylabel("Alt")
    plt.savefig("AzAlt.png")

ra, dec = moby2.pointing.get_coords(tod.ctime,tod.az,tod.alt)
