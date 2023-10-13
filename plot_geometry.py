#!/usr/bin/env python

import numpy as np
from glob import glob
import os

from netCDF4 import Dataset

from decide_version import decide_version

def get_geometry(dirname):
    ncfile, version = decide_version(dirname)
    if ncfile is None:
        print("GX output does not exist in '" + dirname + "'")
    
    #ncfile = dirname + "/" + ncfile
    if version == 1:
        theta_str = "theta"
    else:
        theta_str = "/Grids/theta"

    ret = []
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            theta = f[theta_str][()]
            ntheta = len(theta)
            outputs = ["bmag", "gradpar", "gbdrift", "cvdrift", "cvdrift0", "gds2", "gds21", "gds22"]
            for output in outputs:
                y = f["/Geometry/" + output][()]
                print(output)
                if output == 'gradpar':
                    y = np.array([y] * ntheta)
                ret.append(y)
        ret = [theta] + ret
        ret = tuple(ret)
    else:
        print("File does not exist '" + ncfile + "'")
        ret = (np.nan,) * 9
    return ret

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    #fig, ax = plt.subplots()
    fig, axes = plt.subplots(4,2, sharex=True)
    axes = axes.flatten()
    Ndirs = len(ds)
    cmap = plt.get_cmap("brg",Ndirs)
    legend = []

    
    
    for i,d in enumerate(ds):
        theta,bmag, gradpar, gbdrift, cvdrift, cvdrift0, gds2, gds21, gds22 = get_geometry(d)
        outputs = [bmag, gradpar, gbdrift, cvdrift, cvdrift0, gds2, gds21, gds22]
        for j,o in enumerate(outputs):
            if not np.isnan(o).all():
                print(j)
                print(o.shape)
                axes[j].plot(theta,o,color=cmap(i))
                axes[j].set_title("")
        legend.append(d)

    axes[-2].set_xlabel(r"$\theta/\pi$")
    axes[-1].set_xlabel(r"$\theta/\pi$")
    #ax.set_ylabel(r'bmag')
    axes[1].legend(legend)
    plt.show()
