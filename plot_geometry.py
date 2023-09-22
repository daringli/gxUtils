#!/usr/bin/env python

import numpy as np
from glob import glob
import os

from netCDF4 import Dataset

def get_geometry(dirname,ncfile='gx.nc'):
    ncfile = dirname + "/" + ncfile
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            bmag = f["/Geometry/bmag"][()]
            theta = f.variables["theta"][()]
        ret = (np.array(theta),np.array(bmag))
    else:
        print("File does not exist '" + ncfile + "'")
        ret = (np.nan,np.nan)
    return ret

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    fig, ax = plt.subplots()

    Ndirs = len(ds)
    cmap = plt.get_cmap("brg",Ndirs)
    legend = []
    
    for i,d in enumerate(ds):    
        theta,bmag = get_geometry(d)
        if not np.isnan(bmag).all():
            ax.plot(theta,bmag,color=cmap(i))
            legend.append(d)
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'bmag')
    ax.legend(legend)
    plt.show()
