#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os

from netCDF4 import Dataset

from gx_input import Gx_input

def get_Phi2_vs_z(dirname,ncfile='gx.nc'):
    ncfile = dirname + "/" + ncfile
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            t = f.variables["time"][()]
            i = 0
            if ma.is_masked(t):
                while t.mask[-1 - i] == True:
                    i = i + 1
                    
            Phi2 = f["/Spectra/Phi2zt"][()][-1 -i]
            theta = f.variables["theta"][()]
            npol = Gx_input(dirname).npol
    
        ret = (npol*np.array(theta),np.array(Phi2))
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
        theta,Phi2 = get_Phi2_vs_z(d)
        if not np.isnan(Phi2).all():
            ax.plot(theta,Phi2,color=cmap(i))
            legend.append(d)
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$|\Phi|^2$')
    ax.legend(legend)
    plt.show()
