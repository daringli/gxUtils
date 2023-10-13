#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os

from netCDF4 import Dataset


from decide_version import decide_version

from gx_input import Gx_input

def get_Phi2_vs_z(dirname):

    ncfile, version = decide_version(dirname)
    if ncfile is None:
        print("GX output does not exist in '" + dirname + "'")
        return (np.nan,np.nan)

    if version == 1:
        time_str = "time"
        theta_str = "theta"
        Phi2_str = "/Spectra/Phi2zt"
        
    else:
        time_str = "Grids/time"
        theta_str = "Grids/theta"
        Phi2_str = "/Diagnostics/Phi2_zt"
            
        
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            t = f[time_str][()]
            i = 0
            if ma.is_masked(t):
                while t.mask[-1 - i] == True:
                    i = i + 1
                    
            Phi2 = f[Phi2_str][()][:-1 -i]
            Nt = Phi2.shape[0]
            bignum = 1e20
            i = 0
            for i in range(Nt):
                if Phi2[i,0] > bignum:
                    Phi2 = Phi2[i-1,:]
                    break
            else:
                Phi2 = Phi2[Nt-1,:]
                    
            theta = f[theta_str][()]
            gi = Gx_input(dirname)
            if gi.geo_option == 'vmec':
                npol = gi.npol
            else:
                npol = 1
                
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
