#!/usr/bin/env python

from plot_Phi2_vs_z import get_Phi2_vs_z

import os
import numpy as np
from netCDF4 import Dataset


def get_space_filling(dirname, ncfile='gx.nc'):

    theta, Phi2 = get_Phi2_vs_z(dirname,ncfile)

    
    ncfile = dirname + "/" + ncfile
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            gradpar = f["/Geometry/gradpar"][()]  # \nabla_\| z = \nabla_\| \theta
            bmag = f["/Geometry/bmag"][()]

        I = np.trapz(Phi2/(gradpar * bmag), theta)
        I2 = np.trapz(1/(gradpar * bmag), theta)
        ret = np.sqrt(I/I2/np.max(Phi2))
    
    
    else:
        print("File does not exist '" + ncfile + "'")
        ret = np.nan

    return ret

if __name__=="__main__":
    import sys
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    for i,d in enumerate(ds):
        ret = get_space_filling(d)
        print(d)
        print(ret)
