#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os

from netCDF4 import Dataset

def get_ky_nc(dirname,ncfile='gx.nc', tfit = 85):
    ncfile = dirname + "/" + ncfile
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            #f.replace('--', np.nan)
            ky = np.sqrt(2)*f.variables["ky"][()]
    return ky

