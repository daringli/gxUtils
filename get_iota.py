from netCDF4 import Dataset
import numpy as np


from decide_version import decide_version

def get_iota(dirname):
    
    ncfile, version = decide_version(dirname)

    if version == 1:
        with Dataset(ncfile, mode='r') as f:
            q = f.groups['Geometry'].variables['q'][()]

    elif version == 2:
        with Dataset(ncfile, mode='r') as f:
            q = f.groups['Geometry'].variables['q'][()]
    else:
        q = np.nan
    return 1/q

def get_shat(dirname):
    ncfile, version = decide_version(dirname)

    if version == 1:
        with Dataset(ncfile, mode='r') as f:
            shat = f.groups['Geometry'].variables['shat'][()]

    elif version == 2:
        with Dataset(ncfile, mode='r') as f:
            shat = f.groups['Geometry'].variables['shat'][()]
    else:
        shat = np.nan
    return shat
