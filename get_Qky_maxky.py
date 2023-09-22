#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

sqrt2 = np.sqrt(2)
from netCDF4 import Dataset

def get_maxky(d, fname = "gx.nc"):
    fname = d + "/" + fname
    data = Dataset(fname, mode='r')
    t = sqrt2 * data.variables['time'][:]
    ky = sqrt2 * data.variables['ky'][:]
    Qkyt = data.groups['Spectra'].variables['Qkyst'][:,0,:]/(2*sqrt2)
    Qky = np.mean(Qkyt[int(len(t)/2):], axis=0)

    i = np.argmax(Qky)
    return ky[i]
