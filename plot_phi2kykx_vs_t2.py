#!/usr/bin/env python

from netcdf_util import netcdf_file


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

from plot_phi2kykx_vs_t import stellaClickHandler

import sys

BIGNUM = 1e200

class gxClickHandler(stellaClickHandler):
    def __init__(self, gx_output = 'gx.nc', full_kx = False, Ntheta0 = 0, sqrt2 = np.sqrt(2)):
        self.imax_pressed = False
        self.imin_pressed = False
        self.maxlines = []
        self.minlines = []
        self.imin = 0
        self.imax = None
        self.full_kx = full_kx
        self.Ntheta0 = Ntheta0
        with netcdf_file(gx_output,'r',mmap=False) as f:
            self.phi2_vs_kxky  = f['/Spectra/Phi2kxkyt'][()] # time, ky, kx
            self.ky = sqrt2*f.variables['ky'][()]
            self.kx = sqrt2*f.variables['kx'][()]
            self.t = sqrt2*f.variables['time'][()]
            self.shat = f['/Geometry/shat'][()]

        self.Nkx = len(self.kx)
        self.Nky = len(self.ky)
        self.Nt = len(self.t)
        if not self.full_kx:
            self.phi2s = np.transpose(self.phi2_vs_kxky)[self.Nkx//2:,:,:] # kx, ky, time
            self.kx = self.kx[self.Nkx//2:]
            self.Nkx = self.Nkx//2 + 1
            self.kx0_index = 0
        else:
            self.phi2s = np.transpose(self.phi2_vs_kxky)[:,:,:] # kx, ky, time
            self.kx0_index = self.Nkx//2 + 1
        ##self.phi2s = np.reshape(self.phi2s,(self.Nky*self.Nkx,self.Nt))
        self.gamma = np.zeros((self.Nkx,self.Nky))
        self.iinfs = np.full((self.Nkx,self.Nky),None)
        
        for ikx, phi2_x in enumerate(self.phi2s):
            for iky, phi2 in enumerate(phi2_x):
                iinf = np.where(phi2 > BIGNUM)[0]
                if len(iinf) == 0:
                    iinf = np.where(np.isnan(phi2))[0]
                if len(iinf) > 0:
                    self.iinfs[ikx,iky] = iinf[0]
        # fix this
 

if __name__=="__main__":
    import sys
    argc = len(sys.argv)

    sqrt2 = np.sqrt(2)
    if argc > 1:
        dir = sys.argv[1]
        if argc > 2:
            sqrt2 = 1.0
            
            
    else:
        dir = '.'
    ch = gxClickHandler(dir + '/gx.nc', full_kx = False, Ntheta0 = 3, sqrt2=sqrt2)
    ch.plot()
