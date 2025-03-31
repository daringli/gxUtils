#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os

from netCDF4 import Dataset


from decide_version import decide_version

from gx_input import Gx_input

def get_Phi2kxky_vs_z(dirname):

    ncfile, version = decide_version(dirname)
    ncfile = dirname + '/gx.big.nc'
    # Phi(time, ky, kx, theta, ri)
    
    if ncfile is None:
        print("GX output does not exist in '" + dirname + "'")
        return (np.nan,np.nan)

    time_str = "Grids/time"
    theta_str = "Grids/theta"
    ky_str = "Grids/ky"
    kx_str = "Grids/kx"
    
    Phi_str = "/Diagnostics/Phi"
            
        
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            t = f[time_str][()]
            theta = f[theta_str][()]
            kx = f[kx_str][()]
            ky = f[ky_str][()]
            
            i = 0
            if ma.is_masked(t):
                while t.mask[-1 - i] == True:
                    i = i + 1
                    
            Phi = f[Phi_str][()][:-1 -i]
            t = t[:-1-i]
            # sum real and imaginary part
            Phi2 = np.sum(Phi**2,axis=-1) 
            print(Phi2.shape)
            print(len(t))
            print(len(ky))
            print(len(kx))
            print(len(theta))
            
            
            Nt = Phi2.shape[0]
            bignum = 1e20
            i = 0
            for i in range(Nt):
                if np.any(Phi2[i,0] > bignum):
                    Phi2 = Phi2[i-1,:]
                    break
            else:
                Phi2 = Phi2[Nt-1,:]
                    
            gi = Gx_input(dirname)
            if gi.geo_option == 'vmec':
                npol = gi.npol
            else:
                npol = 1
                
        ret = (np.array(ky), np.array(kx), npol*np.array(theta), np.array(Phi2))
    else:
        print("File does not exist '" + ncfile + "'")
        ret = (np.nan,np.nan,np.nan,np.nan)
    return ret

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) ==1:
        d = sys.argv[1]
    else:
        d = '.'

    fig, ax = plt.subplots()

    legend = []
    
    ky,kx, theta,Phi2 = get_Phi2kxky_vs_z(d)
    Nky = len(ky)
    Nkx = len(kx)

    Nkx_plot = 3
    Nky_plot = 10 
    plot_kx = np.linspace(0, np.max(kx), Nkx_plot)
    plot_ky = np.linspace(ky[1], np.max(ky), Nky_plot)
    
    cmap = plt.get_cmap("brg", Nky_plot)
    linestyles = ['solid','dashed','dotted']
    if not np.isnan(Phi2).all():
        for i in range(Nky_plot):
            iky = np.argmin(np.abs(ky  - plot_ky[i]))
            for j in range(Nkx_plot):
                ikx = np.argmin(np.abs(kx  - plot_kx[j]))
                print("(ky, kx) = " + str((ky[iky], kx[ikx])))
                ax.plot(theta,Phi2[iky,ikx], color=cmap(i), linestyle=linestyles[j])
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$|\Phi|^2$')
    plt.savefig('Phi_vs_z.png')
