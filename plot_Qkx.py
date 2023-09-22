#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

from netCDF4 import Dataset

sqrt2 = np.sqrt(2)

def get_Qkx_max(d, ispec=0, navgfac=0.5, output = "gx.nc"):
    kx, Qkx = get_Qkx(d, ispec=ispec, navgfac=navgfac, output=output)
    return np.max(Qkx)

def get_maxkx(d, ispec=0, navgfac=0.5, output = "gx.nc"):
    kx, Qkx = get_Qkx(d, ispec=ispec, navgfac=navgfac, output=output)
    i = np.argmax(Qkx)
    return kx[i]

def get_kx_width(d, factor = 0.01, ispec=0, navgfac=0.5, output = "gx.nc"):
    kx, Qkx = get_Qkx(d, ispec=ispec, navgfac=navgfac, output=output)
    print(kx)
    maxQkx = np.max(Qkx)
    i = np.argmin(np.fabs(Qkx/maxQkx - factor))
    return kx[i]

  
    
def get_Qkx(d, ispec=0, navgfac=0.5, label=None, plot=False, ax=None, Lref="a", refsp=None, output="gx.nc"):
    if d[-3:] != ".nc":
        # default filename for outputs
        fname = d + "/" + output
    else:
        fname = d

    data = Dataset(fname, mode='r')
    t = sqrt2 * data.variables['time'][:]
    kx = sqrt2 * data.variables['kx'][:]
    try:
        Qkxt = data.groups['Spectra'].variables['Qkxst'][:,0,:]/(2*sqrt2)
    except KeyError:
        print("error for '"  +fname +"', Skipping.")
        Qkx = np.nan * np.zeros(len(kx))
    else:
        istart_avg = int(len(t)*navgfac)
        Qkx = np.mean(Qkxt[istart_avg:], axis=0)
        
    if plot:
        if ax is None:
            fig, ax  =plt.subplots(1)
        ax.plot(kx,Qkx,'o-')
    
        
    return kx, Qkx

    

if __name__ == "__main__":
    
    print("Plotting Qkx fluxes.....")

    fig, ax = plt.subplots(1)
    
    for fname in sys.argv[1:]:
        if fname[-3:] != ".nc":
            # default filename for outputs
            fname = fname + "/gx.nc"
        try:
            kx, Qkx = get_Qkx(fname, ax=ax, plot=True)
        except:
            print(' usage: python plot_Qkx.py [list of .nc files or dirs with gx.nc files]')

    ax.set_yscale('log')
    #ax.set_xscale('log')
    
    #ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    refsp = 'i'
    ax.set_xlabel(r'$k_y \rho_{%s}$' % refsp)
    ax.set_ylabel(r"$Q/Q_\mathrm{GB}$")
    #plt.xscale('log')
    
    
    plt.tight_layout()
    plt.legend(sys.argv[1:])
    plt.show()
