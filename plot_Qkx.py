#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

from netCDF4 import Dataset
from decide_version import decide_version

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

  
    
def get_Qkx(d, ispec=0, navgfac=0.5, label=None, plot=False, ax=None, Lref="a", refsp=None):
    if d[-3:] == ".nc":
        d = d.rsplit('/',1)[0]
    ncfile, version = decide_version(d)

    if version == 1:
        time_str = "time"
        kx_str = "kx"
        Qkx_str = "/Spectra/Qkxst"
        
    else:
        time_str = "Grids/time"
        kx_str = "Grids/kx"
        Qkx_str = "/Diagnostics/HeatFlux_kxst"



    try: 
        data = Dataset(ncfile, mode='r')
    except FileNotFoundError:
        kx = np.array([np.nan])
        Qkx = np.array([np.nan])
        return kx, Qkx
    t = sqrt2 * data[time_str][:]
    kx = sqrt2 * data[kx_str][:]
    #print(t)
    #print(version)
    try:
        Qkxt = data[Qkx_str][:,0,:]/(2*sqrt2)
    except (KeyError, IndexError):
        print("error for '"  +ncfile +"', Skipping.")
        Qkx = np.nan * np.zeros(len(kx))
    else:
        istart_avg = int(len(t)*navgfac)
        Qkx = np.mean(Qkxt[istart_avg:], axis=0)


    
    if plot:
        if ax is None:
            fig, ax  =plt.subplots(1)
        print("!!! " + str(Qkx))
        ax.plot(kx,Qkx,'o-', markersize=0.3)
    return kx, Qkx

    

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("simdirs", nargs="+", metavar='simdir', help='Simulation directories from which to read GX outputs from.', default=['.'])
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    
    args = parser.parse_args()

    
    
    print("Plotting Qkx fluxes.....")

    fig, ax = plt.subplots(1)
    
    for d in args.simdirs:
        kx, Qkx = get_Qkx(d, ax=ax, plot=True)

    ax.set_yscale('log')
    #ax.set_xscale('log')
    
    #ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    refsp = 'i'
    ax.set_xlabel(r'$k_x \rho_{%s}$' % refsp)
    ax.set_ylabel(r"$Q/Q_\mathrm{GB}$")
    #plt.xscale('log')
    
    
    plt.tight_layout()
    plt.legend(args.simdirs)
    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()

    
