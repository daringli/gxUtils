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
    
    if isinstance(ncfile, list):
        # if we have restarted and got multiple output files
        # use the latest file
        ncfile = ncfile[-1]
    
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

    if np.all(np.isnan(Qkx)):
        print("All nan for '"  +ncfile +"', Skipping.")
        return kx, Qkx
    
    
    if plot:
        if ax is None:
            fig, ax  =plt.subplots(1)
        if label == None:
            label = d
        try:
            ax.plot(kx,Qkx,marker='.', markersize=0.3, label=label)
        except:
            pass
            
    return kx, Qkx

    

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("simdirs", nargs="+", metavar='simdir', help='Simulation directories from which to read GX outputs from.', default=['.'])
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    parser.add_argument("-l","--legend", nargs="*", action='store', metavar='label',  help='Optional list of labels for legends', default = [])
    
    args = parser.parse_args()
    Nl = len(args.legend)
    
    
    refsp = 'i'
    ispec = 0
    print("Plotting Qkx fluxes.....")

    fig, ax = plt.subplots(1)
    
    for i, d in enumerate(args.simdirs):
        if i < Nl:
            label = args.legend[i]
        else:
            label = None
        kx, Qkx = get_Qkx(d, ax=ax, plot=True, label = label, ispec = ispec, refsp = refsp)

    ax.set_yscale('log')
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$k_x \rho_{%s}$' % refsp)
    ax.set_ylabel(r"$Q/Q_\mathrm{GB}$")
    ax.set_ylim(bottom=0)
    legend = plt.legend(loc='upper right')
    legend.set_in_layout(False)
    plt.tight_layout()

    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()

    
