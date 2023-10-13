#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

from netCDF4 import Dataset
from decide_version import decide_version


sqrt2 = np.sqrt(2)


def get_Qky_max(d, ispec=0, navgfac=0.5):
    ky, Qky = get_Qky(d, ispec=ispec, navgfac=navgfac, output=output)
    return np.max(Qky)

def get_maxky(d, ispec=0, navgfac=0.5):
    ky, Qky = get_Qky(d, ispec=ispec, navgfac=navgfac)
    i = np.argmax(Qky)
    return ky[i]


    
def get_Qky(d, ispec=0, navgfac=0.5, label=None, plot=False, ax=None, Lref="a", refsp=None):

    
    if d[-3:] == ".nc":
        d = d.rsplit('/',1)[0]
    ncfile, version = decide_version(d)
    
    if version == 1:
        time_str = "time"
        ky_str = "ky"
        Qky_str = "/Spectra/Qkyst"
        
    else:
        time_str = "Grids/time"
        ky_str = "Grids/ky"
        Qky_str = "/Diagnostics/HeatFlux_kyst"
    
    data = Dataset(ncfile, mode='r')
    t = sqrt2 * data[time_str][:]
    ky = sqrt2 * data[ky_str][:]
    #print(t)
    #print(version)
    try:
        Qkyt = data[Qky_str][:,0,:]/(2*sqrt2)
    except KeyError:
        print("error for '"  +ncfile +"', Skipping.")
        Qky = np.nan * np.zeros(len(ky))
    else:
        istart_avg = int(len(t)*navgfac)
        Qky = np.mean(Qkyt[istart_avg:], axis=0)


    
    if plot:
        if ax is None:
            fig, ax  =plt.subplots(1)
        print("!!! " + str(Qky))
        ax.plot(ky,Qky,'o-')
    
        
    return ky, Qky

    

if __name__ == "__main__":
    
    print("Plotting Qky fluxes.....")

    fig, ax = plt.subplots(1)
    
    for fname in sys.argv[1:]:
        ky, Qky = get_Qky(fname, ax=ax, plot=True)
        
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    refsp = 'i'
    ax.set_xlabel(r'$k_y \rho_{%s}$' % refsp)
    ax.set_ylabel(r"$Q/Q_\mathrm{GB}$")
    #plt.xscale('log')
    
    
    plt.tight_layout()
    plt.legend(sys.argv[1:])
    plt.show()
