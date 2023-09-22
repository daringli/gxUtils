#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

from netCDF4 import Dataset

sqrt2 = np.sqrt(2)


def get_Qky(d, ispec=0, navgfac=0.5, label=None, plot=True, ax=None, Lref="a", refsp=None, output="gx.nc"):

    if d[-3:] != ".nc":
        # default filename for outputs
        fname = fname + "/" + output
    else:
        fname = d

    data = Dataset(fname, mode='r')
    t = sqrt2 * data.variables['time'][:]
    ky = sqrt2 * data.variables['ky'][:]
    try:
        Qkyt = data.groups['Spectra'].variables['Qkyst'][:,0,:]/(2*sqrt2)
    except KeyError:
        print("error for '"  +fname +"', Skipping.")
        Qkyt = np.nan * np.zeros(len(ky))

    istart_avg = int(len(t)*navgfac)
    Qky = np.mean(Qkyt[istart_avg:], axis=0)
        
    if plot:
        if ax is None:
            fig, ax  =plt.subplots(1)
        ax.plot(ky,Qky,'o-')
    
        
    return ky, Qky

    

if __name__ == "__main__":
    
    print("Plotting Qky fluxes.....")

    fig, ax = plt.subplots(1)
    
    for fname in sys.argv[1:]:
        if fname[-3:] != ".nc":
            # default filename for outputs
            fname = fname + "/gx.nc"
        try:
            ky, Qky = get_Qky(fname, ax=ax, plot=True)
        except:
            print(' usage: python plot-Qky.py [list of .nc files or dirs with gx.nc files]')

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
