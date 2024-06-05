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

    if isinstance(ncfile, list):
        # if we have restarted and got multiple output files
        # use the latest file
        ncfile = ncfile[-1]
    
    if version == 1:
        time_str = "time"
        ky_str = "ky"
        Qky_str = "/Spectra/Qkyst"
        
    else:
        time_str = "Grids/time"
        ky_str = "Grids/ky"
        Qky_str = "/Diagnostics/HeatFlux_kyst"

    try: 
        data = Dataset(ncfile, mode='r')
    except FileNotFoundError:
        ky = np.array([np.nan])
        Qky = np.array([np.nan])
        return ky, Qky

    if label == None:
        label = d
        
    t = sqrt2 * data[time_str][:]
    ky = sqrt2 * data[ky_str][:]
    #print(t)
    #print(version)
    try:
        Qkyt = data[Qky_str][:,0,:]/(2*sqrt2)
    except (KeyError, IndexError):
        print("error for '"  +ncfile +"', Skipping.")
        Qky = np.nan * np.zeros(len(ky))
        return ky, Qky
    else:
        istart_avg = int(len(t)*navgfac)
        Qky = np.mean(Qkyt[istart_avg:], axis=0)

    if np.all(np.isnan(Qky)):
        print("All nan for '"  +ncfile +"', Skipping.")
        return ky, Qky
    
    if plot:
        if ax is None:
            fig, ax  =plt.subplots(1)
        try:
            ax.plot(ky,Qky,marker='.', label=label,markersize=0.3)
        except:
            pass
            
        
    return ky, Qky

    

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("simdirs", nargs="+", metavar='simdir', help='Simulation directories from which to read GX outputs from.', default=['.'])
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    parser.add_argument("-l","--legend", nargs="*", action='store', metavar='label',  help='Optional list of labels for legends', default = [])
    
    print("Plotting Qky fluxes.....")

    fig, ax = plt.subplots(1)

    refsp = 'i'
    ispec = 0
    args = parser.parse_args()
    Nl = len(args.legend)
    
    for i, d in enumerate(args.simdirs):
        if i < Nl:
            label = args.legend[i]
        else:
            label = None
        get_Qky(d, ax=ax, plot=True, label=label, refsp = refsp, ispec=ispec)
        
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$k_y \rho_{%s}$' % refsp)
    ax.set_ylabel(r"$Q/Q_\mathrm{GB}$")
    legend = plt.legend(loc='upper right')
    legend.set_in_layout(False)
    plt.tight_layout()
    
    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()
