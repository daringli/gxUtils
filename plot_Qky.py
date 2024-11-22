#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

from netCDF4 import Dataset
from decide_version import decide_version


from matplotlib import colormaps


from get_Qi_t2 import get_Qi_t_nc


sqrt2 = np.sqrt(2)

def get_colors(colormap_name, n):
    return colormaps[colormap_name].resampled(n)(np.linspace(0.0, 1.0, n))


def get_Qky_max(d, ispec=0, navgfac=0.5):
    ky, Qky = get_Qky(d, ispec=ispec, navgfac=navgfac, output=output)
    return np.max(Qky)

def get_maxky(d, ispec=0, navgfac=0.5):
    ky, Qky = get_Qky(d, ispec=ispec, navgfac=navgfac)
    i = np.argmax(Qky)
    return ky[i]


    
def get_Qky(d, ispec=0, navgfac=0.5, label=None, plot=False, ax=None, Lref="a", refsp=None, color=None):

    
    if d[-3:] == ".nc":
        tmp = d.rsplit('/',1)
        d = tmp[0]
        _ncfile = tmp[1]
    else:
        _ncfile = None
    ncfile, version = decide_version(d)

    if isinstance(ncfile, list):
        # if we have restarted and got multiple output files
        # use the latest file
        ncfile = ncfile[-1]
    
    if _ncfile is not None:
        ncfile = d + "/" + _ncfile
    
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
    dky = ky[1] - ky[0]
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
            ax.plot(ky,Qky/dky,marker='.', label=label,markersize=0.3, color=color)
        except:
            pass

    # was used for debugging
    # try:
    #     q, t = get_Qi_t_nc(d)
    # except:
    #     q = np.array([np.nan])
    #     t = np.array([np.nan])
    #     print("failed to load q")
    # print(np.max(q - np.sum(Qkyt,axis=1)))
    # print(np.mean(q[istart_avg:]))
    
        
    return ky, Qky

    

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("simdirs", nargs="+", metavar='simdir', help='Simulation directories from which to read GX outputs from.', default=['.'])
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    parser.add_argument("-l","--legend", nargs="*", action='store', metavar='label',  help='Optional list of labels for legends', default = [])
    parser.add_argument("--colormap","--cm", choices=list(colormaps), metavar='cm',  help='Name of matplotlib colormap', default = None)
    parser.add_argument("--liny", action='store_const', metavar='lineary',  help='Plot with a linear scale on y-axis (default False)', default = False, const=True)
    parser.add_argument("--linx", action='store_const', metavar='linearx',  help='Plot with a linear scale on x-axis (default False)', default = False, const=True)
    
    
    args = parser.parse_args()

    liny = args.liny
    linx = args.linx
    
    N = len(args.simdirs)
    if args.colormap is not None:
        colors=get_colors(args.colormap, N)
    else:
        colors = [None] * N
        
    print("Plotting Qky fluxes.....")

    fig, ax = plt.subplots(1)

    refsp = 'i'
    ispec = 0
    Nl = len(args.legend)
    
    for i, d in enumerate(args.simdirs):
        color = colors[i]
        if i < Nl:
            label = args.legend[i]
        else:
            label = None
        get_Qky(d, ax=ax, plot=True, label=label, refsp = refsp, ispec=ispec, color=color)
        
    if not liny:
        ax.set_yscale('log')
    if not linx:
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
