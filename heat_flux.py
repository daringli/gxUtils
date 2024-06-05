#!/usr/bin/env python

# module to compute and plot time-average heat flux
# can be imported into another script or run as a standalone script with
# > python heat_flux.py [list of .nc files]

import numpy as np
import matplotlib.pyplot as plt
import sys
from netCDF4 import Dataset

from get_Qi_t2 import get_Qi_t_nc

from matplotlib import colormaps


def get_colors(colormap_name, n):
    return colormaps[colormap_name].resampled(n)(np.linspace(0.0, 1.0, n))

def heat_flux(dirname, ispec=0, navgfac=0.5, label=None, plot=True, fig=None, Lref="a", refsp=None, color=None):
    try:
        q, t = get_Qi_t_nc(dirname)
    except:
        q = np.array([np.nan])
        t = np.array([np.nan])
    print(t)
    #species_type = data.groups['Inputs'].groups['Species'].variables['species_type'][ispec]
    #if species_type == 0:
    #    species_tag = "i"
    #elif species_type == 1:
    #    species_tag = "e"
    #if refsp == None:
    species_tag='i'

    # compute time-average and std dev
    try:
        istart_avg = int(len(t)*navgfac)
    except TypeError:
        return
    else:
        qavg = np.mean(q[istart_avg:])
        qstd = np.std(q[istart_avg:])
    if label == None:
        label = dirname
    print(r"%s: Q_%s/Q_GB = %.5g +/- %.5g" % (label, species_tag, qavg, qstd))

    # make a Q vs time plot
    if plot:
        if fig == None:
            fig = plt.figure(0)
        try:
            plt.plot(t,q,'-',label=r"%s: $Q_%s/Q_\mathrm{GB}$ = %.5g"%(label, species_tag, qavg), color=color)
        except:
            pass
        else:
            plt.ylabel(r"$Q/Q_\mathrm{GB}$")
            plt.xlabel(r"$t\ (v_{t%s}/%s)$"%(refsp, Lref))
            legend = plt.legend(loc='upper right')
            legend.set_in_layout(False)
            plt.tight_layout()

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("simdirs", nargs="+", metavar='simdir', help='Simulation directories from which to read GX outputs from.', default=['.'])
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    parser.add_argument("-l","--legend", nargs="*", action='store', metavar='label',  help='Optional list of labels for legends', default = [])
    parser.add_argument("--colormap","--cm", choices=list(colormaps), metavar='cm',  help='Name of matplotlib colormap', default = None)
    
    args = parser.parse_args()

    
    N = len(args.simdirs)
    if args.colormap is not None:
        colors=get_colors(args.colormap, N)
    else:
        colors = [None] * N
    
    print("Plotting heat fluxes.....")
    ispec=0
    Nl = len(args.legend)
    for i, d in enumerate(args.simdirs):
        color = colors[i]
        if i < Nl:
            label = args.legend[i]
        else:
            label = None
        heat_flux(d, ispec=ispec, refsp="i", label=label, color=color)
    
    plt.xlim(0)
    plt.ylim(0)
    # uncomment this line to save a PNG image of the plot
    #plt.savefig("heat_flux.png")
    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()
