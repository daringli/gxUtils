#!/usr/bin/env python

# module to compute and plot time-average heat flux
# can be imported into another script or run as a standalone script with
# > python heat_flux.py [list of .nc files]

import numpy as np
import matplotlib.pyplot as plt
import sys
from netCDF4 import Dataset

from get_Qi_t2 import get_Qi_t_nc


def heat_flux(dirname, ispec=0, navgfac=0.5, label=None, plot=True, fig=None, Lref="a", refsp=None):

    q, t = get_Qi_t_nc(dirname)
    print(t)
    #species_type = data.groups['Inputs'].groups['Species'].variables['species_type'][ispec]
    #if species_type == 0:
    #    species_tag = "i"
    #elif species_type == 1:
    #    species_tag = "e"
    #if refsp == None:
    species_tag='i'

    # compute time-average and std dev
    istart_avg = int(len(t)*navgfac)
    qavg = np.mean(q[istart_avg:])
    qstd = np.std(q[istart_avg:])
    if label == None:
        label = dirname
    print(r"%s: Q_%s/Q_GB = %.5g +/- %.5g" % (label, species_tag, qavg, qstd))

    # make a Q vs time plot
    if plot:
        if fig == None:
            fig = plt.figure(0)
        plt.plot(t,q,'-',label=r"%s: $Q_%s/Q_\mathrm{GB}$ = %.5g"%(label, species_tag, qavg))
        plt.ylabel(r"$Q/Q_\mathrm{GB}$")
        plt.xlabel(r"$t\ (v_{t%s}/%s)$"%(refsp, Lref))
        legend = plt.legend(loc='upper right')
        legend.set_in_layout(False)
        plt.tight_layout()

if __name__ == "__main__":
    
    print("Plotting heat fluxes.....")
    ispec=0
    for fname in sys.argv[1:]:
        heat_flux(fname, ispec=ispec, refsp="i")
    
    plt.xlim(0)
    plt.ylim(0)
    # uncomment this line to save a PNG image of the plot
    #plt.savefig("heat_flux.png")
    plt.show()
