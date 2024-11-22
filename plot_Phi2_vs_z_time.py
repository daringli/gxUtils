#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os

from netCDF4 import Dataset

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colormaps
from mpl_toolkits.axes_grid1 import make_axes_locatable


from decide_version import decide_version

from gx_input import Gx_input


sqrt2 = np.sqrt(2)

def get_colors(colormap_name, n):
    return colormaps[colormap_name].resampled(n)(np.linspace(0.0, 1.0, n))

class Phi2_data(object):
    def __init__(self, dirname, bignum = 1e20):
        """Read the netcdf gx output file in directory 'dirname'"""
        self.bignum = bignum
                
        
        ncfile, version = decide_version(dirname)
        if ncfile is None:                            
            raise ValueError("File not found")

        if isinstance(ncfile, list):
            # if we have restarted and got multiple output files
            # use the latest file
            # TODO: do better
            ncfile = ncfile[-1]


        
        if version == 1:
            raise ValueError("Non-next GX output is not supported")
        else:
            time_str = "Grids/time"
            theta_str = "Grids/theta"
            Phi2_str = "/Diagnostics/Phi2_zt"
            Phi2_zonal_str = "/Diagnostics/Phi2_zonal_zt"

        if os.path.isfile(ncfile):
            with Dataset(ncfile,'r') as f:
                self.t = f[time_str][()]
                self.theta = f[theta_str][()]
                self.Phi2 = f[Phi2_str][()]
                self.Phi2_zonal = f[Phi2_zonal_str][()]

                self.Nt = self.Phi2.shape[0]
                self.Ntheta = self.Phi2.shape[1]

                gi = Gx_input(dirname)
                if gi.geo_option == 'vmec':
                    npol = gi.npol
                else:
                    npol = 1
                self.theta = npol * self.theta
                
    def plot(self, ax=None, colormap = 'rainbow'):
        if ax is None:
            fig, ax  =plt.subplots(1)

        N = len(self.t)
        if colormap is not None:
            colors=get_colors(colormap, N)
        else:
            colors = [None] * N

        ax.set_prop_cycle('color', colors)

        try:
            _theta = np.array([self.theta] * N)
            ax.plot(_theta.T,self.Phi2.T,marker='.', markersize=0.3)
        except:
            pass
        norm = mpl.colors.Normalize(vmin=self.t[0], vmax=self.t[-1])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=colormap,
                                            norm=norm,
                                            orientation='vertical')

        #cb1.set_label(r'$k_y \rho$', labelpad=-40, y=1.05, rotation=0)
        cax.set_title(r'$t$')

    
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("simdir", metavar='simdir', help='Simulation directory from which to read GX outputs from.', default=['.'])
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    parser.add_argument("--colormap","--cm", choices=list(colormaps), metavar='cm',  help='Name of matplotlib colormap', default = 'rainbow')
    
    args = parser.parse_args()

    print("Plotting Phiky(t)")

    fig, ax = plt.subplots(1)

    refsp = 'i'
    ispec = 0
    
    d = args.simdir
    label = ''
    
    data =  Phi2_data(d)
    fig, ax = plt.subplots(1)
    data.plot(ax=ax)

    
    ax.set_xlabel(r'$k_y \rho_{%s}$' % refsp)
    ax.set_ylabel(r"$\Phi$")

    
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'$|\Phi|^2$')
    plt.tight_layout()

    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()
