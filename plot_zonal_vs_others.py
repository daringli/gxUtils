#!/usr/bin/env python

import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colormaps
from mpl_toolkits.axes_grid1 import make_axes_locatable

from netCDF4 import Dataset


from decide_version import decide_version


sqrt2 = np.sqrt(2)

def get_colors(colormap_name, n):
    return colormaps[colormap_name].resampled(n)(np.linspace(0.0, 1.0, n))


def get_Phiky(d, ispec=0, navgfac=0.5, plot=False, ax=None, Lref="a", refsp=None, colormap=None):
    if d[-3:] == ".nc":
        d = d.rsplit('/',1)[0]
    ncfile, version = decide_version(d)
    
    if isinstance(ncfile, list):
        # if we have restarted and got multiple output files
        # use the latest file
        ncfile = ncfile[-1]
    
    if version == 1:
        raise ValueError("non-next GX outputs not supported.")
        
    else:
        time_str = "Grids/time"
        ky_str = "Grids/ky"
        Phikyt_str = "/Diagnostics/Phi2_kyt"

    data = Dataset(ncfile, mode='r')
    t = sqrt2 * data[time_str][:]
    ky = sqrt2 * data[ky_str][:]
    #print(t)
    #print(version)
    try:
        Phikyt = data[Phikyt_str][:,:]
    except (KeyError, IndexError):
        print("error for '"  +ncfile +"', Skipping.")
        Phikyt = np.nan * np.zeros(len(ky), len(t))


    N = len(t)
    if colormap is not None:
        colors=get_colors(colormap, N)
    else:
        colors = [None] * N

    ax.set_prop_cycle('color', colors)
        
    if plot:
        if ax is None:
            fig, ax  =plt.subplots(1)
        try:
            _ky = np.array([ky] * len(t))
            ax.plot(_ky.T,Phikyt.T,marker='.', markersize=0.3)
        except:
            pass
        norm = mpl.colors.Normalize(vmin=t[0], vmax=t[-1])
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cb1 = mpl.colorbar.ColorbarBase(cax, cmap=colormap,
                                            norm=norm,
                                            orientation='vertical')

        #cb1.set_label(r'$k_y \rho$', labelpad=-40, y=1.05, rotation=0)
        cax.set_title(r'$t$')
        

    
    return ky, t, Phikyt

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("simdir", metavar='simdir', help='Simulation directory from which to read GX outputs from.', default=['.'])
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    parser.add_argument("-l","--legend", nargs="*", action='store', metavar='label',  help='Optional list of labels for legends', default = [])
    parser.add_argument("--colormap","--cm", choices=list(colormaps), metavar='cm',  help='Name of matplotlib colormap', default = 'rainbow')
    
    args = parser.parse_args()

    print("Plotting Phiky(t)")

    fig, ax = plt.subplots(1)

    refsp = 'i'
    ispec = 0
    Nl = len(args.legend)

    d = args.simdir
    label = ''
    
    get_Phiky(d, ax=ax, plot=True, refsp = refsp, ispec=ispec, colormap=args.colormap)
        
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r'$k_y \rho_{%s}$' % refsp)
    ax.set_ylabel(r"$\Phi$")
    plt.tight_layout()
    
    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()
