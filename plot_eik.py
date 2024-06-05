#!/usr/bin/env python

from eik_tools import EikFile
import sys
import matplotlib.pyplot as plt
import numpy as np


def plot_eik(files, outputs  = ["bmag", "gradpar", "gbdrift", "cvdrift", "cvdrift0", "gds2", "gds21", "gds22"], axes=None):
    
    nplots = len(outputs)
    if axes is None:
        if nplots <= 4:
            ncolsMax = 2
        else:
            ncolsMax = 4
        ncols = np.min((nplots,ncolsMax))
        nrows = (nplots-1)//ncols + 1
        fig,axes = plt.subplots(nrows=nrows,ncols=ncols, squeeze=False,sharex=True)
        axes = axes.flatten()
    else:
        axes = axes.flatten()
        assert len(axes) >= nplots

    
    for f in files:
        ef = EikFile(f)
        theta = ef.theta
        ntheta = len(theta)

        ys = np.zeros((8,ntheta))
        x = ef.scaled_theta/(np.pi)
        
        for i, output in enumerate(outputs):
            ys[i] = getattr(ef,output)
            axes[i].plot(x,ys[i])
            axes[i].set_title(output)

    return axes, ncols

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("eikfiles", nargs="+", metavar='eikfile', help='Simulation directories from which to take geometries.', default=['.'])
    parser.add_argument("--plot", nargs="+", action='store', metavar='variables', default=['bmag', 'gradpar', 'gbdrift', 'cvdrift', 'cvdrift0', 'gds2', 'gds21', 'gds22'], help='Can be used to specify exactly which variables (modB, gds2, etc) to plot.')
    parser.add_argument("-p","--preset", action='store', nargs='?', metavar='preset', default=0, help='Useful presets for plotting all variables that enter in the gyrokinetic equation, etc.', type=int)
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    parser.add_argument("-l","--legend", nargs="*", action='store', metavar='label',  help='Optional list of labels for legends', default = [])
    
    
    args = parser.parse_args()
    files = args.eikfiles
    preset = args.preset
    
    if preset == 0:
        add_variables = []
    elif preset == 1:
        add_variables = ['gradpar', 'gds2', 'gds21', 'gds22']
    elif preset == 2:
        add_variables = ['bmag', 'gradpar', 'gbdrift', 'cvdrift', 'cvdrift0', 'gds2', 'gds21', 'gds22']
    variables = args.plot
    variables = set(variables + add_variables)
    
    axes = None
    legends = []
    Nl = len(args.legend)
    for i, f in enumerate(files):
        if i < Nl:
            legends.append(args.legend[i])
        else:
            legends.append(None)
        

    axes, ncols = plot_eik(files, variables)
    axes[-1].legend(legends)

    for i in range(ncols):
        axes[-(i+1)].set_xlabel(r"$\theta/\pi$")


    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()
