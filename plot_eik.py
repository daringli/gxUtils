#!/usr/bin/env python

from eik_tools import EikFile
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams.update({'font.size': 6, 'lines.linewidth' : 1})


def plot_file(f, label=None, titles=None):    
    ef = EikFile(f)
    theta = ef.theta
    ntheta = len(theta)

    ys = np.zeros((8,ntheta))
    x = ef.scaled_theta/(np.pi)
    outputs = ["bmag", "gradpar", "gbdrift", "cvdrift", "cvdrift0", "gds2", "gds21", "gds22"]

    if titles == None:
        titles = outputs
    
    if label == None:
        label = f
    
    for i, output in enumerate(outputs):
        ys[i] = getattr(ef,output)
        axes[i].plot(x,ys[i], label=label)
        axes[i].set_title(titles[i])
        axes[i].ticklabel_format(useOffset=False)
    

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("eikfiles", nargs="+", metavar='eikfile', help='Eikfiles from which to take geometries.', default=['.'])
    parser.add_argument("-o","--output", nargs="?", action='store', metavar='filename',  help='Optional filename to save output to.', default = None)
    parser.add_argument("-l","--legend", nargs="*", action='store', metavar='label',  help='Optional list of labels for legends', default = [])
    
    parser.add_argument("-p", "--pretty", action="store_true", default=False, help="Whether to write pretty LaTeX titles. Default False")

    parser.add_argument("-c", "--Ncols", nargs="?", action="store", default=2, help="Whether to write pretty LaTeX titles. Default False")
    

    args = parser.parse_args()

    ncols = int(args.Ncols)
    nrows = int(np.ceil(8/ncols))
    
    fig, axes = plt.subplots(nrows, ncols, sharex=True)
    axes = axes.flatten()

    
    Nl = len(args.legend)

    if args.pretty:
        titles = [r"$B$",r"$\nabla_{\|\|} \theta$", r"$\frac{2}{B^3}(\vec{B} \times \nabla B) \cdot \nabla y$", \
                  r"$\frac{2}{B^3}(\vec{B} \times \nabla B) \cdot \nabla y$", r"$\frac{2\hat{s}}{B^3}(\vec{B} \times \nabla B) \cdot \nabla x$", r"$|\nabla y|^2$", r"$\hat{s}\nabla y \cdot \nabla x$", r"$\hat{s}^2|\nabla x|^2$"]
    else:
        titles = None

    for i, f in enumerate(args.eikfiles):
        if i < Nl:
            label = args.legend[i]
        else:
            label = None
 
        plot_file(f, label=label, titles=titles)

    axes[1].legend()

    # some nice presets for common column numbers
    if ncols == 2:
        axes[-2].set_xlabel(r"$\theta/\pi$")
        axes[-1].set_xlabel(r"$\theta/\pi$")
        plt.subplots_adjust(left=0.1,
                        bottom=0.1, 
                        right=0.9, 
                        top=0.9, 
                        wspace=0.4, 
                        hspace=0.4)
    elif ncols == 4:
        axes[-4].set_xlabel(r"$\theta/\pi$")
        axes[-3].set_xlabel(r"$\theta/\pi$")
        axes[-2].set_xlabel(r"$\theta/\pi$")
        axes[-1].set_xlabel(r"$\theta/\pi$")
        plt.subplots_adjust(left=0.05,
                        bottom=0.3, 
                        right=0.98, 
                        top=0.9, 
                        wspace=0.5, 
                            hspace=0.25)

        
    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()
