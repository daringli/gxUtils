#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

from simsopt.mhd.vmec_diagnostics import vmec_fieldlines
from simsopt.mhd import Vmec

import sys
#from gx_input import Gx_input as GK_input
from stella_input import Stella_input as GK_input


from scipy.integrate import cumtrapz


def plot_GK_input(dirname, axes=None, n=200, coord='L', variables=[], preset = 2, normalization='code', npol=3):

    d = dirname

    if preset == 1:
        add_variables = ['B_sup_phi', 'gds2', 'gds21', 'gds22']
    elif preset == 2:
        add_variables = ['modB', 'B_sup_phi', 'gbdrift', 'cvdrift', 'cvdrift0', 'gds2', 'gds21', 'gds22']
    elif preset == 3:
        add_variables = ['modB', 'B_sup_theta_pest', 'B_sup_phi', 'B_cross_grad_B_dot_grad_alpha', 'B_cross_grad_B_dot_grad_psi',
                     'B_cross_kappa_dot_grad_alpha', 'B_cross_kappa_dot_grad_psi',
                     'grad_alpha_dot_grad_alpha', 'grad_alpha_dot_grad_psi', 'grad_psi_dot_grad_psi',
                     'bmag', 'gradpar_theta_pest', 'gradpar_phi', 'gbdrift', 'gbdrift0', 'cvdrift', 'cvdrift0', 'gds2', 'gds21', 'gds22']

    else:
        add_variables = []
    variables.extend(add_variables)
    variables = list(set(variables))
                  
    #variables = ['B_cross_kappa_dot_grad_alpha', 'grad_alpha_dot_grad_alpha', 'modB', 'grad_psi_dot_grad_psi', 'gds2']
    #variables = ['B_cross_kappa_dot_grad_alpha', 'grad_alpha_dot_grad_alpha', 'modB']
    #variables = ['grad_alpha_dot_grad_alpha']
    #variables = ['gradpar_phi']

    nplots = len(variables)
    if nplots <= 4:
        ncolsMax = 2
    else:
        ncolsMax = 4
    ncols = np.min((nplots,ncolsMax))

    nrows = (nplots-1)//ncols + 1

    if axes is None:
        fig,axes = plt.subplots(nrows=nrows,ncols=ncols, squeeze=False,sharex=True)
        axes = axes.flatten()

   
    gi = GK_input(d)
    wout_filename = d +"/" + gi.vmec_filename

    vmec = Vmec(wout_filename, mpi=None)
    s = gi.torflux
    alpha0 = gi.alpha0

    if n % 2 == 0:
        ntheta = n + 1
    else:
        ntheta = n
    print(ntheta)
    print(npol)
    theta1d =  np.linspace(-npol * np.pi, npol * np.pi, num=ntheta)

    f = vmec_fieldlines(vmec, s, alpha0, theta1d=theta1d, phi_center=0, plot=False, show=False)

    phi = f.phi[0,0]
    Nphi = len(phi)
    
    i0 = Nphi//2 + 1
    Isort = np.argsort(phi)
    phi = phi[Isort]
    gradpar_phi = f.gradpar_phi[0,0,Isort]
    L = np.zeros(Nphi)
    dphi = phi[1] - phi[0]
    L[(i0-1)::-1] = cumtrapz(1/gradpar_phi[i0::-1],phi[i0::-1])
    L[(i0+1):] = cumtrapz(1/gradpar_phi[i0:],phi[i0:])

    if f.phi[0,0,0] > 0:
        L = L[::-1]

    #print(L)
    #print(L[Nphi//2])

    istart = np.argmin(np.fabs(L + 1.25))
    istop = np.argmin(np.fabs(L - 1.25))

    #np.sum(f.B_cross_kappa_dot_grad_alpha[istart:istop])/L[istart:istop]
    if istart > istop:
        istart,istop = istop,istart
    avg_curv = np.mean(f.B_cross_kappa_dot_grad_alpha[0,0,istart:istop])
    print(d)
    print(-avg_curv)


    for j, variable in enumerate(variables):
        irow = j//ncols
        if variable == 'B_cross_kappa_dot_grad_alpha':
            setattr(f,variable, -eval("f." + variable))
        if coord=='phi':
            axes[j].plot(f.phi[0, 0, :], eval("f." + variable + '[0, 0, :]'))
            if irow == nrows-1:
                axes[j].set_xlabel('Cylindrical angle $\phi$')
        elif coord =='L':
            axes[j].plot(L, eval("f." + variable + '[0, 0, :]'))
            if irow == nrows-1:
                axes[j].set_xlabel('Arclenght $l(\phi)$')
        elif coord =='theta' or coord=='theta_pest':
            axes[j].plot(f.theta_pest[0,0,:], eval("f." + variable + '[0, 0, :]'))
            if irow == nrows-1:
                axes[j].set_xlabel(r'Pest angle $\theta$')
        axes[j].set_title(variable)

    #plt.figtext(0.5, 0.995, f's={f.s[0]}, alpha={f.alpha[0]}', ha='center', va='top')

    #plt.tight_layout()
    #axes[-1].legend(dirnames)   
    #plt.show()
    return axes

    
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("simdirs", nargs="+", metavar='simdir', help='Simulation directories from which to take geometries.', default=['.'])
    parser.add_argument("--plot", nargs="+", action='store', metavar='variables', default=['modB', 'B_sup_phi', 'gbdrift', 'cvdrift', 'cvdrift0', 'gds2', 'gds21', 'gds22'], help='Can be used to specify exactly which variables (modB, gds2, etc) to plot.')
    parser.add_argument("-p","--preset", action='store', nargs='?', metavar='preset', default=1, help='Useful presets for plotting all variables that enter in the gyrokinetic equation, etc.', type=int)
    parser.add_argument("-n","--nz", nargs="+", action='store', metavar='n', default=999, help='List or integer. Parallel resolution to plot with.', dest ='n', type=int)
    parser.add_argument("--normalization", nargs="?", action='store', metavar='normalization', default='code', choices=["code","si"], help='Whether to use "code" normalization and names (gds2) or SI units (|nabla alpha|^2).')
    parser.add_argument("-c","--coord", nargs="?", action='store', metavar='preset', default="L", help='Parallel coordinate to plot against.', choices = ['L','phi','theta'])
    parser.add_argument("-r","--resolution", action='store_true', dest='resol', help='Whether to plot a resolution scan.')
    parser.add_argument("--npol", "--np", action='store', dest='npol', help='Number of poloidal turns to plot.', default=6, type=int)
    
    args = parser.parse_args()
    if args.resol:
        ns = []
        dirnames = []
        for d in args.simdirs:
            for n in args.n:
                print(n)
                dirnames.append(d)
                ns.append(n)

    else:
        dirnames = args.simdirs
        ns = args.n
    if type(ns) == int:
        ns = [ns]
    if type(dirnames) == str:
        dirnames = [dirnames]


    coord = args.coord
    normalization = args.normalization
    preset = args.preset

    
    axes = None
    legends = []
    
    for d, n in zip(dirnames, ns):
        legends.append(d + " ntheta=" + str(n))
        axes = plot_GK_input(d, axes, n=n, coord=coord, preset = preset, normalization=normalization, npol=args.npol)
    axes[-1].legend(legends)   
    plt.show()
