#!/usr/bin/env python


import numpy as np
from netcdf_util import netcdf_file
from new_kperp2 import get_kperp2


BIG_NUM = 1.0

def get_kperp2(dirname):
    gx_output = dirname + "/" + 'gx.nc'
    with netcdf_file(gx_output,'r',mmap=False) as f:
        # phi_vs_t(t, tube, zed, kx, ky, ri)
        phi2_vs_z = f['/Spectra/Phi2zt'][()][-1] # time, theta
        ky = f.variables['ky'][()]
        kx = f.variables['kx'][()]
        z = f.variables['zed'][()]
        kperp2 = f.variables['kperp2'][()]
        kperp2 = kperp2[:,0,0,:]
        gds2 = f['/Geometry/gds2'][()] # theta # |nabla y|^2
        gds22 = f['/Geometry/gds22'][()] # theta # |nabla x|^2 shat^2
        gds21 = f['/Geometry/gds21'][()] # theta # nabla x \cdot nabla y shat
        shat = f['/Geometry/shat'][()] # theta # x/q dq/dx = -x/iota d(iota)/dx
        
        jacob = f.variables['jacobian'][()] # nzed, nalpha    

        nablay2 = gds2
        nablax2 = gds22/shat**2
        nablax_nablay = gds21/shat
        

        phi2_vs_z = phi2_vs_z[:,0,:]
        kperp2 = (ky**2 * nablay2 + kx**2 * nablax2 + 2 * kx * ky * nablax_nablay) * phi2_vs_z * jacob/np.sum(phi2_vs_z * jacob,axis=0)
        kperp2 = np.sum(kperp2,axis=0)
        return kperp2



def get_gamma(dirname, tfit=20.0):

    gamma = np.load('my_growthrate.npy')
    ky = np.load('my_ky.npy')
    kx = np.load('my_kx.npy')
    return kx, ky, gamma


def rogerio_objective(dirname):
    try:
        kx, ky, gamma = get_gamma(dirname)
    except FileNotFoundError:
        print("failed to extract gamma from: " + dirname)
        return BIG_NUM
    kperp2 = get_kperp2(dirname)
    Nky = len(kperp2)
    return np.trapz(gamma/kperp2,ky)



if __name__=="__main__":
    # Print the value of the rogerio objective of a gx simulation.
    print("start")
    import sys
    argc = len(sys.argv)
    if argc > 1:
        dirname = sys.argv[1]
    else:
        dirname = '.'
    val = rogerio_objective(dirname)
    print("2222?")
    print(val)
