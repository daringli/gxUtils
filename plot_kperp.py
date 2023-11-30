#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from decide_version import decide_version
from netCDF4 import Dataset

sqrt2 = np.sqrt(2)

def get_kperp(dirname, against = 'ky'):
    
    ncfile, version = decide_version(dirname)
    
    with Dataset(ncfile,'r') as f:
        gds2 =  f["/Geometry/gds2"][()]
        gds21 = f["/Geometry/gds21"][()]
        gds22 = f["/Geometry/gds22"][()]
        bmag = f["/Geometry/bmag"][()]
        shat = f["/Geometry/shat"][()]

        if version == 2:
            ky = f["Grids/ky"]
            kx = f["Grids/kx"]
            theta = f["Grids/theta"]
        else:
            ky = sqrt2 * f["ky"]
            kx = sqrt2 * f["kx"]
            theta = f["theta"]

        gds2 = np.array(gds2)
        gds21 = np.array(gds21)
        gds22 = np.array(gds22)
        ky = np.array(ky)
        kx = np.array(kx)
        theta = np.array(theta)

        # theta, kx, ky
        gds2 = gds2[:,None,None]
        gds21 = gds21[:,None,None]
        gds22 = gds22[:,None,None]

        ky3d = ky[None, None, :]
        kx3d = kx[None, :, None]
        theta3d = theta[:, None, None]
        
        kperp2 = ky3d**2 * gds2 + kx3d**2 * gds22/shat**2 + 2* kx3d * ky3d * gds21/shat
    

        if against == "kxky":
            theta0 = 0.0
            i = np.argmin(np.abs(theta - theta0))
            plt.pcolor(ky3d[0] ,kx3d[0] ,kperp2[i])
            plt.xlabel(r"$\rho k_x$")
            plt.ylabel(r"$\rho k_y$")
            plt.title(r"$\rho^2 k_\perp^2 (\theta={:.2f})$".format(theta0))

            
        elif against == "thetaky":
            kx0 = 0.0
            i = np.argmin(np.abs(kx - kx0))
            x=theta3d[:,0,0]
            y=ky3d[0,0,:]
            print(x.shape)
            print(y.shape)
            print(kperp2[:,i].shape)
            plt.pcolor(x/np.pi, y  ,kperp2[:,i].T)
            plt.xlabel(r"$\theta/\pi$")
            plt.ylabel(r"$\rho k_y$")
            plt.title(r"$\rho^2 k_\perp^2 (k_x = {:.2f})$".format(kx0))
            
        elif against == "thetakx":
            ky0 = 0.35
            i = np.argmin(np.abs(ky - ky0))
            x=theta3d[:,0,0]
            y=kx3d[0,:,0]
            print(x.shape)
            print(y.shape)
            print(kperp2[:,:,i].shape)
            plt.pcolor(x/np.pi, y  ,kperp2[:,:,i].T)
            plt.xlabel(r"$\theta/\pi$")
            plt.ylabel(r"$\rho k_x$")
            plt.title(r"$\rho^2 k_\perp^2 (k_y = {:.2f})$".format(ky0))


        elif against == "ky":
            x=ky3d[0,0,:]
            y = np.sum(kperp2,axis =(0,1))
            plt.plot(x, y)
            plt.xlabel(r"$\rho k_y$")
            plt.ylabel(r"$\int k_\perp^2$")
            
        #plt.colorbar()
        
       

if __name__ == "__main__":
    
    import sys
    dirnames = sys.argv[1:]
    for dirname in dirnames:
        get_kperp(dirname)

    plt.show()
