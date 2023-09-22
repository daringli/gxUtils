#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os



from netCDF4 import Dataset

def test_Qi_ky_nc(dirname,ncfile='gx.nc', tfit = 85):
    ncfile = dirname + "/" + ncfile
    if os.path.isfile(ncfile):
        with Dataset(ncfile,'r') as f:
            #f.replace('--', np.nan)
            # Qkxkyst(time, s, ky, kx)
            Qikxky = f["/Spectra/Qkxkyst"][()][:,0,:,:]/(2*np.sqrt(2)) # time, species, ky, kx
            Qiky = f["/Spectra/Qkyst"][()][:,0,:]/(2*np.sqrt(2)) # time, species, ky
            Qizst = f["/Spectra/Qzst"][()][:,0,:]/(2*np.sqrt(2)) # time, s, theta
            Qi = f["/Fluxes/qflux"][()][:,0]/(2*np.sqrt(2))
            t = np.sqrt(2)*f.variables["time"][()]
            ky = np.sqrt(2)*f.variables["ky"][()]
            kx = np.sqrt(2)*f.variables["kx"][()]
        # remove invalid entries
        # at the end of vector
        # netcdf4 returns a masked array, not normal array. See
        # https://numpy.org/doc/stable/reference/maskedarray.generic.html
        #print(t.mask)
        if ma.is_masked(t):
            while t.mask[-1] == True:
                t = t[:-1]
                Qikxky = Qikxky[:-1]

    istart = np.argmin(np.fabs(t - tfit))
    T = t[-1] - t[istart]
    Qikxky = np.trapz(Qikxky[istart:],t[istart:], axis=0)/T # ky,kx
    Qiky = np.trapz(Qiky[istart:],t[istart:], axis=0)/T # ky
    Qiky2 = np.sum(Qikxky,axis=1)

    Qi = np.trapz(Qi[istart:],t[istart:], axis=0)/T

    Qizst = np.trapz(Qizst[istart:],t[istart:], axis=0)/T # z
    

    print("Qi = " + str(Qi))
    Qi2 = np.sum(Qiky)
    print("Qi2 = " + str(Qi2))
    Qi3 = np.sum(Qizst)
    print("Qi3 = " + str(Qi3)) # jacobian is already baked into Qizst due to choice of theta coordinates

    Nky = len(ky)
    avg_dky = (ky[-1] - ky[0])/(Nky-1)
    Qi2_trapz = np.trapz(Qiky,ky)/avg_dky + (Qiky[0] + Qiky[-1])/2
    print("Qi2_trapz = " + str(Qi2_trapz))

    Qiky_NU = list(Qiky[:])
    ky_NU = list(ky[:])

    del Qiky_NU[13]
    del ky_NU[13]
    del Qiky_NU[14]
    del ky_NU[14]
    
    Nky_NU = len(ky_NU)
    avg_dky_NU = (ky_NU[-1] - ky_NU[0])/(Nky_NU-1)
    
    Qi2_trapz_NU = np.trapz(Qiky_NU,ky_NU)/avg_dky_NU + (Qiky_NU[0] + Qiky_NU[-1])/2
    print("Qi2_trapz_NU = " + str(Qi2_trapz_NU))
    
    return Qiky, Qiky2, ky

if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    fig, ax = plt.subplots()

    Ndirs = len(ds)
    cmap = plt.get_cmap("brg",2*Ndirs)
    legend = []
    
    for i,d in enumerate(ds):    
        Qi,Qi2,ky = test_Qi_ky_nc(d)
        if not np.isnan(Qi).all():
            ax.plot(ky,Qi,marker='.',color=cmap(2*i))
            ax.plot(ky,Qi2,marker='.',color=cmap(2*i+1),linestyle='dashed',label='_nolegend_')
            legend.append(d)
    ax.set_xlabel(r'$k_y \rho$')
    ax.set_ylabel(r'$Q_{\rm{i}}/[Q_{\rm{gB}}]$')
    ax.legend(legend)
    plt.show()
