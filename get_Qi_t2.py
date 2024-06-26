#!/usr/bin/env python

import numpy as np
import numpy.ma as ma
from glob import glob
import os

from decide_version import decide_version
from netCDF4 import Dataset

def get_Qi_t_stdout(dirname,stdout_file='out.*'):
    stdout_file = sorted(glob(dirname + "/" + stdout_file))[-1]
    
    with open(stdout_file,'r') as f:
        lines = f.readlines()

    t = []
    Qi = []
    for l in lines:
        if "gx: Step" in l:
            ls = l.rsplit(':',1)[-1].split('=')
            for i,entry in enumerate(ls):
                if 'Time' in entry:
                    t.append(float(ls[i+1].split(',',1)[0]))
                if 'Q_i' in entry:
                    Qi.append(float(ls[i+1].split(None,1)[0]))
    return(np.array(Qi),np.array(t))


def get_dataset(ncfile,version):

    if not os.path.isfile(ncfile):
        return np.array([np.nan]), np.array([np.nan])
    
    with Dataset(ncfile,'r') as f:
        if version == 1:
            #f.replace('--', np.nan)
            # Q_gB = n_i T_i v_{Ti} \rho_{Ti}^2/a^2
            try:
                Qi = f["/Fluxes/qflux"][()][:,0]/(2*np.sqrt(2))
            except IndexError:
                Qi = np.array([np.nan])
                t = np.array([np.nan])
            t = np.sqrt(2)*f.variables["time"][()]
        else:
            try:
                Qi = f["/Diagnostics/HeatFlux_st"][()][:,0]/(2*np.sqrt(2))
            except IndexError:
                Qi = np.array([np.nan])
                t = np.array([np.nan])
            t = np.sqrt(2)*f.groups["Grids"].variables["time"][:]
    return Qi,t

def get_Qi_t_nc(dirname, version=None):

    if dirname[-3:] == ".nc":
        tmp = dirname.rsplit('/',1)
        dirname = tmp[0]
        _ncfile = tmp[1]
    else:
        _ncfile = None
    ncfile, version = decide_version(dirname)

    if _ncfile is not None:
        ncfile = dirname + "/" + _ncfile
    
    if ncfile is None:
        print("GX output does not exist in '" + dirname + "'")
        return (np.nan,np.nan)
    
    if isinstance(ncfile, list):
        Qis = []
        ts = []
        for f in ncfile:
            _Qi, _t = get_dataset(f,version)
            Qis.append(_Qi)
            ts.append(_t)
        Qi = np.concatenate(Qis)
        t = np.concatenate(ts)
    else:
        Qi, t = get_dataset(ncfile,version)
    
    # remove invalid entries
    # at the end of vector
    # netcdf4 returns a masked array, not normal array. See
    # https://numpy.org/doc/stable/reference/maskedarray.generic.html
    #print(t.mask)
    if ma.is_masked(t):
        while t.mask[-1] == True:
            t = t[:-1]
            Qi = Qi[:-1]
    ret = (np.array(Qi),np.array(t))  
    return ret



if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    fig, ax = plt.subplots()

    Ndirs = len(ds)
    cmap = plt.get_cmap("brg",Ndirs)
    legend = []
    
    for i,d in enumerate(ds):    
        Qi,t = get_Qi_t_nc(d)
        if not np.isnan(Qi).all():
            ax.plot(t,Qi,marker='.',color=cmap(i))
            legend.append(d)
    ax.set_xlabel(r'$t/[a/v_{\rm{Ti}}]$')
    ax.set_ylabel(r'$Q_{\rm{i}}/[Q_{\rm{gB}}]$')
    ax.legend(legend)
    plt.show()
