#!/usr/bin/env python

import numpy as np
from glob import glob

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


def get_Qi_t_nc(dirname,ncfile='gx.nc'):
    ncfile = dirname + "/" + ncfile
    
    with Dataset(ncfile,'r') as f:
        print(f)
        Qi = f["/Fluxes/qflux"][()]
        print(Qi)
        t = f.variables["time"][()]
    return(np.array(Qi),np.array(t))



if __name__ == "__main__":
    import sys
    import matplotlib.pyplot as plt
    if len(sys.argv) > 1:
        ds = sys.argv[1:]
    else:
        ds = ['.']

    fig, ax = plt.subplots()
    
    for d in ds:    
        Qi,t = get_Qi_t_nc(d)
        ax.plot(t,Qi,marker='.')
    ax.set_xlabel(r'$t/[a/v_{\rm{Ti}}]$')
    ax.set_ylabel(r'$Q_{\rm{i}}/[Q_{\rm{gB}}]$')
    ax.legend(ds)
    plt.show()
