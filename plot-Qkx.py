#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import sys

from netCDF4 import Dataset

sqrt2 = np.sqrt(2)

plt.figure(0)
for fname in sys.argv[1:]:
  if fname[-3:] != ".nc":
    # default filename for outputs
    fname = fname + "/gx.nc"
  
  data = Dataset(fname, mode='r')
  t = sqrt2 * data.variables['time'][:]
  kx = sqrt2 * data.variables['kx'][:]
  Qkxt = data.groups['Spectra'].variables['Qkxst'][:,0,:]/(2*sqrt2)
  Qkx = np.mean(Qkxt[int(len(t)/2):], axis=0)
  plt.plot(kx, Qkx, 'o-')
  print(np.sum(Qkx))

refsp = 'i'
plt.xlabel(r'$k_y \rho_{%s}$' % refsp)
plt.ylabel(r"$Q/Q_\mathrm{GB}$")
#plt.xscale('log')
plt.yscale('log')
plt.xscale('log')
plt.tight_layout()
plt.legend(sys.argv[1:])
plt.show()
