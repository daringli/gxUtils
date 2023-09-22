#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.stats as stats
import sys

from netCDF4 import Dataset

sqrt2 = np.sqrt(2)
fig, ax = plt.subplots(1)

for fname in sys.argv[1:]:
  if fname[-3:] != ".nc":
    # default filename for outputs
    fname = fname + "/gx.nc"
  
  data = Dataset(fname, mode='r')
  t = sqrt2 * data.variables['time'][:]
  ky = sqrt2 * data.variables['ky'][:]
  kx = sqrt2 * data.variables['kx'][:]
  Qkxkyt = data.groups['Spectra'].variables['Qkxkyst'][:,0,:]/(2*sqrt2)
  Qkxky = np.mean(Qkxkyt[int(len(t)/2):], axis=0)

  vmax = Qkxky.max()
  vmin = np.max((Qkxky.min(),vmax * 1e-3))
  print(vmin,vmax)
  pcm = ax.pcolormesh(kx,ky, Qkxky, norm=colors.LogNorm(vmin=vmin, vmax=vmax),
                   cmap='magma', shading='auto')
  print(np.sum(Qkxky))

refsp = 'i'
plt.xlabel(r'$k_x \rho_{%s}$' % refsp)
plt.ylabel(r'$k_y \rho_{%s}$' % refsp)
#plt.ylabel(r"$Q/Q_\mathrm{GB}$")
#plt.xscale('log')

cb = fig.colorbar(pcm, ax=ax, extend='max')
cb.set_label(r"$Q/Q_\mathrm{GB}$", rotation=0)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim([0,np.max(ky)])
ax.set_xlim([0,np.max(kx)])

plt.tight_layout()
#plt.legend(sys.argv[1:])
plt.show()
