#!/usr/bin/env python

from eik_tools import EikFile
import sys
import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(4,2, sharex=True)
axes = axes.flatten()

files = sys.argv[1:]
for f in files:
    
    ef = EikFile(f)
    theta = ef.theta
    ntheta = len(theta)
    
    ys = np.zeros((8,ntheta))
    x = ef.scaled_theta/(np.pi)
    outputs = ["bmag", "gradpar", "gbdrift", "cvdrift", "cvdrift0", "gds2", "gds21", "gds22"]

    for i, output in enumerate(outputs):
        ys[i] = getattr(ef,output)
        axes[i].plot(x,ys[i])
        axes[i].set_title(output)

axes[1].legend(files)

axes[-2].set_xlabel(r"$\theta/\pi$")
axes[-1].set_xlabel(r"$\theta/\pi$")

plt.show()
