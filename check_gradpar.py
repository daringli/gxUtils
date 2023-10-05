from scipy.io import netcdf_file
import numpy as np

if __name__ == "__main__":
    import sys
    ncfile = sys.argv[1]


    with netcdf_file(ncfile) as f:
        gradpar = f.variables["gradpar"][()]


    np.set_printoptions(32)
    print(gradpar - np.mean(gradpar))
