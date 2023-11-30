from decide_version import decide_version
from netCDF4 import Dataset


def get_dataset(ncfile,version):
    with Dataset(ncfile,'r') as f:
        if version == 1:
            #f.replace('--', np.nan)
            # Q_gB = n_i T_i v_{Ti} \rho_{Ti}^2/a^2
            gradpar = f["/Geometry/gradpar"][()]
        else:
            gradpar = f["/Geometry/gradpar"][()]
            
    return gradpar

def get_gradpar(dirname, version=None):

    ncfile, version = decide_version(dirname)
    if ncfile is None:
        print("GX output does not exist in '" + dirname + "'")
        return (np.nan,np.nan)
    gradpar = get_dataset(ncfile,version)
    return gradpar
