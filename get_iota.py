from netCDF4 import Dataset

def get_iota(d, fname = 'gx.nc'):
    
    
    data = Dataset(d + "/" + fname, mode='r')

    q = data.groups['Geometry'].variables['q'][()]
    return 1/q
