import numpy as np
import toml

groupnames = ["qspectra", "species", "dimensions", 'domain', "physics", "time", "initialization", "geometry", "species", "boltzmann", "dissipation", "expert", "forcing"]


def val2str(val):
    if type(val) == bool:
        if val:
            return 'true'
        else:
            return 'false'
    else:
        return str(val)

def set_toml_variable(filename, group,var,value):
    with open(filename,'r') as f:
        lines = f.readlines()
    group = group.lower()
    startGroup = None
    endGroup = None
    step = 0
    searchstrings = ["["+group+"]"]
    #print(searchstring)
    # preprocess to extract boundary of group
    for i,l in enumerate(lines):
        l = l.split('#',1)[0] # strip comments
        l = l.lower()
        match = any([searchstring in l for searchstring in searchstrings])
        if match:
            if step == 0:
                startGroup = i
                searchstrings = ["[" + gn + "]" for gn in groupnames]
                step = 1
            elif step == 1:
                endGroup = i
                step = 2
                break

    #print(startGroup,endGroup)
    for i,l in enumerate(lines[startGroup:endGroup]):
        ls = l.split('#',1)
        l = ls[0]
        if var in l:
            if len(ls) > 1:
                lcomment = ls[1]
            else:
                lcomment = ''
            l = l.split('=')[0] + '= ' + val2str(value)
            #print(l + lcomment)
            lines[startGroup + i] = l + " # " + lcomment.strip() + "\n"

    with open(filename,'w') as f:
        f.write(''.join(lines))


class Gx_input(object):

    
    
    # Default values from gx/src/parameters.cu
    # some of the defaults are sentinel values with unclear meaning.
    # Currently does not contain everything, add more as needed.
    defaults = {}

    defaults["Dimensions"] = {
        "ntheta" : 32, \
        "ny" : 0, \
        "nx" : 0, \
        "nky" : 0, \
        "nkx" : 0, \
        "nhermite" : 4, \
        "nlaguerre" : 2, \
        "nspecies" : 1, \
        "nperiod" : 1, \
    }

    defaults["Domain"] = {
        "y0" : 10.0, \
        "x0" : -1.0, \
        "jtwist" : -1, \
        "zp" : -1.0, \
        "boundary" : "linked", \
        "long_wavelength_GK" : False, \
    }
        
    defaults["Time"] = {
        "dt" : 0.05, \
        "nstep" : 1e20, \
        "scheme" : "sspx3", \
        "cfl" : 0.9, \
        "stages" : 10, \
        "t_max" : 1e20, \
        "t_add" : -1.0, \
    }
    
    def __init__(self,dirname,input_name = "gx.in"):
        self.dirname = dirname
        self.input_name = input_name
        self.full_input_name = dirname + "/" + input_name
        with open(self.full_input_name, "r") as f:
            self.data = toml.load(f)

    def get_value_from_input_or_defaults(self,groupname,varname):
        varname = varname
        groupname = groupname
        #parser = f90nml.Parser()
        #parser.comment_tokens = "!"
        #inputs = parser.read(self.input_name)
        if not varname in self.data[groupname].keys():
            return Gx_input.defaults[groupname][varname]
        else:
            return self.data[groupname][varname]#.split('!',1)[0].strip()

    def get_group(self,var, species = None):
        """Return the toml group of the variable (if unambigious)"""
        ret = None
        for key in Gx_input.defaults:
            if var in Gx_input.defaults[key]:
                ret = key
                break
        if ret is None:
            raise ValueError("Cannot get group corresponding to variable '" + var + "'.")
        return ret

    def changevar(self, group, var, val):
        set_toml_variable(self.full_input_name,group,var,val)
        with open(self.full_input_name, "r") as f:
            self.data = toml.load(f)



    # below are setters and getters for the individual input parameters
    @property
    def nky(self):
        return self.get_value_from_input_or_defaults("Dimensions","nky")

    @nky.setter
    def nky(self,val):
        self.changevar("Dimensions","nky",val)

        
    @property
    def nkx(self):
        return self.get_value_from_input_or_defaults("Dimensions","nkx")

    @nkx.setter
    def nkx(self,val):
        self.changevar("Dimensions","nkx",val)


        
    @property
    def nx(self):
        return self.get_value_from_input_or_defaults("Dimensions","nx")

    @nx.setter
    def nx(self,val):
        self.changevar("Dimensions","nx",val)

            
    @property
    def ny(self):
        return self.get_value_from_input_or_defaults("Dimensions","ny")

    @ny.setter
    def ny(self,val):
        self.changevar("Dimensions","ny",val)


    

        
    @property
    def dt(self):
        return self.get_value_from_input_or_defaults("Time","dt")

    @dt.setter
    def dt(self,val):
        self.changevar("Time","dt",val)

    @property
    def alpha(self):
        return self.get_value_from_input_or_defaults("Geometry","alpha")

    @alpha.setter
    def alpha(self,val):
        self.changevar("Geometry","alpha",val)

    @property
    def alpha0(self):
        return self.alpha # for stella compatibility

    @alpha0.setter
    def alpha0(self,val): # for stella compatibility
        self.alpha = val

    @property
    def desired_normalized_toroidal_flux(self):
        return self.get_value_from_input_or_defaults("Geometry","desired_normalized_toroidal_flux")
    
    @desired_normalized_toroidal_flux.setter
    def desired_normalized_toroidal_flux(self,val):
        self.changevar("Geometry","desired_normalized_toroidal_flux",val)

    @property
    def torflux(self): #stella compatibility
        return self.desired_normalized_toroidal_flux

        
    @torflux.setter
    def torflux(self,val): #stella compatibility
        self.desired_normalized_toroidal_flux = val
    

    @property
    def ntheta(self):
        return self.get_value_from_input_or_defaults("Dimensions","ntheta")

    @ntheta.setter
    def ntheta(self,val):
        self.changevar("Dimensions","ntheta",val)

    
    @property
    def npol(self):
        return self.get_value_from_input_or_defaults("Geometry","npol")

    @npol.setter
    def npol(self,val):
        self.changevar("Geometry","npol",val)


    @property
    def y0(self):
        return self.get_value_from_input_or_defaults("Domain","y0")

    @y0.setter
    def y0(self,val):
        self.changevar("Domain","y0",val)

    
    @property
    def x0(self):
        return self.get_value_from_input_or_defaults("Domain","x0")

    @x0.setter
    def x0(self,val):
        self.changevar("Domain","x0",val)



    @property
    def vmec_file(self):
        return self.get_value_from_input_or_defaults("Geometry","vmec_file")
    
    @vmec_file.setter
    def vmec_file(self,val):
        self.changevar("Geometry","vmec_file",val)

    @property 
    def vmec_filename(self): #stella compatibility
        return self.vmec_file
    
    @vmec_filename.setter
    def vmec_filename(self,val): # stella compatibility
        self.vmec_file = val
    
    @property
    def nspecies(self):
        return self.get_value_from_input_or_defaults("Dimensions","nspecies")
        
    @nspecies.setter
    def nspecies(self,val):
        self.changevar("Dimensions","nspecies",val)
        
    @property
    def tprim(self):
        return np.array(self.get_value_from_input_or_defaults("species","tprim"))
    
    @tprim.setter
    def tprim(self, val):
        val = list(val)
        #print(val)
        if len(val) != self.nspecies:
            print("Warning, setting tprim array to an array that doesn't have nspecies entries")
        self.changevar("species","tprim",val)

    @property
    def fprim(self):
        return np.array(self.get_value_from_input_or_defaults("species","fprim"))
    
    @fprim.setter
    def fprim(self, val):
        val = list(val)
        #print(val)
        if len(val) != self.nspecies:
            print("Warning, setting fprim array to an array that doesn't have nspecies entries")
        self.changevar("species","fprim",val)


    @property
    def z(self):
        return np.array(self.get_value_from_input_or_defaults("species","z"))
    
    @z.setter
    def z(self, val):
        val = list(val)
        #print(val)
        if len(val) != self.nspecies:
            print("Warning, setting z array to an array that doesn't have nspecies entries")
        self.changevar("species","z",val)

    
    @property
    def mass(self):
        return np.array(self.get_value_from_input_or_defaults("species","mass"))
    
    @mass.setter
    def mass(self, val):
        val = list(val)
        #print(val)
        if len(val) != self.nspecies:
            print("Warning, setting mass array to an array that doesn't have nspecies entries")
        self.changevar("species","mass",val)


    
    @property
    def Qspectra_kx(self):
        return self.get_value_from_input_or_defaults("Qspectra","kx")
    
    @Qspectra_kx.setter
    def Qspectra_kx(self, val):
        self.changevar("Qspectra","kx",val)
    
    @property
    def Qspectra_ky(self):
        return self.get_value_from_input_or_defaults("Qspectra","ky")
    
    @Qspectra_ky.setter
    def Qspectra_ky(self, val):
        self.changevar("Qspectra","ky",val)

    @property
    def Qspectra_z(self):
        return self.get_value_from_input_or_defaults("Qspectra","z")
    
    @Qspectra_z.setter
    def Qspectra_z(self, val):
        self.changevar("Qspectra","z",val)
        
    @property
    def Qspectra_kxky(self):
        return self.get_value_from_input_or_defaults("Qspectra","kxky")
    
    @Qspectra_kxky.setter
    def Qspectra_kxky(self, val):
        self.changevar("Qspectra","kxky",val)


    @property
    def all_zonal_scalars(self):
        return self.get_value_from_input_or_defaults("Diagnostics","all_zonal_scalars")
    
    @all_zonal_scalars.setter
    def all_zonal_scalars(self, val):
        self.changevar("Diagnostics","all_zonal_scalars",val)
    
    @property
    def all_zonal(self):
        return self.get_value_from_input_or_defaults("Diagnostics","all_zonal")
    
    @all_zonal.setter
    def all_zonal(self, val):
        self.changevar("Diagnostics","all_zonal",val)

       
    def ny(self):
        return self.get_value_from_input_or_defaults("Dimensions","ny")

    @ny.setter
    def ny(self,val):
        self.changevar("Dimensions","ny",val)

    
    @property
    def nx(self):
        return self.get_value_from_input_or_defaults("Dimensions","nx")

    @nx.setter
    def nx(self,val):
        self.changevar("Dimensions","nx",val)

    
    @property
    def x0(self):
        return self.get_value_from_input_or_defaults("Domain","x0")

    @x0.setter
    def x0(self,val):
        self.changevar("Domain","x0",val)

        
    @property
    def y0(self):
        return self.get_value_from_input_or_defaults("Domain","y0")

    @y0.setter
    def y0(self,val):
        self.changevar("Domain","y0",val)



if __name__=="__main__":

    import sys

    if len(sys.argv) > 1:
        inputname = sys.argv[1]
        tmp = inputname.rsplit('/',1)
        dirname = tmp[0]
        inputname = tmp[1]
    else:
        inputname = 'gx.in'
        dirname = '.'
    print(dirname)
    gi = Gx_input(dirname,inputname,debug=True)
    print(gi.nky)
    gi.nky = 1
    
