import numpy as np
import toml

def set_toml_variable(filename, group,var,value):
    with open(filename,'r') as f:
        lines = f.readlines()
    startGroup = None
    endGroup = None
    step = 0
    searchstring = "["+group+"]"
    print(searchstring)
    # preprocess to extract boundary of group
    for i,l in enumerate(lines):
        l = l.split('#',1)[0] # strip comments
        if searchstring in l:
            if step == 0:
                startGroup = i
                searchstring = "["
                step = 1
            elif step == 1:
                endGroup = i
                step = 2
                break

    print(startGroup,endGroup)
    for i,l in enumerate(lines[startGroup:endGroup]): 
        if var in l:
            ls = l.split('#',1)
            if len(ls) > 1:
                lcomment = ls[1]
            else:
                lcomment = ''
            l = ls[0]
            l = l.split('=')[0] + '= ' + str(value)
            print(l + lcomment)
            lines[startGroup + i] = l + " #" + lcomment + "\n"

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
    def dt(self):
        return self.get_value_from_input_or_defaults("Time","dt")

    @dt.setter
    def dt(self,val):
        self.changevar("Time","dt",val)
    
    
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
    
