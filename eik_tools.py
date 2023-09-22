#!/usr/bin/env python

from scipy.io import netcdf_file
from scipy.interpolate import interp1d
import numpy as np
#from netcdf_util import netcdf_file

class EikFile(object):

    extrapol_fudge = 1e-10
    attribs = ['bmag', 'gbdrift', 'gbdrift0', 'cvdrift', 'cvdrift0', 'gds2', 'gds21', 'gds22']
    
    def __init__(self, eikfile = 'eik.out'):
        self.eikfile = eikfile

        if self.eikfile[-2:] == 'nc':
            self.eiktype = 2
        else:
            self.eiktype = 1
        
        self._read_geometry()
            

    def _read_geometry_nc(self):
        with netcdf_file(self.eikfile) as f:
            self.theta =    f.variables['theta'][()]
            self.bmag =     f.variables['bmag'][()]
            self.gradpar =  f.variables['gradpar'][()]
            self.grho  =    f.variables['grho'][()]
            self.gbdrift =  f.variables['gbdrift'][()]
            self.gbdrift0 = f.variables['gbdrift0'][()]
            self.cvdrift =  f.variables['cvdrift'][()]
            self.cvdrift0 = f.variables['cvdrift0'][()]
            self.gds2 =     f.variables['gds2'][()]
            self.gds21 =    f.variables['gds21'][()]
            self.gds22 =    f.variables['gds22'][()]
            # scalars
            self.kxfac =    f.variables['kxfac'][()]
            self.q =        f.variables['q'][()]
            self.scale =    float(f.variables['scale'][()])
            if np.fabs(self.scale - np.round(self.scale)) < 1e-3:
                # probably meant to be an integer
                self.scale = int(np.round(self.scale))
            self.Rmaj =     f.variables['Rmaj'][()]
            self.shat =     f.variables['shat'][()]
            self.drhodpsi = f.variables['drhodpsi'][()]
            # metadata
            self.wout =     f.variables['wout'].filename

    def _read_geometry_eik(self):
        with open(self.eikfile, 'r') as f:
            lines = f.readlines()

        # this metadata is missing in the eikfile
        self.wout = None
        
        state = 0 # find header
        for l in lines:
            l = l.strip()
            ls = l.split()
            if len(ls) == 0:
                continue
            else:
                if state == 0:
                    if ls[0] == 'ntgrid':
                        state = 1 # next non-empty line is scalars
                elif state == 1:
                    # we assume the same order as always
                    # so the order is not set by the header
                    self.ntheta = int(ls[2])
                    self.drhodpsi = float(ls[3])
                    self.Rmaj = float(ls[4])
                    self.shat = float(ls[5])
                    self.kxfac = float(ls[6])
                    self.q = float(ls[7])
                    self.scale = float(ls[8])
                    if np.fabs(self.scale - np.round(self.scale)) < 1e-3:
                        # probably meant to be an integer
                        self.scale = int(np.round(self.scale))
                    print("scale" + str(self.scale))
                    
                    # we can init arrays since we know ntheta
                    self.gbdrift = np.zeros(self.ntheta)
                    self.gradpar = np.zeros(self.ntheta)
                    self.grho = np.zeros(self.ntheta)
                    self.theta = np.zeros(self.ntheta)
                    self.cvdrift = np.zeros(self.ntheta)
                    self.gds2 = np.zeros(self.ntheta)
                    self.bmag = np.zeros(self.ntheta)
                    self.gds21 = np.zeros(self.ntheta)
                    self.gds22 = np.zeros(self.ntheta)
                    self.cvdrift0 = np.zeros(self.ntheta)
                    self.gbdrift0 = np.zeros(self.ntheta)
                    
                    
                    state = 2 # search for next header
                elif state == 2:
                    if ls[0] == 'gbdrift':
                        state = 3 # next self.ntheta non-empty lines are arrays
                        i = 0
                        
                elif state == 3:
                    self.gbdrift[i] = float(ls[0])
                    self.gradpar[i] = float(ls[1])
                    self.grho[i] = float(ls[2])
                    self.theta[i] = float(ls[3])
                    i = i + 1
                    if i == self.ntheta:
                        state = 4 # search for next header
                        
                elif state == 4:
                    if ls[0] == 'cvdrift':
                        state = 5 # next self.ntheta non-empty lines are arrays
                        i = 0
                        
                elif state == 5:
                    self.cvdrift[i] = float(ls[0])
                    self.gds2[i] = float(ls[1])
                    self.bmag[i] = float(ls[2])
                    i = i + 1
                    if i == self.ntheta:
                        state = 6 # search for next header

                elif state == 6:
                    if ls[0] == 'gds21':
                        state = 7 # next self.ntheta non-empty lines are arrays
                        i = 0
                        
                elif state == 7:
                    self.gds21[i] = float(ls[0])
                    self.gds22[i] = float(ls[1])
                    i = i + 1
                    if i == self.ntheta:
                        state = 8

                elif state == 8:
                    if ls[0] == 'cvdrift0':
                        state = 9 # next self.ntheta non-empty lines are arrays
                        i = 0
                        
                elif state == 9:
                    self.cvdrift0[i] = float(ls[0])
                    self.gbdrift0[i] = float(ls[1])
                    i = i + 1
                    if i == self.ntheta:
                        break # we are done
        
        
    def _read_geometry(self):
        if self.eiktype == 2:
            self._read_geometry_nc()
        elif self.eiktype == 1:
            self._read_geometry_eik()
            
        self.ntheta = len(self.theta)

    def create_arclength_variable(self):
        # self.attrib_interpolator('b_dot_grad_z', uncoord=False)
        Ntheta = len(self.theta)
        if Ntheta % 2 == 0:
            raise NotImplementedException("Eik files always use an odd number of parallel grid points?")
        else:
            i0 = Ntheta//2 + 1 #location of middle
            Isort = np.argsort(self.theta)
            theta = self.theta[Isort]
            gradpar_theta = self.gradpar[Isort]
            L = np.zeros(Ntheta)
            dtheta = theta[1] - theta[0]
            L[(i0-1)::-1] = cumtrapz(1/gradpar_theta[i0::-1],theta[i0::-1],axis=0)
            L[(i0+1):] = cumtrapz(1/gradpar_theta[i0:],theta[i0:],axis=0)

            if self.theta[0] > 0:
                L = L[::-1]

        self.L = L
        
    @property
    def scaled_theta(self):
        return self.theta * self.scale

    # stella compatiblity
    @property
    def exb_nonlin(self):
        return self.kxfac/2.0

    @exb_nonlin.setter
    def exb_nonlin(self, val):
        self.kxfac = 2 * val

    
    
    # want a lowercase alias since
    # the netcdf and plain text file have different conventions
    @property
    def rmaj(self):
        return self.Rmaj

    @rmaj.setter
    def rmaj(self, val):
        self.Rmaj = val

    @property
    def z(self):
        return self.theta

    @z.setter
    def z(self, val):
        self.theta = val

    # FOR STELLA COMPATIBILITY
    # except not really, since zed refers to a different parallel coordinate
    @property
    def zed(self):
        return self.theta

    @zed.setter
    def zed(self, val):
        self.theta = val


    def _write_wout(self,f,variable):
        f.createVariable(variable, 'i', ())
        var = f.variables[variable]
        #var.assignValue()
        var.filename = self.wout
        
    def _write_scalar(self,f,variable):
        f.createVariable(variable, 'f', ())
        var = f.variables[variable]
        var.assignValue(getattr(self,variable))
        
        
    def _write_array(self,f,variable):
        f.createVariable(variable, 'f', ('z',))
        var = f.variables[variable]
        var[:] = getattr(self,variable)


    

    

    def avg_attrib(self, attrib):
        attr = getattr(self, attrib)
        # \int dl X/L - not a flux-surface average.
        return np.trapz(attr/self.gradpar, self.theta)/np.trapz(1/self.gradpar,self.theta)

    def avg_bmag(self):
        return self.avg_attrib('bmag')

    def avg_gbdrift(self):
        return self.avg_attrib('gbdrift')

    def avg_gbdrift0(self):
        return self.avg_attrib('gbdrift0')
    
    def avg_cvdrift(self):
        return self.avg_attrib('cvdrift')
    
    def avg_cvdrift0(self):
        return self.avg_attrib('cvdrift0')

    def avg_gds2(self):
        return self.avg_attrib('gds2')

    def avg_gds21(self):
        return self.avg_attrib('gds21')

    def avg_gds22(self):
        return self.avg_attrib('gds22')
    
    
    def attrib_interpolator(self, attrib, coord='L'):
        if coord=='L':
            # arclength coordinate
            assert(self.Nalpha==1)
            x = np.array(list(self.L))
        elif (coord=='zed' or coord=='z'or coord=='theta'):
            # always from -pi to pi
            x = np.array(list(self.theta))
        elif coord=='scaled_theta':
            # covers the range of the poloidal angle of the flux tube
            x = np.array(list(self.scaled_theta))
        
        # to allow for very slight extrapolation
        # so that we can copy from files with the same theta range
        x[0] =   x[0] * (1 + EikFile.extrapol_fudge)
        x[-1] = x[-1] * (1 + EikFile.extrapol_fudge)
        ret = interp1d(x, getattr(self,attrib))
        return ret
        
    def copy_attrib(self, copyFromEikFile, attrib, coord='L'):
        tmp = copyFromEikFile.attrib_interpolator(attrib, coord)
        if coord=='L':
            setattr(self,attrib,tmp(self.L))
        elif (coord=='zed' or coord=='z' or coord=='theta'):
            setattr(self,attrib,tmp(self.theta))
        elif coord=='scaled_theta':
            setattr(self,attrib,tmp(self.scaled_theta))

    def copy_bmag(self, copyFromEikFile, coord='L'):
        self.copy_attrib(copyFromEikFile, 'bmag', coord)

    def copy_gbdrift(self, copyFromEikFile, coord='L'):
        self.copy_attrib(copyFromEikFile, 'gbdrift', coord)

    def copy_gbdrift0(self, copyFromEikFile, coord='L'):
        self.copy_attrib(copyFromEikFile, 'gbdrift0', coord)

    def copy_cvdrift(self, copyFromEikFile, coord='L'):
        self.copy_attrib(copyFromEikFile, 'cvdrift', coord)

    def copy_cvdrift0(self, copyFromEikFile, coord='L'):
        self.copy_attrib(copyFromEikFile, 'cvdrift0', coord)

    def copy_gds2(self, copyFromEikFile, coord='L'):
        self.copy_attrib(copyFromEikFile, 'gds2', coord)


    def copy_gds21(self, copyFromEikFile, coord='L'):
        self.copy_attrib(copyFromEikFile, 'gds21', coord)

    def copy_gds22(self, copyFromEikFile, coord='L'):
        self.copy_attrib(copyFromEikFile, 'gds22', coord)

    def write_output(self,output_filename, netcdf = None):
        if netcdf is None:
            if self.eiktype == 2:
                netcdf = True
            else:
                netcdf = False
                
        if netcdf:
            # since the input may have been a text eik file, we cannot copy a potentially existing eik file.
            # we copy every dataset individually
            
            with netcdf_file(output_filename,'w', maskandscale=False) as f:
                f.createDimension('z', self.ntheta+1) #YEAH! +1

                #f.createVariable("theta",'f',('z',))
                #theta[:] = self.theta[:]
                #f.theta = theta

                self._write_array(f,'theta')
                self._write_array(f,'bmag')
                self._write_array(f,'gradpar')
                self._write_array(f,'grho')
                self._write_array(f,'gbdrift')
                self._write_array(f,'gbdrift0')
                self._write_array(f,'cvdrift')
                self._write_array(f,'cvdrift0')
                self._write_array(f,'gds2')
                self._write_array(f,'gds21')
                self._write_array(f,'gds22')

                self._write_scalar(f, 'kxfac')
                self._write_scalar(f, 'q')
                self._write_scalar(f, 'scale')
                self._write_scalar(f, 'Rmaj')
                self._write_scalar(f, 'shat')
                self._write_scalar(f, 'drhodpsi')


                self._write_wout(f, 'wout')
        else:
            with open(output_filename,'w') as f:
                f.write(str(self))
                
    def __str__(self):
        # based on gx/geometry_modules/vmec/src/geometric_coefficients.cpp
        def header_format(s):
            return "{:.6e}".format(s)

        def table_format(s):
            return "{:>20.10e}".format(s)

        
        ret = "ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q scale\n"
        ret = ret + str(self.ntheta//2) + " 1 " + str(self.ntheta) + " 1.0 " +  header_format(self.Rmaj) + " " + header_format(self.shat) + " 1.0 " + header_format(self.q) + " "  + str(self.scale) + " \n"
        
        ret = ret + "gbdrift\t gradpar\t grho\t tgrid\n"
        for i in range(self.ntheta):
            ret = ret + table_format(self.gbdrift[i]) + "\t" + table_format(self.gradpar[i]) + "\t" + table_format(self.grho[i]) + "\t" + table_format(self.theta[i]) + "\n"

        ret = ret + "cvdrift\t gds2\t bmag\t tgrid\n"
        for i in range(self.ntheta):
            ret = ret + table_format(self.cvdrift[i]) + "\t" + table_format(self.gds2[i]) + "\t" + table_format(self.bmag[i]) + "\t" + table_format(self.theta[i]) + "\n"

        ret = ret + "gds21\t gds22\t tgrid\n"
        for i in range(self.ntheta):
            ret = ret + table_format(self.gds21[i]) + "\t" + table_format(self.gds22[i]) + "\t" + table_format(self.theta[i]) + "\n"

        ret = ret + "cvdrift0\t gbdrift0\t tgrid\n"
        for i in range(self.ntheta):
            ret = ret + table_format(self.cvdrift0[i]) + "\t" + table_format(self.gbdrift0[i]) + "\t" + table_format(self.theta[i]) + "\n"

            
        return ret
    

    
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()
    
    parser.add_argument("eikfiles",nargs="+", metavar='eikfile')

    parser.add_argument("-o", "--output", action="store", help="Output files to write to. If none, will print to stdout.", dest='outputs', nargs="+", default =[])
    
    parser.add_argument("-z","--zero",action='store', nargs="+", choices=['bmag','gbdrift','gbdrift0','cvdrift','cvdrift0','gds2','gds21','gds22'], help="List of geometric coefficients to set to zero.", metavar='geometric coefficient',dest='zeros', default=[])

    parser.add_argument("-a","--avg",action='store', nargs="+", choices=['bmag','gbdrift','gbdrift0','cvdrift','cvdrift0','gds2','gds21','gds22'], help="List of geometric coefficients to set to its average.", metavar='geometric coefficient',dest='avgs', default=[])
    
    parser.add_argument("-t", "--text", action="store_true", default=False, help="Whether to write output as plain text file, as opposed to netcdf. Only used if --output is specified. Default False", dest='text')

    parser.add_argument("-c", "--coord", action='store', choices=['L','theta','scaled_theta'], help="Which parallel coordinate to use for interpolation when copying quantities from a different geo file.", dest='coord', default='theta')
    

    for attrib in EikFile.attribs:
        parser.add_argument("--" + attrib, action='store', help="An eikfile to copy " + attrib +" from.", metavar='eikfile',dest=attrib)
    
    args = parser.parse_args()

    nfilenames = len(args.eikfiles)
    noutputs = len(args.outputs)
    nzeros = len(args.zeros)
    fatal_exit = False
    
    if (noutputs > 0) and (noutputs != nfilenames):
        print("Output must be unspecified or match the number of inputs.")
        exit(1)

    override_eikfiles = {}
    for attrib in EikFile.attribs:
        inzero = attrib in args.zeros
        inavg = attrib in args.avgs
        copy_eikfile = getattr(args,attrib)
        override = copy_eikfile is not None
        if inzero and inavg and override:
            print("Cannot both set " + attrib + " to zero AND its average AND copy it from eikfile '" + copy_eikfile + "'.")
            exit(1)
        elif inzero and inavg:
            print("Cannot both set " + attrib + " to zero AND its average.")
            exit(1)
        elif  inzero and override:
            print("Cannot both copy " + attrib + " from eikfile '" + copy_eikfile + "' and set it to zero.")
            exit(1)
        elif  inavg and override:
            print("Cannot both copy " + attrib + " from eikfile '" + copy_eikfile + "' and set it to its average.")
            exit(1)

        if override:
            override_eikfiles[attrib] = EikFile(copy_eikfile)
            
            

    for i,filename in enumerate(args.eikfiles):
        eik = EikFile(filename)
        ntheta = len(eik.theta)
        z = np.zeros(ntheta)
        o = np.ones(ntheta)

        for attrib_to_zero in args.zeros:
            setattr(eik, attrib_to_zero, z)

        for attrib_to_avg in args.avgs:
            setattr(eik, attrib_to_avg, eik.avg_attrib(attrib_to_avg) * o)

        for attrib_to_override in override_eikfiles:
            eik.copy_attrib(override_eikfiles[attrib_to_override], attrib_to_override, coord=args.coord)
            
        if noutputs == 0:
            print(eik)
        else:
            eik.write_output(args.outputs[i], netcdf = not args.text)
    #eik.write_output('eik.nc')
