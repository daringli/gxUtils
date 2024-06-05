#!/usr/bin/env python

import os
from glob import glob

def decide_version(dirname):

    ncfile1 = dirname + "/gx.nc"
    ncfile2 = dirname + "/gx.out.nc"

    if os.path.isfile(ncfile1):
        version = 1
        ncfile = ncfile1
        
    elif os.path.isfile(ncfile2):
        version = 2
        ncfile = ncfile2
    else:
        ncfile = None
        version = None

    # check for multiple files
    ncfiles = sorted(glob(dirname + "/gx.out[0-9].nc"))
    if len(ncfiles) > 0:
        version = 2
        ncfile = ncfiles
    
    return ncfile, version
