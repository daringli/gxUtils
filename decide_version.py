#!/usr/bin/env python

import os

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
    
    return ncfile, version
