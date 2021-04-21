#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

import sys
import numpy as nmp
from netCDF4 import Dataset, num2date, default_fillvals
from datetime import datetime, timedelta
from calendar import timegm
from geopy import distance

from gonzag.ncio import ToEpochTime

narg = len(sys.argv)
if narg != 2:
    print('Usage: '+sys.argv[0]+' <satellite_file>')
    sys.exit(0)
cf_sat  = sys.argv[1]


id_sat = Dataset(cf_sat)
#
clndr = id_sat.variables['time']
vtime = clndr[:]
vepoch = ToEpochTime( vtime, clndr.units, clndr.calendar )
#
vlat  = id_sat.variables['latitude'][:]
vlon  = id_sat.variables['longitude'][:]
#
id_sat.close()

Nt = len(vtime)

vdt_e = nmp.zeros(Nt-1); #, dtype=nmp.int64)
vdt_o = nmp.zeros(Nt-1)
vdr = nmp.zeros(Nt-1)

for jt in range(Nt-1):
    vdt_o[jt] = 86400.*vtime[jt+1] - 86400.*vtime[jt]
    vdt_e[jt] = vepoch[jt+1] - vepoch[jt]
    vdr[jt] = distance.distance( (vlat[jt+1],vlon[jt+1]) , (vlat[jt],vlon[jt])).km
    print(' jt, Dt [s], Dr [km] =', jt+1, round(vdt_e[jt],2),'s ', round(vdr[jt],2),'km ', round(vdt_o[jt],2),'s ' )


