#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

#import numpy as nmp
#from netCDF4 import Dataset
#import time ; # to check the execution speed...
#
import gonzag as gz
from gonzag.config import ldebug

# Should be command-line arguments:
#np_box_radius = 4 ; # in number of points... (should give km and get this one based on ocean model res...)
#dist_found = 100. ; # threshold distance for found [km] ! ""    ""
#l_dump_np_track_on_model_grid = True
#if_talk = 500 ; # verbose frequency: we chat every if_talk time steps !!!

def argument_parsing():
    '''
    ARGUMENT PARSING / USAGE
    '''
    import argparse as ap
    #
    parser = ap.ArgumentParser(description='Generate pixel maps of a given scalar.')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-s', '--fsat', required=True,                help='specify the NEMO netCDF file to read from...')
    requiredNamed.add_argument('-n', '--vsat', required=True, default="ssh", help='specify the field/diagnostic to plot (ex: SST)')
    requiredNamed.add_argument('-m', '--fmod', required=True,                help='specify the field/diagnostic to plot (ex: SST)')
    requiredNamed.add_argument('-v', '--vmod', required=True, default="ssh", help='specify the field/diagnostic to plot (ex: SST)')
    #
    parser.add_argument('-p', '--ewpr', type=int, default=-1, help='East-West periodicity of the model grid (-1 => None ; >=0 => periodicity with overlap of "ewper" points!))')
    parser.add_argument('-l', '--fmsk' , default="0",         help='specify model grid land-sea mask (<ncfile> or "0" for _FillValue)')
    parser.add_argument('-k', '--vmsk' , default="tmask",     help='variable name of model grid land-sea mask in')
    #
    args = parser.parse_args()
    print('')
    print(' *** Satellinte file and variable:\n', '  => ', args.fsat, '"'+args.vsat+'"')
    print(' ***    Model   file and variable:\n', '  => ', args.fmod, '"'+args.vmod,'"\n')
    #
    flsm = args.fmsk
    vlsm = args.vmsk
    ewpr = args.ewpr
    if flsm == "0": flsm = args.fmod ; vlsm = '_FillValue'
    #
    return args.fsat, args.vsat, args.fmod, args.vmod, flsm, vlsm, ewpr



if __name__ == '__main__':

    file_sat,  name_ssh_sat, file_mod, name_ssh_mod, file_lsm, name_lsm, kEWper = argument_parsing()






    # If Model2SatTrack were to be made independant of NetCDF interface, would provide it with:
    # * itime_s, itime_m
    # * xlat_m, xlon_m, mask_m
    # * vlat_s, vlon_s
    # * vssh_s[t]       => ok, small enough to be loaded in memory
    # * xssh_s[t,:,:]   => NOT ok, way too big!

    ii = gz.Model2SatTrack( file_sat, name_ssh_sat, file_mod, name_ssh_mod, file_lsm, name_lsm, \
                            ew_prd_mod=kEWper, file_out='test.nc' )

