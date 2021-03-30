#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

import gonzag as gz
from gonzag.config import ldebug

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
    #parser.add_argument('-p', '--ewpr', type=int, default=-1, help='East-West periodicity of the model grid (-1 => None ; >=0 => periodicity with overlap of "ewper" points!))')
    parser.add_argument('-l', '--fmsk' , default="0",         help='specify model grid land-sea mask (<ncfile> or "0" for _FillValue)')
    parser.add_argument('-k', '--vmsk' , default="tmask",     help='variable name of model grid land-sea mask in')
    #
    args = parser.parse_args()
    print('')
    print(' *** Satellite file and variable:\n', '  => ', args.fsat, '"'+args.vsat+'"')
    print(' ***   Model   file and variable:\n', '  => ', args.fmod, '"'+args.vmod,'"\n')
    #
    flsm = args.fmsk
    vlsm = args.vmsk
    #ewpr = args.ewpr
    if flsm == "0": flsm = args.fmod ; vlsm = '_FillValue'
    #
    return args.fsat, args.vsat, args.fmod, args.vmod, flsm, vlsm



if __name__ == '__main__':

    file_sat,  name_ssh_sat, file_mod, name_ssh_mod, file_lsm_mod, name_lsm_mod = argument_parsing()

    # Time overlap between model and satellite data ?
    (it1,it2), (Nts,Ntm) = gz.GetEpochTimeOverlap( file_sat , file_mod )
    print(' *** Time overlap between model and satellite in UNIX epoch time: it1, it2', it1,'--',it2)
    print('   => UTC: "'+gz.EpochT2Str(it1)+'" to "'+gz.EpochT2Str(it2)+'"\n')


    print('\n\n\n #####   M O D E L   2 D + T   D O M A I N   a.k.a.  S O U R C E   #####\n')

    clsm = name_lsm_mod
    if name_lsm_mod=='_FillValue': clsm = name_lsm_mod+'@'+name_ssh_mod
    
    ModelGrid = gz.ModGrid( file_mod, it1, it2, file_lsm_mod, clsm )
    
    
    print('\n\n\n #####   S A T E L L I T E   1 D   T R A C K   a.k.a.  T A R G E T   #####\n')

    SatelliteTrack = gz.SatTrack( file_sat, it1, it2, Np=Nts, domain_bounds=ModelGrid.domain_bounds, l_0_360=ModelGrid.l360 )


    # Mapping and interpolation:
    
    ii = gz.Model2SatTrack( ModelGrid, name_ssh_mod, SatelliteTrack, name_ssh_sat, file_out='result.nc' )

