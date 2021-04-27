#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

import gonzag as gz
from gonzag.config import IsZarr,ldebug, rmissval, l_save_track_on_model_grid


def argument_parsing():
    '''
    ARGUMENT PARSING / USAGE
    '''
    import argparse as ap
    #
    parser = ap.ArgumentParser(description='I interpolate, in space and time, the SSH (or any other 2D field) from a non-regular structured gridded domain of an OGCM onto a satellite track.')
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument('-s', '--fsat', required=True,                help='satellite alongtrack data NetCDF file to read from')
    requiredNamed.add_argument('-n', '--vsat', required=True, default="ssh", help='name of field of interest in satellite file (default="ssh")')
    requiredNamed.add_argument('-m', '--fmod', required=True,                help='model NetCDF file to read from')
    requiredNamed.add_argument('-v', '--vmod', required=True, default="ssh", help='name of field of interest in model file (default="ssh")')
    #
    parser.add_argument('-l', '--fmsk' , default="0",         help='land-sea mask on model grid: "-l <ncfile>" or "-l 0" to use "_FillValue" attribute of field of interest')
    parser.add_argument('-k', '--vmsk' , default="tmask",     help='name of land-sea mask in <ncfile> if "-l <ncfile>" used')
    #
    parser.add_argument('-D','--distgrid', dest='distgrid', action='store_true', help='take into acount possible strong local distorsion of the model grid')
    parser.set_defaults(distgrid=False)
    #
    args = parser.parse_args()
    print('')
    print(' *** Satellite file and variable:\n', '  => ', args.fsat, '"'+args.vsat+'"')
    print(' ***   Model   file and variable:\n', '  => ', args.fmod, '"'+args.vmod,'"\n')
    #
    flsm = args.fmsk
    vlsm = args.vmsk
    if flsm == "0": flsm = args.fmod ; vlsm = '_FillValue'
    #
    return args.fsat, args.vsat, args.fmod, args.vmod, flsm, vlsm, args.distgrid



if __name__ == '__main__':

    file_sat,  name_ssh_sat, file_mod, name_ssh_mod, file_lsm_mod, name_lsm_mod, l_griddist = argument_parsing()

    # Time overlap between model and satellite data ?
    (it1,it2), (Nts,Ntm) = gz.GetEpochTimeOverlap( file_sat , file_mod )
    print(' *** Time overlap between model and satellite in UNIX epoch time: it1, it2', it1,'--',it2)
    print('   => UTC: "'+gz.EpochT2Str(it1)+'" to "'+gz.EpochT2Str(it2)+'"\n')


    print('\n\n\n #####   M O D E L   2 D + T   D O M A I N   a.k.a.  S O U R C E   #####\n')

    clsm = name_lsm_mod
    if name_lsm_mod=='_FillValue': clsm = name_lsm_mod+'@'+name_ssh_mod

    ModelGrid = gz.ModGrid( file_mod, it1, it2, file_lsm_mod, clsm, distorded_grid=l_griddist )



    print('\n\n\n #####   S A T E L L I T E   1 D   T R A C K   a.k.a.  T A R G E T   #####\n')

    SatelliteTrack = gz.SatTrack( file_sat, it1, it2, Np=Nts, \
                                  domain_bounds=ModelGrid.domain_bounds, l_0_360=ModelGrid.l360 )



    # Mapping and interpolation:

    RES = gz.Model2SatTrack( ModelGrid, name_ssh_mod, SatelliteTrack, name_ssh_sat )


    # Save the result into a NetCDF file:
    import numpy as nmp
    if not IsZarr:
        from gonzag.ncio import SaveTimeSeries, Save2Dfield
    else:
        from gonzag.zarrio import SaveTimeSeries, Save2Dfield

    c1     = 'Model SSH interpolated in space (' ; c2=') and time on satellite track'
    vvar   = [ 'latitude', 'longitude', name_ssh_mod+'_np'   , name_ssh_mod+'_bl' , name_ssh_sat          , 'distance'                            ]
    vunits = [ 'deg.N'   , 'deg.E'    ,          'm'         ,     'm'            ,    'm'                ,    'km'                               ]
    vlongN = [ 'Latitude', 'Longitude', c1+'nearest-point'+c2,  c1+'bilinear'+c2  , 'Input satellite data', 'Cumulated distance from first point' ]




    iw = SaveTimeSeries( RES.time, \
                         nmp.array( [RES.lat, RES.lon, RES.ssh_mod_np, RES.ssh_mod, RES.ssh_sat, RES.distance] ), \
                         vvar, 'result.nc', \
                         time_units='seconds since 1970-01-01 00:00:00', \
                         vunits=vunits, vlnm=vlongN, missing_val=rmissval )



    if l_save_track_on_model_grid:
        xmsk = 1 - nmp.ma.getmask(RES.XNPtrack).astype(nmp.int8)
        Save2Dfield( 'xnp_msk.nc', RES.XNPtrack, xlon=ModelGrid.lon, xlat=ModelGrid.lat, name='track', mask=xmsk )
