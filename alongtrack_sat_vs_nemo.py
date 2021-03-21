#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

from sys import exit
import numpy as nmp
from netCDF4 import Dataset
import time ; # to check the execution speed...
#
import gonzag as gz
from gonzag.config import ldebug

# Should be command-line arguments:
np_box_radius = 4 ; # in number of points... (should give km and get this one based on ocean model res...)
dist_found = 100. ; # threshold distance for found [km] ! ""    ""
l_dump_np_track_on_model_grid = True
ew_per=2 ; # east-west periodicity of source/model grid...
if_talk = 500 ; # verbose frequency: we chat every if_talk time steps !!!



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
    parser.add_argument('-l', '--fmsk' , default="0",              help='specify model grid land-sea mask (<ncfile> or "0" for _FillValue)')
    parser.add_argument('-k', '--vmsk' , default="tmask",          help='variable name of model grid land-sea mask in')
    #
    args = parser.parse_args()
    print('')
    print(' *** Satellinte file and variable:\n', '  => ', args.fsat, '"'+args.vsat+'"')
    print(' ***    Model   file and variable:\n', '  => ', args.fmod, '"'+args.vmod,'"\n')
    #
    flsm = args.fmsk
    vlsm = args.vmsk
    if flsm == "0": flsm = args.fmod ; vlsm = '_FillValue'
    #
    return args.fsat, args.vsat, args.fmod, args.vmod, flsm, vlsm















if __name__ == '__main__':

    file_sat,  name_ssh_sat, file_mod, name_ssh_mod, file_lsm, name_lsm = argument_parsing()

    # Checking existence of input files:
    gz.chck4f(file_sat)
    gz.chck4f(file_mod)
    if name_lsm=='_FillValue':
        clsm = name_lsm+'@'+name_ssh_mod
    else:
        clsm = name_lsm
        gz.chck4f(file_lsm)
    print('\n Satellite: '+file_sat+'\n Model: '+file_mod+'\n')

    # Calendar/time vectors:
    vdate_s, itime_s = gz.GetTimeVector( file_sat )
    vdate_m, itime_m = gz.GetTimeVector( file_mod )
    
    # Work period, i.e. the time overlap for the two time-series:
    it_min, it_max = gz.GetTimeOverlapBounds( itime_s, itime_m )
    print(' *** Time overlap for Model and Track data:', it_min, it_max,'\n')


    # Get model coordinates and land-sea mask:
    xlat_m = gz.GetModelCoor( file_mod,  'latitude' )
    xlon_m = gz.GetModelCoor( file_mod, 'longitude' )
    mask_m = gz.GetModelLSM( file_lsm, clsm ) ; # land-sea mask...
    (Nj,Ni) = xlat_m.shape
    if xlon_m.shape != (Nj,Ni) or mask_m.shape != (Nj,Ni): MsgExit('shape disagreement for model arrays')
    xlon_m = nmp.mod(xlon_m, 360.) ; # forces longitude to be in the [0,360] range...
    #if ldebug: Save2Dfield( 'mask_model.nc', mask_m, name='mask' ) #lolodbg

    # Get distortion angle of model grid:
    cf_angle = 'grid_angle_model.nc'
    xangle_m = gz.GridAngle( xlat_m, xlon_m )
    #xangle_m = nmp.zeros(xlat_m.shape) ; print(' NOOOOOOT    => done!\n'); #lolodbg
    if ldebug: Save2Dfield( cf_angle, xangle_m, name='angle', mask=mask_m )


    # Relevant satellite time slice:
    jts_1, jts_2 = gz.scan_idx_sat( itime_s, it_min, it_max )
    Nt_s = jts_2 - jts_1 + 1
    if ldebug:
        print(' jts_1 =', jts_1, '  ==> itime_s[jts_1] =', itime_s[jts_1] )
        print(' jts_2 =', jts_2, '  ==> itime_s[jts_2] =', itime_s[jts_2],'\n' )

    # Get satellite track lat,lon for the relevant time slice:
    vlat_s = gz.GetSatCoord( file_sat, jts_1,jts_2, 'latitude'  )
    vlon_s = gz.GetSatCoord( file_sat, jts_1,jts_2, 'longitude' )

    #print('BEFORE: itime_s[0]=',itime_s[0])
    itime_s = itime_s[jts_1:jts_2+1] ; # trimming satellite time vector to usefull period
    #print('AFTER: itime_s[0]=',itime_s[0])    
    if len(itime_s)!=Nt_s or len(vlat_s)!=Nt_s or len(vlon_s)!=Nt_s: MsgExit('problem with satellite record length')
    
    
    ##LOLOdbg: test SaveTimeSeriesNC :
    #vd1 = nmp.zeros(Nt_s) ; vd1[:] = nmp.sin(itime_s.astype(nmp.float32))
    #vd2 = nmp.zeros(Nt_s) ; vd2[:] = nmp.cos(itime_s.astype(nmp.float32))
    #vvar   = ['sin','cos']
    #vunits = ['bullshit', 'bullshit']
    #vlongN = ['Sine of unix time','Cosine of unix time']
    #iw = SaveTimeSeriesNC( itime_s, nmp.array([vd1,vd2]), vvar, 'test.nc', time_units='seconds since 1970-01-01 00:00:00', \
    #                       vunits=vunits, vlnm=vlongN, missing_val=-9999. )
    #
    # Getting satellite data to include in file to write:
    #id_sat   = Dataset(file_sat)
    #vssh_s = id_sat.variables[name_ssh_sat][jts_1:jts_2+1]
    #id_sat.close()
    #
    #vssh_m_bl = nmp.zeros(Nt_s) ; # vector to store the model data interpolated in time and space (bilinear) on the satellite track...
    #vssh_m_np = nmp.zeros(Nt_s) ; # vector to store the model data interpolated in time and space (nearest-point) on the satellite track...
    #
    #c1 = 'Model SSH interpolated in space (' ; c2=') and time on satellite track'
    #vvar   = [ name_ssh_mod+'_bl', name_ssh_mod+'_np'    , name_ssh_sat           ]
    #vunits = [     'm'           ,          'm'          ,    'm'                 ]
    #vlongN = [ c1+'bilinear'+c2  , c1+'nearest-point'+c2 , 'Input satellite data' ] 
    #
    #iw = SaveTimeSeriesNC( itime_s, nmp.array( [vssh_m_bl, vssh_m_np, vssh_s] ), vvar, 'test.nc', \
    #                       time_units='seconds since 1970-01-01 00:00:00', \
    #                       vunits=vunits, vlnm=vlongN, missing_val=-9999. )
    #exit(0)

    
    jt_stop = 10000 ; Nt_s = len(vlat_s[:jt_stop]) ; # lolodbg
    #jt_stop = len(vlat_s)

    
    startTime = time.time()

    BT = gz.BilinTrack( vlat_s[:jt_stop], vlon_s[:jt_stop], xlat_m, xlon_m, src_grid_local_angle=xangle_m, \
                        k_ew_per=2, rd_found_km=dist_found, np_box_r=np_box_radius, freq_talk=if_talk )

    print('\n *** Execution time to build the mapping: ', time.time() - startTime, '\n')

    #################################################################################################
    # Debug part to see if mapping/weights was correctly done:
    if ldebug:
        print('\n *** ploting meshes...')
        for jt in range(len(vlat_s[:jt_stop])):
            if jt%if_talk==0:
                print('   ==> plot for jt =', jt, ' / ', len(vlat_s[:jt_stop]))
                gz.PlotMesh( (vlon_s[jt],vlat_s[jt]), xlat_m, xlon_m, BT.SM[jt,:,:], BT.WB[jt,:], \
                             fig_name='mesh_jt'+'%5.5i'%(jt)+'.png' )
    #################################################################################################

    if l_dump_np_track_on_model_grid:
        # Show the satellite track on the model grid:
        xnp_msk = nmp.zeros((Nj,Ni)) ; xnp_msk[:,:] = -9999.
        for jt in range(len(vlat_s[:jt_stop])):
            [jj,ji] = BT.NP[jt,:]
            xnp_msk[jj,ji] = float(jt)
        xnp_msk[nmp.where(mask_m==0)] = -100.
        xmsk_tmp = nmp.zeros((Nj,Ni))
        xmsk_tmp[nmp.where(xnp_msk>-110.)] = 1
        Save2Dfield( 'xnp_msk.nc', xnp_msk, xlon=xlon_m, xlat=xlat_m, name='nb', mask=xmsk_tmp )
        del xmsk_tmp



    # All bi-lin mapping stuff is done it's time
    # => read satellite SSH at time t_s

    vssh_m_np = nmp.zeros(Nt_s) ; # vector to store the model data interpolated in time and space (nearest-point) on the satellite track...
    vssh_m_bl = nmp.zeros(Nt_s) ; # vector to store the model data interpolated in time and space (bilinear) on the satellite track...


    
    # Time increment on the satellite time:
    id_mod = Dataset(file_mod)    
    jtm1   = 0   ; jtm2   = 0
    jtm1_o = -10 ; jtm2_o = -10
    
    for jts in range(Nt_s):

        rtd = vdate_s[jts] ; # formatted
        itt = itime_s[jts] ; # unix time

        # Get surrounding records for model:
        jt = jtm1
        while not (itime_m[jt]<=itt and itime_m[jt+1]>itt): jt=jt+1
        jtm1 = jt ; jtm2 = jt+1
        
        print('\n *** Satelite time at jt = '+'%5.5i'%(jts)+' ==> ',rtd, itt)
        print('   => surounding jt for model: ', jtm1, jtm2, '(',vdate_m[jtm1],vdate_m[jtm2],') / ', \
              itime_m[jtm1],itime_m[jtm2] )

        # If first time we have these jtm1 & jtm2, getting the two surrounding fields:
        if (jtm1>jtm1_o) and (jtm2>jtm2_o):
            if (jtm1_o == -10) or (jtm1 > jtm2_o):
                print(' *** Reading '+name_ssh_mod+' in '+file_mod)
                print('    => at jtm1=', jtm1)
                #lilo
                xvar1 = id_mod.variables[name_ssh_mod][jtm1,:,:]
                #if l_use_anomaly: xvar1 = xvar1 - xmean
            else:
                xvar1[:,:] = xvar2[:,:]
            #
            print(' *** Reading '+name_ssh_mod+' in '+file_mod)
            print('    => at jtm2=', jtm2)
            xvar2 = id_mod.variables[name_ssh_mod][jtm2,:,:]
            #if l_use_anomaly ) xvar2 = xvar2 - xmean

            # slope only needs to be calculated when xvar2 and xvar1 have been updated:
            xslope = (xvar2 - xvar1) / float(itime_m[jtm2] - itime_m[jtm1])

        # Linear interpolation of field at time itt:
        print('   => Model data is interpolated at current time out of model records '+str(jtm1+1)+' & '+str(jtm2+1))
        xvar = xvar1[:,:] + xslope[:,:]*float(itt - itime_m[jtm1])


        [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = BT.SM[jts,:,:]
        [w1, w2, w3, w4]                    = BT.WB[jts,:]

        # Nearest-point "interpolation":
        vssh_m_np[jts] = xvar[j1,i1]
        
        # Bilinear interpolation:
        Sw = nmp.sum([w1, w2, w3, w4])
        if abs(Sw-1.)> 0.001:
            print('    FLAGGING MISSING VALUE at jts = '+str(jts)+' !!!')
            vssh_m_bl[jts] = -9999.
            #MsgExit('Sum of weights not equal to 1.: '+str(Sw))
        else:
            vssh_m_bl[jts] = w1*xvar[j1,i1] + w2*xvar[j2,i2] + w3*xvar[j3,i3] + w4*xvar[j4,i4]
        
        jtm1_o = jtm1 ; jtm2_o = jtm2

    # end of loop on jts
    id_mod.close()




    # Getting satellite data for reference purposes:
    id_sat = Dataset(file_sat)
    vssh_s = id_sat.variables[name_ssh_sat][jts_1:jts_1+Nt_s]
    id_sat.close()

    #print('nmp.shape(vssh_s) =', nmp.shape(vssh_s))    
    #exit(0)





    print('LOLO: nmp.shape(vssh_m_bl) = ', nmp.shape(vssh_m_bl))
    print('LOLO: nmp.shape(vssh_m_np) = ', nmp.shape(vssh_m_np))
    print('LOLO: nmp.shape(vssh_s) = ', nmp.shape(vssh_s))

        
    c1 = 'Model SSH interpolated in space (' ; c2=') and time on satellite track'
    vvar   = [ name_ssh_mod+'_bl', name_ssh_mod+'_np'    , name_ssh_sat           ]
    vunits = [     'm'           ,          'm'          ,    'm'                 ]
    vlongN = [ c1+'bilinear'+c2  , c1+'nearest-point'+c2 , 'Input satellite data' ] 


    print('LOLO: nmp.shape(nmp.array( [vssh_m_bl, vssh_m_np, vssh_s] )) = ', nmp.shape(nmp.array( [vssh_m_bl, vssh_m_np, vssh_s] )) )
    
    iw = SaveTimeSeriesNC( itime_s[:Nt_s], nmp.array( [vssh_m_bl, vssh_m_np, vssh_s] ), vvar, 'test.nc', \
                           time_units='seconds since 1970-01-01 00:00:00', \
                           vunits=vunits, vlnm=vlongN, missing_val=-9999. )



    exit(0)
