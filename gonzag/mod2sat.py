#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

#from sys import exit
#import numpy as nmp
#from netCDF4 import Dataset
import time ; # to check the execution speed...
#
#import gonzag as gz
from .config import deg2km, ldebug, rmissval
from .utils  import *
from .ncio   import *
from .bilin_mapping import BilinTrack

# Should be command-line arguments:
np_box_radius = 4 ; # in number of points... (should give km and get this one based on ocean model res...)
l_dump_np_track_on_model_grid = True
if_talk = 500 ; # verbose frequency: we chat every if_talk time steps !!!



def Model2SatTrack( file_sat,  name_ssh_sat, file_mod, name_ssh_mod, file_lsm_mod, name_lsm_mod, \
                    ew_prd_mod=-1, file_out='mod2sat.nc' ):

    # Checking existence of input files:
    chck4f(file_sat)
    chck4f(file_mod)
    if name_lsm_mod=='_FillValue':
        clsm = name_lsm_mod+'@'+name_ssh_mod
    else:
        clsm = name_lsm_mod
        chck4f(file_lsm_mod)
    print('\n Satellite: '+file_sat+'\n Model: '+file_mod+'\n')


    # LIMITATION #1: TIME => overlap time-periode beteen model and satellite data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    vdate_s, itime_s = GetTimeVector( file_sat )
    vdate_m, itime_m = GetTimeVector( file_mod )
    
    # Work period, i.e. the time overlap for the two time-series, if any...
    it_min, it_max = GetTimeOverlapBounds( itime_s, itime_m )
    print(' *** Time overlap for Model and Track data:', it_min, it_max,'\n')

    # Relevant period in terms of record index for satellite data:
    jts_1, jts_2 = scan_idx_sat( itime_s, it_min, it_max )
    Nt_s = jts_2 - jts_1 + 1
    if ldebug:
        print(' jts_1 =', jts_1, '  ==> itime_s[jts_1] =', itime_s[jts_1] )
        print(' jts_2 =', jts_2, '  ==> itime_s[jts_2] =', itime_s[jts_2],'\n' )



    # LIMITATION #2: SPACE => if model is a regional box then a lot can be removed from satellite data
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Get model coordinates and land-sea mask:
    xlat_m = GetModelCoor( file_mod,  'latitude' )
    xlon_m = GetModelCoor( file_mod, 'longitude' )
    mask_m = GetModelLSM( file_lsm_mod, clsm ) ; # land-sea mask...
    (Nj,Ni) = xlat_m.shape
    if xlon_m.shape != (Nj,Ni) or mask_m.shape != (Nj,Ni): MsgExit('shape disagreement for model arrays')
    xlon_m = nmp.mod(xlon_m, 360.) ; # forces longitude to be in the [0,360] range...
    #if ldebug: Save2Dfield( 'mask_model.nc', mask_m, name='mask' ) #lolodbg

    # Rough estimate of the resolution of the model grid:
    res_model_dg = GetModelResolution( xlon_m )
    if res_model_dg>5. or res_model_dg<0.001: MsgExit('Model resolution found is surprising, prefer to stop => check "GetModelResolution()" in utils.py')
    d_found_km = 0.65*res_model_dg*deg2km
    print('   ==> "found" distance criterion when searching for nearest point is ', d_found_km, ' km\n')

    
    # Global or regional config ?
    l_glob_lon_wize, l360, lon_min, lon_max = IsGlobalLongitudeWise( xlon_m, resd=res_model_dg )
    lat_min = nmp.amin(xlat_m) ; lat_max = nmp.amax(xlat_m)
    #
    cw = 'regional'
    if l_glob_lon_wize: cw = 'global'
    print(' *** Seems like the model domain is '+cw+' (in terms of longitude span)...')
    print('     => lon_min, lon_max = ', lon_min, lon_max, '\n')

    if ew_prd_mod>=0 and not l_glob_lon_wize:
        print('\n  WARNING: forcing East-West periodicity to NONE (ew_prd_mod=-1) because regional domain!\n')
    
    # Get satellite data during the relevant time slice:
    vdate_s = vdate_s[jts_1:jts_2+1] ; # trimming satellite time vector to usefull period
    itime_s = itime_s[jts_1:jts_2+1] ; # trimming satellite time vector to usefull period
    #
    vlat_s  = GetSatCoord( file_sat, jts_1,jts_2, 'latitude'  )
    vlon_s  = GetSatCoord( file_sat, jts_1,jts_2, 'longitude' )
    if l360:
        vlon_s = nmp.mod(vlon_s, 360.)
    else:
        vlon_s = degE_to_degWE( vlon_s )
    #
    id_sat = Dataset(file_sat)
    vssh_s = id_sat.variables[name_ssh_sat][jts_1:jts_1+Nt_s]
    id_sat.close()
    #
    if itime_s.shape!=(Nt_s,) or vlat_s.shape!=(Nt_s,) or vlon_s.shape!=(Nt_s,) or vssh_s.shape!=(Nt_s,):
        MsgExt('satellite arrays do not agree in shape after time slicing')
               
    print(' *** Track size before removing points outside of model domain: '+str(len(vlon_s)))
    keep_lat = nmp.where( (vlat_s[:]>=lat_min) & (vlat_s[:]<=lat_max) & (vlon_s[:]>=lon_min) & (vlon_s[:]<=lon_max) )    
    itime_s = itime_s[keep_lat]
    vlat_s  = vlat_s[keep_lat]
    vlon_s  = vlon_s[keep_lat]
    vssh_s  = vssh_s[keep_lat]
    (Nt_s,) = vssh_s.shape

    print('   => Track size AFTER removing points outside of model domain: '+str(len(vlon_s))+'\n')

    
    
    # Get distortion angle of model grid:
    cf_angle = 'grid_angle_model.nc'
    xangle_m = GridAngle( xlat_m, xlon_m )    
    if ldebug: Save2Dfield( cf_angle, xangle_m, name='angle', mask=mask_m )
    #xangle_m = nmp.zeros(xlat_m.shape) ; print(' NOOOOOOT    => done!\n'); #lolodbg

    
    #jt_stop = 10000 ; Nt_s = len(vlat_s[:jt_stop]) ; # lolodbg
    jt_stop = Nt_s
    
    startTime = time.time()

    BT = BilinTrack( vlat_s[:jt_stop], vlon_s[:jt_stop], xlat_m, xlon_m, src_grid_local_angle=xangle_m, \
                        k_ew_per=ew_prd_mod, rd_found_km=d_found_km, np_box_r=np_box_radius, freq_talk=if_talk )

    print('\n *** Execution time to build the mapping: ', time.time() - startTime, '\n')

    #################################################################################################
    # Debug part to see if mapping/weights was correctly done:
    if ldebug:
        print('\n *** ploting meshes...')
        for jt in range(len(vlat_s[:jt_stop])):
            if jt%if_talk==0:
                print('   ==> plot for jt =', jt, ' / ', len(vlat_s[:jt_stop]))
                PlotMesh( (vlon_s[jt],vlat_s[jt]), xlat_m, xlon_m, BT.SM[jt,:,:], BT.WB[jt,:], \
                             fig_name='mesh_jt'+'%5.5i'%(jt)+'.png' )
    #################################################################################################
    
    if l_dump_np_track_on_model_grid:
        # Show the satellite track on the model grid:
        xnp_msk = nmp.zeros((Nj,Ni)) ; xnp_msk[:,:] = rmissval
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

    vssh_m_np[:] = rmissval
    vssh_m_bl[:] = rmissval

    
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

        Sm = mask_m[j1,i1] + mask_m[j2,i2] + mask_m[j3,i3] + mask_m[j4,i4]

        if Sm == 4:
            # All 4 model points are ocean point !
            
            # Nearest-point "interpolation":
            vssh_m_np[jts] = xvar[j1,i1]
            
            # Bilinear interpolation:
            Sw = nmp.sum([w1, w2, w3, w4])
            if abs(Sw-1.)> 0.001:
                print('    FLAGGING MISSING VALUE at jts = '+str(jts)+' !!!')
            else:
                vssh_m_bl[jts] = w1*xvar[j1,i1] + w2*xvar[j2,i2] + w3*xvar[j3,i3] + w4*xvar[j4,i4]
        
        jtm1_o = jtm1 ; jtm2_o = jtm2

    # end of loop on jts
    id_mod.close()




    # Getting satellite data for reference purposes:

    #print('nmp.shape(vssh_s) =', nmp.shape(vssh_s))    
    #exit(0)


    c1 = 'Model SSH interpolated in space (' ; c2=') and time on satellite track'
    vvar   = [ name_ssh_mod+'_bl', name_ssh_mod+'_np'    , name_ssh_sat           ]
    vunits = [     'm'           ,          'm'          ,    'm'                 ]
    vlongN = [ c1+'bilinear'+c2  , c1+'nearest-point'+c2 , 'Input satellite data' ] 


    if nmp.ma.is_masked(vssh_s):
        idxma = nmp.where( nmp.ma.getmask(vssh_s) )
        vssh_s[idxma] = rmissval

    iw = SaveTimeSeriesNC( itime_s[:Nt_s], nmp.array( [vssh_m_bl, vssh_m_np, vssh_s] ), vvar, file_out, \
                              time_units='seconds since 1970-01-01 00:00:00', \
                              vunits=vunits, vlnm=vlongN, missing_val=rmissval )



    return 1
