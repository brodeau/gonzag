#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

import time ; # to report execution speed of certain parts of the code...
#
from .config import ldebug, if_talk, l_plot_meshes, deg2km, rfactor, search_box_w_km, l_dump_np_track_on_model_grid, l_plot_meshes, rmissval
from .utils  import *
from .ncio   import *
from .bilin_mapping import BilinTrack

# Should be command-line arguments:
#np_box_radius = 4 ; # in number of points... (should give km and get this one based on ocean model res...)
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
    jt_1, jt_2 = scan_idx_sat( itime_s, it_min, it_max )
    Nt_s = jt_2 - jt_1 + 1
    if ldebug:
        print(' jt_1 =', jt_1, '  ==> itime_s[jt_1] =', itime_s[jt_1] )
        print(' jt_2 =', jt_2, '  ==> itime_s[jt_2] =', itime_s[jt_2],'\n' )



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
    d_found_km = rfactor*res_model_dg*deg2km
    print('   ==> "found" distance criterion when searching for nearest point is ', d_found_km, ' km\n')


    # Global or regional config ?
    l_glob_lon_wize, l360, lon_min, lon_max = IsGlobalLongitudeWise( xlon_m, resd=res_model_dg )
    lat_min = nmp.amin(xlat_m) ; lat_max = nmp.amax(xlat_m)
    cw = 'regional'
    if l_glob_lon_wize: cw = 'global'
    print('     => lat_min, lat_max = ', lat_min, lat_max)
    print('     => lon_min, lon_max = ', lon_min, lon_max, '\n')
    print(' *** Seems like the model domain is '+cw+' (in terms of longitude span)...')
    if ew_prd_mod>=0 and not l_glob_lon_wize:
        print('\n  WARNING: forcing East-West periodicity to NONE (ew_prd_mod=-1) because regional domain!\n')


    # Size of the search zoom box:
    np_box_radius = SearchBoxSize( res_model_dg*deg2km, search_box_w_km )
    print(' *** Will use zoom boxes of width of '+str(2*np_box_radius+1)+' points for 1st attempts of nearest-point location...\n')


    # Get satellite data during the relevant time slice:
    vdate_s = vdate_s[jt_1:jt_2+1] ; # trimming satellite time vector to usefull period
    itime_s = itime_s[jt_1:jt_2+1] ; # trimming satellite time vector to usefull period
    #
    vlat_s  = GetSatCoord( file_sat, jt_1,jt_2, 'latitude'  )
    vlon_s  = GetSatCoord( file_sat, jt_1,jt_2, 'longitude' )
    if l360:
        vlon_s = nmp.mod(vlon_s, 360.)
    else:
        vlon_s = degE_to_degWE( vlon_s )
    #
    id_sat = Dataset(file_sat)
    vssh_s = id_sat.variables[name_ssh_sat][jt_1:jt_1+Nt_s]
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



    # The BIG GUY:
    startTime = time.time()

    BT = BilinTrack( vlat_s[:jt_stop], vlon_s[:jt_stop], xlat_m, xlon_m, src_grid_local_angle=xangle_m, \
                     k_ew_per=ew_prd_mod, rd_found_km=d_found_km, np_box_r=np_box_radius, freq_talk=if_talk )

    time_bl_mapping = time.time() - startTime



    #################################################################################################
    # Debug part to see if mapping/weights was correctly done:
    if ldebug and l_plot_meshes:
        print('\n *** ploting meshes...')
        for jt in range(len(vlat_s[:jt_stop])):
            if jt%if_talk==0:
                print('   ==> plot for jt =', jt, ' / ', len(vlat_s[:jt_stop]))
                [jj,ji] = BT.NP[jt,:]
                if (jj,ji) == (-1,-1):
                    print('      ===> NO! Was not found!')
                else:
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
        Save2Dfield( 'xnp_msk.nc', xnp_msk, xlon=xlon_m, xlat=xlat_m, name='track', mask=xmsk_tmp )
        del xmsk_tmp



    # All bi-linear mapping stuff is done it's time
    # => read satellite SSH at time t_s

    vssh_m_np = nmp.zeros(Nt_s) ; vssh_m_np[:] = rmissval; # vector to store the model data interpolated in time and space (nearest-point) on the satellite track...
    vssh_m_bl = nmp.zeros(Nt_s) ; vssh_m_bl[:] = rmissval; # vector to store the model data interpolated in time and space (bilinear) on the satellite track...
    vdistance = nmp.zeros(Nt_s)

    # Time increment on the satellite time:
    id_mod = Dataset(file_mod)
    ktm1   = 0   ; ktm2   = 0
    ktm1_o = -10 ; ktm2_o = -10

    startTime = time.time()

    for jt in range(Nt_s):

        # kt* : index for model
        # jt* : index for sat. track...

        rtd = vdate_s[jt] ; # formatted
        itt = itime_s[jt] ; # unix time

        # Get surrounding records for model:
        kt = ktm1
        while not (itime_m[kt]<=itt and itime_m[kt+1]>itt): kt=kt+1
        ktm1 = kt ; ktm2 = kt+1

        if jt%if_talk==0:
            print('\n *** Satelite time at jt = '+'%5.5i'%(jt)+' ==> ',rtd, itt)
            print('   => surounding kt for model: ', ktm1, ktm2, '(',vdate_m[ktm1],vdate_m[ktm2],') / ', \
                  itime_m[ktm1],itime_m[ktm2] )

        # If first time we have these ktm1 & ktm2, getting the two surrounding fields:
        if (ktm1>ktm1_o) and (ktm2>ktm2_o):
            if (ktm1_o == -10) or (ktm1 > ktm2_o):
                print(' *** Reading '+name_ssh_mod+' in '+file_mod+'\n    => at ktm1=', ktm1)
                Xm1 = id_mod.variables[name_ssh_mod][ktm1,:,:]
            else:
                Xm1[:,:] = Xm2[:,:]
            #
            print(' *** Reading '+name_ssh_mod+' in '+file_mod+'\n    => at ktm2=', ktm2)
            Xm2 = id_mod.variables[name_ssh_mod][ktm2,:,:]

            # slope only needs to be calculated when Xm2 and Xm1 have been updated:
            Xa = (Xm2 - Xm1) / float(itime_m[ktm2] - itime_m[ktm1])

        # Linear interpolation of field at time itt:
        if jt%if_talk==0: print('   => Model data is interpolated at current time out of model records '+str(ktm1+1)+' & '+str(ktm2+1))
        Xm = Xm1[:,:] + Xa[:,:]*float(itt - itime_m[ktm1])

        [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = BT.SM[jt,:,:]
        [w1, w2, w3, w4]                    = BT.WB[jt,:]

        Sm = mask_m[j1,i1] + mask_m[j2,i2] + mask_m[j3,i3] + mask_m[j4,i4]
        if Sm == 4:
            # All 4 model points are ocean point !

            vssh_m_np[jt] = Xm[j1,i1] ; # Nearest-point "interpolation"

            # Bilinear interpolation:
            Sw = nmp.sum([w1, w2, w3, w4])
            if abs(Sw-1.)> 0.001:
                print('    FLAGGING MISSING VALUE at jt = '+str(jt)+' !!!')
            else:
                vssh_m_bl[jt] = w1*Xm[j1,i1] + w2*Xm[j2,i2] + w3*Xm[j3,i3] + w4*Xm[j4,i4]

        ktm1_o = ktm1 ; ktm2_o = ktm2

    # end of loop on jt
    id_mod.close()

    time_bl_interp = time.time() - startTime

    # Distance parcourue since first point:
    for jt in range(1,Nt_s):
        vdistance[jt] = vdistance[jt-1] + haversine_sclr( vlat_s[jt], vlon_s[jt], vlat_s[jt-1], vlon_s[jt-1] )


    # Output file

    c1 = 'Model SSH interpolated in space (' ; c2=') and time on satellite track'
    vvar   = [ 'latitude', 'longitude', name_ssh_mod+'_bl', name_ssh_mod+'_np'    , name_ssh_sat           , 'distance' ]
    vunits = [ 'deg.N'   , 'deg.E'    ,   'm'           ,          'm'          ,    'm'                   ,    'km'    ]
    vlongN = [ 'Latitude', 'Longitude',  c1+'bilinear'+c2  , c1+'nearest-point'+c2 , 'Input satellite data', 'Cumulated distance since first point' ]

    if nmp.ma.is_masked(vssh_s):
        idxma = nmp.where( nmp.ma.getmask(vssh_s) )
        vssh_s[idxma] = rmissval

    iw = SaveTimeSeries( itime_s[:Nt_s], nmp.array( [vlat_s, vlon_s, vssh_m_bl, vssh_m_np, vssh_s, vdistance] ), vvar, file_out, \
                         time_units='seconds since 1970-01-01 00:00:00', \
                         vunits=vunits, vlnm=vlongN, missing_val=rmissval )

    print('\n *** Time report:')
    print('     - Construction of the source-target bilinear mapping took: '+str(round(time_bl_mapping,0))+' s')
    print('     - Interpolation of model data on the '+str(Nt_s)+' satellite points took: '+str(round(time_bl_interp,0))+' s \n')

    return 1


