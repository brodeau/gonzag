#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#   ///// https://github.com/brodeau/gonzag \\\\\
#
#       L. Brodeau, 2021
#
############################################################################

import time ; # to report execution speed of certain parts of the code...
#
from .config import IsZarr, ldebug, ivrb, nb_talk, l_plot_meshes, deg2km, rfactor, search_box_w_km, l_save_track_on_model_grid, l_plot_meshes, rmissval
from .utils  import *
from .bilin_mapping import BilinTrack



class Model2SatTrack:

    def __init__( self, MG, name_ssh_mod, ST, name_ssh_sat ):
        '''
        # MG: Model grid => "ModelGrid" object (util.py)
        # ST: Satellite track rid of unneeded points => "SatTrack" object (util.py)
        #
        '''
        # To be replaced with test on file extension of MG.file & ST.file:
        if not IsZarr:        
            from .ncio   import GetModel2DVar, GetSatSSH
            #Save2Dfield, SaveTimeSeries        
            
        (Nj,Ni) = MG.shape
    
        d_found_km = rfactor*MG.HResDeg*deg2km
        print(' *** "found" distance criterion when searching for nearest point on model grid is ', d_found_km, ' km\n')
    
        # Size of the search zoom box:
        np_box_radius = SearchBoxSize( MG.HResDeg*deg2km, search_box_w_km )
        print(' *** Will use zoom boxes of width of '+str(2*np_box_radius+1)+' points for 1st attempts of nearest-point location...\n')
        
        Nt = ST.size ; # number of satellit observation point to work with here...
    
        if_talk = Nt//nb_talk
        
        # The BIG GUY:
        startTime = time.time()
    
        BT = BilinTrack( ST.lat, ST.lon, MG.lat, MG.lon, src_grid_local_angle=MG.xangle, \
                         k_ew_per=MG.EWPer, rd_found_km=d_found_km, np_box_r=np_box_radius, freq_talk=if_talk )
        
        time_bl_mapping = time.time() - startTime
    
    
    
        #################################################################################################
        # Debug part to see if mapping/weights was correctly done:
        if ldebug and l_plot_meshes:
            print('\n *** ploting meshes...')
            for jt in range(Nt):
                if jt%if_talk==0:
                    print('   ==> plot for jt =', jt, ' / ', Nt)
                    [jj,ji] = BT.NP[jt,:]
                    if (jj,ji) == (-1,-1):
                        print('      ===> NO! Was not found!')
                    else:
                        PlotMesh( (ST.lon[jt],ST.lat[jt]), MG.lat, MG.lon, BT.SM[jt,:,:], BT.WB[jt,:], \
                                  fig_name='mesh_jt'+'%5.5i'%(jt)+'.png' )
        #################################################################################################
    
        # All bi-linear mapping stuff is done it's time
        # => read satellite SSH at time t_s
    
        vssh_m_np = nmp.zeros(Nt) ; vssh_m_np[:] = rmissval; # vector to store the model data interpolated in time and space (nearest-point) on the satellite track...
        vssh_m_bl = nmp.zeros(Nt) ; vssh_m_bl[:] = rmissval; # vector to store the model data interpolated in time and space (bilinear) on the satellite track...
        vdistance = nmp.zeros(Nt)
    
        # Time increment on the satellite time:
        ktm1   = 0   ; ktm2   = 0
        ktm1_o = -10 ; ktm2_o = -10

        print('\n *** Starting space-time interpolation of model data onto the '+str(Nt)+' selected track points...')
        
        startTime = time.time()
    
        for jt in range(Nt):
    
            # kt* : index for model
            # jt* : index for sat. track...
    
            itt = ST.time[jt] ; # unix time
    
            # Get surrounding records for model:
            kt = ktm1
            while not (MG.time[kt]<=itt and MG.time[kt+1]>itt): kt=kt+1
            ktm1 = kt ; ktm2 = kt+1
    
            if jt%if_talk==0:
                print('      jt = '+'%5.5i'%(jt)+' => satelite time = '+EpochT2Str(itt))
                if ivrb>0: print('   => surounding kt for model: ', ktm1, ktm2, '(',EpochT2Str(MG.time[ktm1]),EpochT2Str(MG.time[ktm2]),') / ', \
                      MG.time[ktm1],MG.time[ktm2] )
    
            # If first time we have these ktm1 & ktm2, getting the two surrounding fields:
            if (ktm1>ktm1_o) and (ktm2>ktm2_o):
                if (ktm1_o == -10) or (ktm1 > ktm2_o):
                    if ivrb>0: print(' *** Reading '+name_ssh_mod+' in '+MG.file+'\n    => at ktm1=', ktm1)
                    Xm1 = GetModel2DVar( MG.file, name_ssh_mod, kt=ktm1 )
                else:
                    Xm1[:,:] = Xm2[:,:]
                #
                if ivrb>0: print(' *** Reading '+name_ssh_mod+' in '+MG.file+'\n    => at ktm2=', ktm2)
                Xm2 = GetModel2DVar( MG.file, name_ssh_mod, kt=ktm2 )
    
                # slope only needs to be calculated when Xm2 and Xm1 have been updated:
                Xa = (Xm2 - Xm1) / float(MG.time[ktm2] - MG.time[ktm1])
    
            # Linear interpolation of field at time itt:
            if ivrb>0 and jt%if_talk==0: print('   => Model data is interpolated at current time out of model records '+str(ktm1)+' & '+str(ktm2))
            Xm = Xm1[:,:] + Xa[:,:]*float(itt - MG.time[ktm1])
    
            [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = BT.SM[jt,:,:]
            [w1, w2, w3, w4]                    = BT.WB[jt,:]
    
            Sm = MG.mask[j1,i1] + MG.mask[j2,i2] + MG.mask[j3,i3] + MG.mask[j4,i4]
            if Sm == 4:
                # All 4 model points are ocean point !
    
                vssh_m_np[jt] = Xm[j1,i1] ; # Nearest-point "interpolation"
    
                # Bilinear interpolation:
                Sw = nmp.sum([w1, w2, w3, w4])
                if abs(Sw-1.)> 0.001:
                    if ivrb>0: print('    FLAGGING MISSING VALUE at jt = '+str(jt)+' !!!')
                else:
                    vssh_m_bl[jt] = w1*Xm[j1,i1] + w2*Xm[j2,i2] + w3*Xm[j3,i3] + w4*Xm[j4,i4]
    
            ktm1_o = ktm1 ; ktm2_o = ktm2
    
        # end of loop on jt
        
        time_bl_interp = time.time() - startTime

        print(' *** Space-time interpolation done!\n')

        
        # Distance parcourue since first point:
        for jt in range(1,Nt):
            vdistance[jt] = vdistance[jt-1] + haversine_sclr( ST.lat[jt], ST.lon[jt], ST.lat[jt-1], ST.lon[jt-1] )
    
        # Satellite SSH:
        vssh_s = GetSatSSH( ST.file, name_ssh_sat,  kt1=ST.jt1, kt2=ST.jt2, ikeep=ST.keepit )
                
        self.time       = ST.time
        self.lat        = ST.lat        
        self.lon        = ST.lon
        self.ssh_mod_np = vssh_m_np
        self.ssh_mod_bl = vssh_m_bl
        self.ssh_sat    = vssh_s
        self.distance   = vdistance
        
        print('\n *** Time report:')
        print('     - Construction of the source-target bilinear mapping took: '+str(round(time_bl_mapping,0))+' s')
        print('     - Interpolation of model data on the '+str(Nt)+' satellite points took: '+str(round(time_bl_interp,0))+' s \n')


        if l_save_track_on_model_grid:
            # Save the satellite nearest-point track on the model grid:
            xnp = nmp.zeros((Nj,Ni)) ; xnp[:,:] = rmissval
            for jt in range(Nt):
                [jj,ji] = BT.NP[jt,:]
                xnp[jj,ji] = float(jt)
            xnp[nmp.where(MG.mask==0)] = -100.
            xmsk = nmp.zeros((Nj,Ni))
            xmsk[nmp.where(xnp>-110.)] = 1
            #xnp = nmp.ma.masked_where(xnp>-110.)
            #
            self.XNPtrack = xnp
            self.XNPmask  = xmsk
