#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#       L. Brodeau, 2021
############################################################################

from sys import exit
from math import radians, cos, sin, asin, sqrt, pi, tan, log, atan2, copysign
import numpy as nmp


def MsgExit( cmsg ):
    print('\n ERROR: '+cmsg+' !\n')
    exit(0)

def chck4f( ncfile ):
    from os.path import exists
    if not exists(ncfile):
        MsgExit('File '+ncfile+' does not exist')

def degE_to_degWE( x ):
    '''
    # From longitude in 0 -- 360 frame to -180 -- +180 frame...
    '''
    return copysign(1., 180.-x)*min(x, abs(x-360.))

def degE_to_degWE_1d( vx ):
    '''
    # From longitude in 0 -- 360 frame to -180 -- +180 frame...
    '''
    return nmp.copysign(1., 180.-vx[:])*nmp.minimum(vx[:], nmp.abs(vx[:]-360.))

def find_j_i_min(x):
    '''
    # Yes, reinventing the wheel here, but it turns out
    # it is faster this way!
    '''
    k = x.argmin()
    ny = x.shape[1]
    return k//ny, k%ny




def GetTimeOverlapBounds( itsat, itmod ):
    kt1_s = itsat[0]
    kt2_s = itsat[-1]
    kt1_m = itmod[0]
    kt2_m = itmod[-1]
    #print('\n *** Earliest/latest dates:\n => for satellite data:',kt1_s,kt2_s,'\n => for model     data:',kt1_m,kt2_m)
    if (kt1_m >= kt2_s) or (kt1_s >= kt2_m) or (kt2_m <= kt1_s) or (kt2_s <= kt1_m):
        MsgExit('No time overlap for Model and Track file')
    return max(kt1_s, kt1_m), min(kt2_s, kt2_m)


def scan_idx_sat( itsat, it1, it2 ):
    '''
    # Finding indices when we can start and stop when scanning the track file:
    '''
    Nt0 = len(itsat)
    for kt1 in range(  0, Nt0-1):
        if  (itsat[kt1] <= it1) and (itsat[kt1+1] > it1): break
    for kt2 in range(kt1, Nt0-1):
         if (itsat[kt2] <= it2) and (itsat[kt2+1] > it2): break
    kt2 = kt2 + 1
    #
    return kt1, kt2






def GridAngle( xlat, xlon ):
    ''' To be used with a NEMO ORCA-type of grid, to get an idea of the local distortion (rotation)
    #   of the grid
    # Returns local distortion of the grid in degrees [-180,180]
    '''
    to_rad = pi/180.
    pio4   = pi/4.
    (Ny,Nx) = nmp.shape(xlat)
    #
    print(' *** Computing angle distortion of the model grid...')
    angle = nmp.zeros((Ny,Nx))
    for ji in range(Nx):
        for jj in range(1,Ny-1):

            if abs( xlon[jj+1,ji]%360. - xlon[jj-1,ji]% 360. ) < 1.e-8:
                sint = 0.
                cost = 1.
                #
            else:
                zt0 = tan( pio4 - to_rad*xlat[jj,ji]/2. )
                # North pole direction & modulous (at t-point)
                zxnpt = 0. - 2. * cos( to_rad*xlon[jj,ji] ) * zt0
                zynpt = 0. - 2. * sin( to_rad*xlon[jj,ji] ) * zt0
                znnpt            = zxnpt*zxnpt + zynpt*zynpt
                #
                # j-direction: v-point segment direction (around t-point)
                zt1 = tan( pio4 - to_rad*xlat[jj+1,ji]/2. )
                zt2 = tan( pio4 - to_rad*xlat[jj-1,ji]/2. )
                zxvvt =  2.*( cos(to_rad*xlon[jj+1,ji])*zt1 - cos(to_rad*xlon[jj-1,ji])*zt2 )
                zyvvt =  2.*( sin(to_rad*xlon[jj+1,ji])*zt1 - sin(to_rad*xlon[jj-1,ji])*zt2 )
                znvvt = sqrt( znnpt * ( zxvvt*zxvvt + zyvvt*zyvvt )  )
                znvvt = max( znvvt, 1.e-12 )
                #
                sint = ( zxnpt*zyvvt - zynpt*zxvvt ) / znvvt
                cost = ( zxnpt*zxvvt + zynpt*zyvvt ) / znvvt
                angle[jj,ji] = atan2( sint, cost ) * 180./pi
        #
    return angle[:,:]




def RadiusEarth( lat ):
    '''
    Returns the radius of Earth in km as a function of latitude provided in degree N.
    '''
    latr = radians(lat) #converting into radians
    #R_eq    = 6378.137  #Radius at sea level at equator
    #R_pl    = 6356.752  #Radius at poles
    c    = (R_eq**2*cos(latr))**2
    d    = (R_pl**2*sin(latr))**2
    e    = (R_eq*cos(latr))**2
    f    = (R_pl*sin(latr))**2
    R    = sqrt((c+d)/(e+f))
    #print('\nRadius of Earth at '+str(round(lat,2))+'N = '+str(round(R,2))+' km\n')
    return R


def haversine_sclr( lat1, lon1, lat2, lon2 ):
    '''
    Returns the distance in km at the surface of the earth
    between two GPS points (degreesN, degreesE)
    '''
    R = RadiusEarth( 0.5*(lat1 + lat2) )
    dLat = radians(lat2 - lat1)
    dLon = radians(lon2 - lon1)
    lat1 = radians(lat1)
    lat2 = radians(lat2)
    a = sin(dLat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    c = 2*asin(sqrt(a))
    return R * c


def Haversine( plat, plon, xlat, xlon ):
    '''
    # Returns the distance in km at the surface of the earth
    # between two GPS points (degreesN, degreesE)
    # (plat,plon)  : a point
    # xlat, xlon : 2D arrays
    #
    # Here we do not need accuracy on Earth radius, since the call
    # to this function is suposely made to find nearest point
    '''
    to_rad = pi/180.
    #
    #R = RadiusEarth( plat ) ; # it's the target location that matters...
    R = 6360.
    #
    a1 = nmp.sin( 0.5 * ((xlat[:,:] - plat)*to_rad) )
    a2 = nmp.sin( 0.5 * ((xlon[:,:] - plon)*to_rad) )
    a3 = nmp.cos( xlat[:,:]*to_rad ) * cos(plat*to_rad)
    #
    return 2.*R*nmp.arcsin( nmp.sqrt( a1*a1 + a3 * a2*a2 ) )


def PlotMesh( pcoor_trg, Ys, Xs, isrc_msh, wghts, fig_name='mesh.png' ):
    '''
    isrc_msh: 2D integer array of shape (4,2)
    wghts:    1D real array of shape (4,)
    '''
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    #
    (yT,xT)                             = pcoor_trg
    [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = isrc_msh[:,:]
    [ wb1, wb2, wb3, wb4 ]              = wghts[:]
    #
    fig = plt.figure(num = 1, figsize=[7,5], facecolor='w', edgecolor='k')
    ax1 = plt.axes([0.09, 0.07, 0.6, 0.9])
    plt.plot( [  yT ] , [ xT ]  , marker='o', ms=15, color='k', label='target point' ) ; # target point !
    plt.plot( [ Xs[j1,i1] ], [ Ys[j1,i1] ], marker='o', ms=10, label='P1: w='+str(round(wb1,3)) ) ; # nearest point !
    plt.plot( [ Xs[j2,i2] ], [ Ys[j2,i2] ], marker='o', ms=10, label='P2: w='+str(round(wb2,3)) ) ; #
    plt.plot( [ Xs[j3,i3] ], [ Ys[j3,i3] ], marker='o', ms=10, label='P3: w='+str(round(wb3,3)) ) ; # target point !
    plt.plot( [ Xs[j4,i4] ], [ Ys[j4,i4] ], marker='o', ms=10, label='P4: w='+str(round(wb4,3)) ) ; # target point !
    ax1.legend(loc='center left', bbox_to_anchor=(1.07, 0.5), fancybox=True)
    plt.savefig(fig_name, dpi=100, transparent=False)
    plt.close(1)