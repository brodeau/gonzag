#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#   ///// https://github.com/brodeau/gonzag \\\\\
#
#       L. Brodeau, 2021
#
############################################################################

from sys import exit
from math import radians, cos, sin, asin, sqrt, pi, tan, log, atan2, copysign
import numpy as nmp
from .config import IsZarr, ldebug, R_eq, R_pl, deg2km

def MsgExit( cmsg ):
    print('\n ERROR: '+cmsg+' !\n')
    exit(0)

def chck4f( ncfile ):
    from os.path import exists
    if not exists(ncfile):
        MsgExit('File '+ncfile+' does not exist')

def degE_to_degWE( X ):
    '''
    # From longitude in 0 -- 360 frame to -180 -- +180 frame...
    '''
    if nmp.shape( X ) == ():
        # X is a scalar
        return     copysign(1., 180.-X)*        min(X,     abs(X-360.))
    else:
        # X is an array
        return nmp.copysign(1., 180.-X)*nmp.minimum(X, nmp.abs(X-360.))
    

def EpochT2Str( time ):
    '''
    # Input: UNIX epoch time (integer or float)
    # Returns: a string of the date understandable by mamals...
    '''
    from datetime import datetime as dtm
    #
    itime = int(round(time,0))
    #
    return dtm.utcfromtimestamp(itime).strftime('%c')


def find_j_i_min(x):
    '''
    # Yes, reinventing the wheel here, but it turns out
    # it is faster this way!
    '''
    k = x.argmin()
    nx = x.shape[1]
    return k//nx, k%nx

def find_j_i_max(x):
    '''
    # Yes, reinventing the wheel here, but it turns out
    # it is faster this way!
    '''
    k = x.argmax()
    nx = x.shape[1]
    return k//nx, k%nx



def GridResolution( X ):
    '''
    # X    : the 2D array of the model grid longitude
    '''
    ny = X.shape[0]
    vx = nmp.abs(X[ny//2,1:] - X[ny//2,:-1])
    res = nmp.mean( vx[nmp.where(vx < 120.)] )
    if ldebug: print(' *** [GridResolution()] Based on the longitude array, the model resolution ~= ', res, ' degrees \n')
    return res


def IsEastWestPeriodic( X ):
    '''
    # X    : the 2D array of the model grid longitude [0:360]
    #  RETURNS: iper: -1 => no E-W periodicity ; iper>=0 => E-W periodicity with iper overlaping points!
    '''
    ny = X.shape[0]
    jj = ny//2 ; # we test at the center...
    iper = -1
    dx = X[jj,1] - X[jj,0]
    lon_last_p1 = (X[jj,-1]+dx)%360.
    for it in range(5):
        if lon_last_p1 == X[jj,it]%360.:
            iper = it
            break
    return iper



def SearchBoxSize( res_mod, width_box ):
    '''
    # Returns half of the width, in number of grid points, of the small zoom-box 
    # of the source (model) domain in which NearestPoint() will initially attempt
    # to locate  the nearest point, before falling back on the whole source (model)
    # domain if unsuccessful.                           
    #  => the smaller the faster the search...
    #
    # * res_mod:   horizontal resolution of the model data, in km
    # * width_box: width of the zoom-box in km
    #
    # TODO: shoud take into account the speed of the satellite
    '''
    return int(0.5*width_box/res_mod)

    

#def IsGlobalLongitudeWise( X, resd=1. ):
#    '''
#    # X    : the 2D array of the model grid longitude
#    # resd : rough order of magnitude of the resolutio of the model grid in degrees
#    # RETURNS: boolean + lon_min and lon_max in the [-180:+180] frame
#    '''
#    X = nmp.mod(X, 360.) ; # no concern, it should already have been done earlier anyway...
#    print( 'X =', X)
#    print( 'X =', degE_to_degWE(X))
#    xmin1 = nmp.amin(degE_to_degWE(X)) ; # in [-180:+180] frame...
#    xmax1 = nmp.amax(degE_to_degWE(X)) ; #     "      "
#    xmin2 = nmp.amin(X) ; # in [0:360] frame...
#    xmax2 = nmp.amax(X) ; #     "     "
#    print(' xmin2, xmax2 =', xmin2, xmax2 )
#    print(' xmin1, xmax1 =', xmin1, xmax1 )
#    l1 = ( xmin1<0. and xmin1>5.*resd-180 ) or ( xmax1>0. and xmax1<180.-5.*resd )
#    l2 = ( xmin2<1.5*resd and xmax2>360.-1.5*resd )
#    #print('l1, l2 =', l1, l2)
#    if not l1:
#        if l2 :
#            lglobal = True
#            print(' *** From model longitude => looks like global setup (longitude-wise)', xmin1, xmax1,'\n')
#        else:
#            MsgExit('cannot find if regional or global source domain...')
##    else:
#        print(' *** From model longitude => looks like regional setup (longitude-wise)', xmin1, xmax1)
#        lglobal = False
#        print('     => will disregard alongtrack points with lon < '+str(xmin1)+' and lon > '+str(xmax1),'\n')
#    #
#    return lglobal, xmin1, xmax1


def IsGlobalLongitudeWise( X, resd=1. ):
    '''
    # LIMITATION: longitude has to increase in the x-direction (second dimension) of X (i increases => lon increases)
    # X    : the 2D array of the model grid longitude
    # resd : rough order of magnitude of the resolutio of the model grid in degrees
    # RETURNS: boolean, boolean, lon_min, lon_max
    '''
    #
    nx   = X.shape[1]
    X    = nmp.mod(X, 360.) ; # no concern, it should already have been done earlier anyway...
    xmin = nmp.amin(X) ; # in [0:360] frame...
    xmax = nmp.amax(X) ; #     "     "
    imin = nmp.argmin(X)%nx
    imax = nmp.argmax(X)%nx
    #
    xb = degE_to_degWE(X)
    xminB = nmp.amin(xb) ; # in [-180:+180] frame...
    xmaxB = nmp.amax(xb) ; #     "      "    
    #
    l360    = True   ; # we'll be in the [0:360] frame...
    lglobal = False
    #
    if (xmin<1.5*resd and xmax>360.-1.5*resd) and (imax > imin):
        # Global longitude domain
        lglobal = True
        xmin=0. ; xmax=360.
        #
    elif (xminB%360. > xmaxB%360.) and (imax < imin):
        # (xminB%360. > xmaxB%360.) is True in the 2 singular cases: domain icludes Greenwhich Meridian or -180:180 transition...
        l360 = False
        xmin = xminB
        xmax = xmaxB
    #
    del xb
    return lglobal, l360, xmin, xmax




def GetEpochTimeOverlap( ncfile_sat, ncfile_mod ):
    '''
    # * irange_sat: time coverage in epoch Unix time for satellite data: (int,int)
    # * irange_sat: time coverage in epoch Unix time for model     data: (int,int)
    '''
    if not IsZarr:
        from .ncio import GetTimeInfo
    else:
        from .zarrio import GetTimeInfo
    #
    nts, range_sat = GetTimeInfo( ncfile_sat )
    ntm, range_mod = GetTimeInfo( ncfile_mod )
    #
    (zt1_s,zt2_s) = range_sat
    (zt1_m,zt2_m) = range_mod
    if ldebug: print('\n *** [GetEpochTimeOverlap()] Earliest/latest dates:\n   => for satellite data:',zt1_s,zt2_s,'\n   => for model     data:',zt1_m,zt2_m,'\n')
    if (zt1_m >= zt2_s) or (zt1_s >= zt2_m) or (zt2_m <= zt1_s) or (zt2_s <= zt1_m):
        MsgExit('No time overlap for Model and Track file')
    return (max(zt1_s, zt1_m), min(zt2_s, zt2_m)), (nts, ntm)


def scan_idx( vt, rt1, rt2 ):
    '''
    # Finding indices when we can start and stop when scanning the track file:
    # * vt: vector containing dates as Epoch UNIX time [float]
    # * rt1, rt2: the 2 dates of interest (first and last) [float]
    # RETURNS: the two corresponding position indices 
    '''
    nt = len(vt)
    for kt1 in range(  0, nt-1):
        if  (vt[kt1] <= rt1) and (vt[kt1+1] > rt1): break
    for kt2 in range(kt1, nt-1):
         if (vt[kt2] <= rt2) and (vt[kt2+1] > rt2): break
    kt2 = kt2 + 1
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
    # ! VECTOR VERSION !
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
    a1 = nmp.sin( 0.5 * ((xlat - plat)*to_rad) )
    a2 = nmp.sin( 0.5 * ((xlon - plon)*to_rad) )
    a3 = nmp.cos( xlat*to_rad ) * cos(plat*to_rad)
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
    plt.plot( [     xT    ], [     yT    ], marker='o', ms=15, color='k', label='target point' ) ; # target point !
    plt.plot( [ Xs[j1,i1] ], [ Ys[j1,i1] ], marker='o', ms=10, label='P1: w='+str(round(wb1,3)) ) ; # nearest point !
    plt.plot( [ Xs[j2,i2] ], [ Ys[j2,i2] ], marker='o', ms=10, label='P2: w='+str(round(wb2,3)) ) ; #
    plt.plot( [ Xs[j3,i3] ], [ Ys[j3,i3] ], marker='o', ms=10, label='P3: w='+str(round(wb3,3)) ) ; # target point !
    plt.plot( [ Xs[j4,i4] ], [ Ys[j4,i4] ], marker='o', ms=10, label='P4: w='+str(round(wb4,3)) ) ; # target point !
    ax1.legend(loc='center left', bbox_to_anchor=(1.07, 0.5), fancybox=True)
    plt.savefig(fig_name, dpi=100, transparent=False)
    plt.close(1)




# C L A S S E S


    
class ModGrid:
    '''
    # Will provide: size=nt, shape=(ny,nx), time[:], lat[:,:], lon[:,:] of Model data...
    # mask
    # domain_bounds (= [ lat_min, lon_min, lat_max, lon_max ])
    '''
    def __init__( self, ncfile, rtu1, rtu2, nclsm, varlsm, distorded_grid=False ):
        '''
        # * ncfile: netCDF file containing satellite track
        # * rtu1: Epoch UNIX time to start getting time from (included) [float]
        # * rtu2: Epoch UNIX time to stop  getting time from (included) [float]
        # * nclsm, varlsm: file and variable to get land-sea mask...
        '''
        if not IsZarr:
            from .ncio import GetTimeEpochVector, GetModelCoor, GetModelLSM, Save2Dfield
        else:
            from .zarrio import GetTimeEpochVector, GetModelCoor, GetModelLSM, Save2Dfield

        chck4f( ncfile )

        self.file = ncfile
        
        rvt = GetTimeEpochVector( ncfile )
        jt1, jt2 = scan_idx( rvt, rtu1, rtu2 )
        self.size = jt2 - jt1 + 1        
        self.time = GetTimeEpochVector( ncfile, kt1=jt1, kt2=jt2 )

        
        zlat =          GetModelCoor( ncfile, 'latitude' )
        zlon = nmp.mod( GetModelCoor( ncfile, 'longitude') , 360. )
        if len(zlat.shape)==1 and len(zlon.shape)==1:
            # Must build the 2D version:
            print(' *** Model latitude and longitude arrays are 1D => building the 2D version!')
            ny = len(zlat) ;  nx = len(zlon)
            xlat = nmp.zeros((ny,nx))
            xlon = nmp.zeros((ny,nx))
            for jx in range(nx): xlat[:,jx] = zlat[:]
            for jy in range(ny): xlon[jy,:] = zlon[:]
            self.lat = xlat
            self.lon = xlon
            del xlat, xlon
        else:
            # Already 2D:
            if zlat.shape != zlon.shape: MsgExit('[SatTrack()] => lat and lon disagree in shape')
            self.lat = zlat
            self.lon = zlon            
        del zlat, zlon
        
        self.shape = self.lat.shape

        # Land-sea mask
        chck4f( nclsm )
        self.mask = GetModelLSM( nclsm, varlsm ) ; 
        if self.mask.shape != self.shape: MsgExit('model land-sea mask has a wrong shape => '+str(self.mask.shape))
        if ldebug: Save2Dfield( 'mask_model.nc', self.mask, name='mask' )
        
        # Horizontal resolution
        self.HResDeg = GridResolution( self.lon )
        if self.HResDeg>5. or self.HResDeg<0.001: MsgExit('Model resolution found is surprising, prefer to stop => check "GetModelResolution()" in utils.py')        
        self.HResKM  = self.HResDeg*deg2km

        # Globality and East-West periodicity ?
        self.IsLonGlobal, self.l360, lon_min, lon_max = IsGlobalLongitudeWise( self.lon, resd=self.HResDeg )
        if self.IsLonGlobal:
            self.EWPer = IsEastWestPeriodic( self.lon )
        else:
            self.EWPer = -1

        lat_min = nmp.amin(self.lat)
        lat_max = nmp.amax(self.lat)

        # Local distortion angle:
        self.IsDistorded = distorded_grid
        if distorded_grid:
            print(' *** Computing angle distortion of the model grid ("-D" option invoked)...')
            self.xangle = GridAngle( self.lat, self.lon )
            if ldebug: Save2Dfield( 'model_grid_disortion.nc', self.xangle, name='angle', mask=self.mask )            
        else:
            print(' *** Skipping computation of angle distortion of the model grid! ("-D" option not invoked)...')
            self.xangle = nmp.zeros(self.shape)    
        
        self.domain_bounds = [ lat_min, lon_min, lat_max, lon_max ]

        if not self.l360: self.lon = degE_to_degWE( self.lon )
                
        # Summary:
        print('\n *** About model gridded (source) domain:')
        print('     * shape = ',self.shape)
        print('     * horizontal resolution: ',self.HResDeg,' degrees or ',self.HResKM,' km')
        print('     * Is this a global domain w.r.t longitude: ', self.IsLonGlobal)
        if self.IsLonGlobal:
            print('       ==> East West periodicity: ', (self.EWPer>=0), ', with an overlap of ',self.EWPer,' points')
        else:
            print('       ==> this is a regional domain')
            if self.l360:
                print('       ==> working in the [0:360] frame...')
            else:
                print('       ==> working in the [-180:180] frame...')
        print('     * lon_min, lon_max = ', round(lon_min,2), round(lon_max,2))
        print('     * lat_min, lat_max = ', round(lat_min,2), round(lat_max,2))
        print('     * should we pay attention to possible STRONG local distorsion in the grid: ', self.IsDistorded)
        print('     * number of time records of interest for the interpolation to come: ', self.size)
        print('       ==> time record indices: '+str(jt1)+' to '+str(jt2)+', included\n')









class SatTrack:
    '''
    # Will provide: size, time[:], lat[:], lon[:] of Satellite track
    '''
    def __init__( self, ncfile, rtu1, rtu2, Np=0, domain_bounds=[-90.,0. , 90.,360.], l_0_360=True ):
        '''
        # *  ncfile: netCDF file containing satellite track
        # *  rtu1: Epoch UNIX time to start getting time from (included) [float]
        # *  rtu2: Epoch UNIX time to stop  getting time from (included) [float]
        # ** Np:     number of points (size) of track in netCDF file...
        # ** domain_bounds: bound of region we are interested in => [ lat_min, lon_min, lat_max, lon_max ]
        '''
        if not IsZarr:
            from .ncio import GetTimeEpochVector, GetSatCoor
        else:
            from .zarrio import GetTimeEpochVector, GetSatCoor

        chck4f( ncfile )

        self.file = ncfile

        print(' *** [SatTrack()] Analyzing the time vector in '+ncfile+' ...')

        #import time
        #startTime = time.time()
        
        if Np<2500:
            # Can afford to read whole time vector, not a problem with such as small of number of records to read
            rvt = GetTimeEpochVector( ncfile, lquiet=True )
            jt1, jt2 = scan_idx( rvt, rtu1, rtu2 )
        else:
            # Subsampling with increment of 500 for first pass...
            kss = 500
            rvt = GetTimeEpochVector( ncfile, isubsamp=kss, lquiet=True ) ; # WAY faster to read a subsampled array with NetCDF-4 !
            j1, j2 = scan_idx( rvt, rtu1, rtu2 )
            j1 = j1*kss ; j2 = j2*kss ; # shorter new range in which to search
            rvt = GetTimeEpochVector( ncfile, kt1=j1, kt2=j2, lquiet=True ) ; # reading without subsampling but a shorter slice now!
            jt1, jt2 = scan_idx( rvt, rtu1, rtu2 )
            jt1 = jt1+j1 ; jt2 = jt2+j1 ; # convert in term of whole length
            del j1, j2, kss
        
        nt = jt2 - jt1 + 1        

        self.jt1   = jt1
        self.jt2   = jt2

        vtime = GetTimeEpochVector( ncfile, kt1=jt1, kt2=jt2 )

        #time_read_time = time.time() - startTime
        #print(' *** ROOM FOR IMPROVEMENT => time_read_time =', time_read_time,'\n')

        vlat  =        GetSatCoor( ncfile, 'latitude' , jt1,jt2 )
        vlon  =        GetSatCoor( ncfile, 'longitude', jt1,jt2 )

        # Make sure we are in the same convention as model data
        # (model data can be in [-180:180] range if regional domain that crosses Greenwhich meridian...
        if l_0_360:
            vlon = nmp.mod( vlon, 360. )
        else:
            vlon = degE_to_degWE( vlon )
            
        #print(' lolo: track size before removing points outside of model domain: '+str(len(vtime)))
        [ ymin,xmin , ymax,xmax ] = domain_bounds
        #print('lolo: lat_min, lat_max =', ymin, ymax)
        #print('lolo: lon_min, lon_max =', xmin, xmax)
        keepit = nmp.where( (vlat[:]>=ymin) & (vlat[:]<=ymax) & (vlon[:]>=xmin) & (vlon[:]<=xmax) )
        #print(' lolo: keepit =', keepit )

        self.time  = vtime[keepit]
        self.lat   =  vlat[keepit]
        self.lon   =  vlon[keepit]
        
        self.size = len(self.time)
        self.keepit = keepit

        del vtime, vlat, vlon
        #print(' lolo: track size AFTER removing points outside of model domain: '+str(self.size))
        
        print('\n *** About satellite track (target) domain:')
        print('     * number of time records of interest for the interpolation to come: ', self.size)
        print('       ==> time record indices: '+str(jt1)+' to '+str(jt2)+', included\n')
        
