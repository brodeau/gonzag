#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#       L. Brodeau, 2021
############################################################################

import sys
#from os import path, getcwd, mkdir
#import numpy as nmp

from math import radians, cos, sin, asin, sqrt, pi, tan, log, atan2, copysign
#from geopy import distance as geopy_distance

import numpy as nmp

from .utils  import Haversine,find_j_i_min


def Heading( plata,plona, platb,plonb ):
    '''
    #   ** Purpose: Compute true heading between point a and b
    #
    #   ** Method : suppose that the 2 points are not too far away from each other
    #               so that heading can be computed with loxodromy
    #
    #      "heading" is the angle (in degrees) that going from point A to point B makes:
    #                 IN A CLOCK-WISE FASHION !!!
    #               * 100% northward =>   0 deg.
    #               * 100% eastward  =>  90 deg.
    #               * 100% southward => 180 deg.
    #
    #  * history:
    #         J.M. Molines for SOSIE, may 2007
    #--------------------------------------------------------------
    '''
    repsilon = 1.E-9
    # There is a problem if the Greenwich meridian is between a and b
    xa = radians(plona)
    xb = radians(plonb)
    rr = max(abs(tan(pi/4.-radians(plata)/2.)), repsilon)  #lolo just to avoid FPE sometimes
    ya = -log(rr)
    rr = max(abs(tan(pi/4.-radians(platb)/2.)), repsilon)  #lolo just to avoid FPE sometimes
    yb = -log(rr)
    xb_xa=(xb-xa)%(2*pi)
    if  xb_xa >=  pi: xb_xa = xb_xa -2*pi
    if  xb_xa <= -pi: xb_xa = xb_xa +2*pi
    angled = atan2(xb_xa, yb - ya)
    rhdg = ( angled*180./pi ) % 360.
    #
    return rhdg


def AlfaBeta( vy, vx ):
    '''
    #----------------------------------------------------------
    #           ***  SUBROUTINE  local_coord    ***
    #
    #  ** Purpose : Compute the local coordinate in a grid cell
    #
    #  ** Method : from N. Daget Web page :
    #       http://aton.cerfacs.fr/~daget/TECHREPORT/TR_CMGC_06_18_html/node8.html
    #
    # * history:
    #      Original : J.M. Molines ( May 2007)
    #
    # INPUT:
    #         vy: vector of longitudes of length 5
    #         vx: vector of longitudes of length 5
    # RETURNS:
    #        alpha:
    #        beta:
    #----------------------------------------------------------
    '''
    nitermax = 100 ; # maximum number of iterations
    zresmax = 0.1
    zA = nmp.zeros((2,2))

    # when near the 0 deg line and we must work in the frame -180 180
    l_s_180 = ( abs(vx[1]-vx[4])>=180. or abs(vx[1]-vx[2])>=180. or abs(vx[1]-vx[3])>=180. )
    if l_s_180:
        zvx    = nmp.zeros(5)
        zvx[:] = vx[:]        ;  # back it up...
        vx[:]  = degE_to_degWE_1d(zvx[:])

    zres=1000. ; zdx=0.5 ; zdy=0.5 ; rA=0. ; rB=0. ; # Initialisation prior to convergence itterative loop
    jiter=0
    while (zres > zresmax) and (jiter < nitermax):
        z1 = vx[2] - vx[1]
        z2 = vx[1] - vx[4]
        z3 = vx[3] - vx[2]
        zA[0,0] =  z1 + (z2 + z3 )*rB
        zA[0,1] = -z2 + (z2 + z3 )*rA
        zA[1,0] = vy[2] - vy[1] + (vy[1] - vy[4] + vy[3] - vy[2])*rB
        zA[1,1] = vy[4] - vy[1] + (vy[1] - vy[4] + vy[3] - vy[2])*rA
        #
        # Determinant
        zdeta = nmp.linalg.det( zA )
        #
        # Solution of
        # |  zdx  |        | zdalp |
        # |       | =  zA .|       |
        # |  zdy  |        | zdbet |
        #zdeta = ( SIGN(1.,zdeta)*MAX(ABS(zdeta), repsilon) )  # just to avoid FPE division by zero sometimes...
        zM = nmp.zeros((2,2))
        zM[:,0] = [zdx,zdy]
        zM[:,1] = zA[:,1]
        zdalp = nmp.linalg.det( zM ) / zdeta
        zM[:,0] = zA[:,0]
        zM[:,1] = [zdx,zdy]
        zdbet = nmp.linalg.det( zM ) / zdeta
        # Update residual ( loop criteria)
        zres = sqrt( zdalp*zdalp + zdbet*zdbet )
        # Update alpha and beta from 1rst guess :
        rA = rA + zdalp
        rB = rB + zdbet
        # Update corresponding lon/lat for this alpha, beta
        z1a = 1.-rA
        z1b = 1.-rB
        zdx = vx[0] - (z1a*z1b*vx[1] + rA*z1b*vx[2] + rA*rB*vx[3] + z1a*rB*vx[4])
        zdy = vy[0] - (z1a*z1b*vy[1] + rA*z1b*vy[2] + rA*rB*vy[3] + z1a*rB*vy[4])
        #
        jiter = jiter + 1  # increment iteration counter
    # end of loop / until zres small enough (or nitermax reach )

    if  jiter >= nitermax:
        # Give up...
        rA = 0.5
        rB = 0.5
    # Problem if the 4 latitudes surrounding 'lati' are equal!
    if  vy[1]==vy[2] and vy[2]==vy[3] and vy[3]==vy[4]:
        rA  = 0.5
        rB  = 0.5
        ipb = 12
    if l_s_180:
        vx[:] = zvx[:] ; # back to values we got in input...
        del zvx
    del zA, zM
    #
    return rA, rB


def NearestPoint( pcoor_trg, Ys, Xs, rd_found_km=100., j_prv=0, i_prv=0, np_box_r=4 ):
    '''
    # * pcoor_trg : GPS coordinates (lat,lon) of target point    ([real],[real])
    # * Ys        : array of source grid latitude            2D numpy.array [real]
    # * Xs        : array of source grid longitude           2D numpy.array [real]
    '''
    (yT,xT) = pcoor_trg
    (Ny,Nx) = Ys.shape
    #
    j1=max(j_prv-np_box_r,0) ; j2=min(j_prv+np_box_r+1,Ny)
    i1=max(i_prv-np_box_r,0) ; i2=min(i_prv+np_box_r+1,Nx)
    lfound = False
    igo = 0
    while not lfound:
        igo = igo + 1
        if igo==2: (j1,i1 , j2,i2) = (0,0 , Ny,Nx) ; # Falling back on whole domain for second pass...
        xd = Haversine( yT, xT,  Ys[j1:j2,i1:i2], Xs[j1:j2,i1:i2] )
        jy, jx = find_j_i_min( xd )
        lfound = ( xd[jy,jx] < rd_found_km )
        if igo==2 and not lfound:
            print('ERROR: NearestPoint() => Fuck up #1!\n ==> maybe your "dist_found" too small?\n')
            sys.exit(0)
    if igo==1: jy=jy+j1 ; jx=jx+i1 ; # Translate to indices in whole domain:
    if jy<0 or jy<0 or jy>=Ny or jx>=Nx: print('ERROR: NearestPoint() => Fuck up #2!'); sys.exit(0)
    return [jy,jx]


def surrounding_j_i( jP, iP, k_ew_per=1, Nx=0 ):
    iPm1 = iP-1
    iPp1 = iP+1
    if k_ew_per>=0:
        if Nx < 1:
            print('ERROR / "surrounding_j_i@sosie_bilin.py": valid "Nx" must be provided when east-west periodicity!')
            sys.exit(0)
    if (iPm1 == -1):  iPm1 = Nx - k_ew_per - 1
    if (iPp1 == Nx):  iPp1 =      k_ew_per
    jPm1 = jP-1
    jPp1 = jP+1
    return jPm1, jPp1, iPm1, iPp1


def Iquadran2SrcMesh( jP, iP, iqd,  k_ew_per=1, Nx=0 ):
    '''
    # Given the value of iquadran 'iqd', this function returns the
    # 3 j,i locations of the 3 source-grid points that, together with
    # the nearest point (jP,iP), surround the target point.
    # In other words (jP,iP) and the 3 other j,i locations returned here
    # form the mesh in which the target point should be found!
    #
    # Input:
    #  * jP,iP    : indices of nearest point on source grid         [integer]
    #  * iqd      : "iquadran" position ( 1 <= iqd <= 4) [integer]
    #  * k_ew_per : east-wet periodicity of source grid: [integer]
    #               k_ew_per == -1 => none
    #               k_ew_per >=  0 => yes please! with an overlap of k_ew_per points
    #  * Nx       : x-size of the source grid, must be provide if k_ew_per >=  0
    '''
    #
    if not iqd in [1,2,3,4]:
        print('ERROR: do not call sosie_bilin.Iquadran2SrcMesh() wit a non valid "iquadran"!')
        sys.exit(0)
    #
    jPm1, jPp1, iPm1, iPp1 = surrounding_j_i( jP, iP,  k_ew_per=k_ew_per, Nx=Nx )
    #
    if   iqd==1:
        # nearest point is the bottom left corner point of local mesh
        i2=iPp1 ; j2 = jP
        i3=iPp1 ; j3 = jPp1
        i4=iP   ; j4 = jPp1
    elif iqd==2:
        # nearest point is the top left corner point of mesh
        i2=iP   ; j2 = jPm1
        i3=iPp1 ; j3 = jPm1
        i4=iPp1 ; j4 = jP
    elif iqd==3:
        # nearest point is the top righ corner point of mesh
        i2=iPm1 ; j2 = jP
        i3=iPm1 ; j3 = jPm1
        i4=iP   ; j4 = jPm1
    elif iqd==4:
        # nearest point is the bottom right corner point of mesh
        i2=iP   ; j2 = jPp1
        i3=iPm1 ; j3 = jPp1
        i4=iPm1 ; j4 = jP
    #
    return (j2,i2), (j3,i3), (j4,i4)



def Iquadran( pcoor_trg, Ys, Xs, jP, iP, k_ew_per=1, grid_s_angle=0. ):
    '''
    # "iquadran" determines which of the 4 "corners" of the source-grid mesh
    # (in wich the target point is located) is the nearest point.
    #  => takes value from 1 to 4, resp. adjacent cells NE, SE, SW, NW on the grid
    #     from a nearest-point point-of-view:
    #
    #      o--o           x--o            o--x            o--o
    # 1 => |  | NE   2 => |  | SE    3 => |  | SW    4 => |  | NW
    #      x--o           o--o            o--o            o--x
    #
    # * pcoor_trg : coordinates (lat,lon) of target point    ([real],[real])
    # * Ys        : array of source grid latitude            2D numpy.array [real]
    # * Xs        : array of source grid longitude           2D numpy.array [real]
    # * jP,iP     : indices of nearest point on source grid  [integer],[integer]
    # * k_ew_per  : east-wet periodicity of source grid:     [integer]
    #                 k_ew_per == -1 => none
    #                 k_ew_per >=  0 => yes please! with an overlap of k_ew_per points
    # * grid_s_angle: local distortion of the source grid = angle [degrees]
    '''
    if ldebug:
        from shapely.geometry import Point
        from shapely.geometry.polygon import Polygon
    #
    (yT,xT) = pcoor_trg
    (Ny,Nx) = Ys.shape
    #
    jPm1, jPp1, iPm1, iPp1 = surrounding_j_i( jP, iP,  k_ew_per=k_ew_per, Nx=Nx )
    iquadran               = -1 ; # "not found" flag value...
    #
    if (iPm1 < 0) or (jPm1 < 0) or (iPp1 > Nx-1):
        print('WARNING: sosie_bilin.Iquadran() => bound problem => ',xT,yT,Nx,Ny,iP,jP)
        print('          iPm1, iPp1, Nx =', iPm1, iPp1, Nx)
        print('          jPm1, jPp1, Ny =', jPm1, jPp1, Ny)
        print('         => ignoring current nearest point for i,j =', ji, jj, '(of target domain)\n')
        #
    else:
        #
        # Getting the 5 point of interest of the source grid:
        lon0 = Xs[jP  ,iP  ]%360. ; lat0 = Ys[jP  ,iP  ] # nearest point
        #
        # Restore target point longitude between 0 and 360
        zT = xT % 360. ; # lolo: zT don't wanna modify: xT
        #
        # The following HEADING stuf is aimed at finding in which source grid cell
        # the target point is comprised with regards to the nearest point
        # (there is actually 4 possible adjacent cells NE, SE, SW and NW)
        #
        # Compute heading of target point and neighbours from the nearest point
        #     => angles in degrees !!!
        ht = Heading(lat0,lon0,   yT,  zT)                          # target point
        hN = Heading(lat0,lon0, Ys[jPp1,iP  ], Xs[jPp1,iP  ]%360.)  # point above
        hE = Heading(lat0,lon0, Ys[jP  ,iPp1], Xs[jP  ,iPp1]%360.)  # point rhs
        hS = Heading(lat0,lon0, Ys[jPm1,iP  ], Xs[jPm1,iP  ]%360.)  # point below
        hW = Heading(lat0,lon0, Ys[jP  ,iPm1], Xs[jP  ,iPm1]%360.)  # 'west' on the grid
        #
        hz = degE_to_degWE(ht)  ;  # same as ht but in the [-180,+180] frame !
        if hN > hE: hN = hN -360.
        #
        if hz > hN and hz <= hE: iquadran=1
        if ht > hE and ht <= hS: iquadran=2
        if ht > hS and ht <= hW: iquadran=3
        if ht > hW and hz <= hN: iquadran=4
    #
    # The above method generally works fine when the grid is not too twisted
    # however whent it's not the case, we need to guess problematic regions
    # here we chose to look at the local rotation of the grid at the nearest point
    if (not iquadran in [1,2,3,4]) or (abs(grid_s_angle)>45.) :
        if ldebug:
            print(' *** / sosie_bilin.Iquadran(): need to go for heavy-duty mode in search for "iquadran"!')
            creason = 'first attempt did not succeed'
            if abs(grid_s_angle)>45.:
                creason = 'local source grid distortion > 45 degrees! (angle = '+str(grid_s_angle)+')'
            print('      ==> because '+creason)
        if not ldebug:
            from shapely.geometry import Point
            from shapely.geometry.polygon import Polygon
        # Using the heavy methods: testing all 4 possibilities, until point found inside the mesh polygon!
        for iq in [4,3,2,1]:
            (j2,i2), (j3,i3), (j4,i4) = Iquadran2SrcMesh( jP, iP, iq,  k_ew_per=k_ew_per, Nx=Nx )
            point = Point(yT,xT)
            polygon = Polygon([(Ys[jP,iP],Xs[jP,iP]), (Ys[j2,i2],Xs[j2,i2]), (Ys[j3,i3],Xs[j3,i3]), (Ys[j4,i4],Xs[j4,i4])])
            if polygon.contains(point):
                if ldebug: print(' ===> FOUND! iq=',iq,' is the guy!')
                iquadran = iq
                break
    #
    # Checking the point is inside the polygon pointed by iquadran:
    if ldebug and (iquadran in [1,2,3,4]):
        # lilo: test if inside polygon:
        (j2,i2), (j3,i3), (j4,i4) = Iquadran2SrcMesh( jP, iP, iquadran,  k_ew_per=k_ew_per, Nx=Nx )
        point = Point(yT,xT)
        polygon = Polygon([(Ys[jP,iP],Xs[jP,iP]), (Ys[j2,i2],Xs[j2,i2]), (Ys[j3,i3],Xs[j3,i3]), (Ys[j4,i4],Xs[j4,i4])])
        if not polygon.contains(point):
            print('WARNING / sosie_bilin.Iquadran(): point not inside polygon!')
            print(' * Model: jP, iP =', jP, iP, ' GPS => ', Ys[jP,iP], Xs[jP,iP])
            print('          local distortion of the grid =>', grid_s_angle,' degrees')
            print(' * Sat:  yT, zT =', yT, zT)
            print('Headings:')
            print(' * direction target:      ht =', ht, '    hz =',hz)
            print(' * direction point above: hN =', hN)
            print(' * direction point rhs:   hE =', hE)
            print(' * direction point below: hS =', hS)
            print(' * direction point lhs:   hW =', hW)
            print('   iquadran to be used =', iquadran,'\n')
    #
    # Final check:
    #if not iquadran in [1,2,3,4]:
    #    # Something is wrong!
    #    print('ERROR / sosie_bilin.Iquadran(): could not find iquadran!'); # => falling back on 0.5,0.5!')
    #    print(' * Model:jP, iP =', jP, iP, ' GPS => ', Ys[jP,iP], Xs[jP,iP])
    #    print(' * Sat:  yT, zT =', yT, zT)
    #    print('Headings:')
    #    print(' * direction target:      ht =', ht, '    hz =',hz)
    #    print(' * direction point above: hN =', hN)
    #    print(' * direction point rhs:   hE =', hE)
    #    print(' * direction point below: hS =', hS)
    #    print(' * direction point lhs:   hW =', hW)
    #    sys.exit(0)
    #
    return iquadran



def IDSourceMesh( pcoor_trg, Ys, Xs, jP, iP, k_ew_per=-1, grid_s_angle=0. ):
    '''
    # * jP,iP    : nearest point on source grid
    # * k_ew_per : E-W periodicity of source grid
    # * pcoor_trg : coordinates (lat,lon) of target point    ([real],[real])
    # * Ys        : array of source grid latitude            2D numpy.array [real]
    # * Xs        : array of source grid longitude           2D numpy.array [real]
    # * jP,iP     : indices of nearest point on source grid  [integer],[integer]
    # * k_ew_per  : east-wet periodicity of source grid:     [integer]
    #                 k_ew_per == -1 => none
    #                 k_ew_per >=  0 => yes please! with an overlap of k_ew_per points
    # * grid_s_angle: local distortion of the source grid = angle [degrees]
    #
    '''
    (yT,xT) = pcoor_trg
    (Ny,Nx) = Ys.shape
    #
    iqdrn = Iquadran( (yT,xT), Ys, Xs, jP, iP, k_ew_per=k_ew_per, grid_s_angle=grid_s_angle )
    #
    if iqdrn in [1,2,3,4]:
        #
        (j1,i1) = (jP,iP)
        (j2,i2), (j3,i3), (j4,i4) = Iquadran2SrcMesh( jP, iP, iqdrn,  k_ew_per=k_ew_per, Nx=Nx )
        #
    else:
        (j1,i1), (j2,i2), (j3,i3), (j4,i4) = (-1,-1), (-1,-1), (-1,-1), (-1,-1) ; # fuck-up flag...
    #
    return nmp.array([ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ], dtype=nmp.int64)



def WeightBL( pcoor_trg, Ys, Xs, isrc_msh ):
    '''
    # * pcoor_trg : coordinates (lat,lon) of target point    ([real],[real])
    # * Ys        : array of source grid latitude            2D numpy.array [real]
    # * Xs        : array of source grid longitude           2D numpy.array [real]
    # * isrc_msh  : the 4 "j,i" coordinates of the source mesh (as found by "IDSourceMesh()")
    #               => 2D numpy.array [int64] of shape (4,2)
    '''
    (yT,xT)                             = pcoor_trg
    [ [j1,i1],[j2,i2],[j3,i3],[j4,i4] ] = isrc_msh[:,:]
    #
    if j1 >= 0:
        vcoor = nmp.zeros((2,5))
        vcoor[:,0] = [    yT    ,    xT    ] ; # 0 => target point
        vcoor[:,1] = [Ys[j1,i1] , Xs[j1,i1]] ; # 1 => nearest point (source grid)
        vcoor[:,2] = [Ys[j2,i2] , Xs[j2,i2]]
        vcoor[:,3] = [Ys[j3,i3] , Xs[j3,i3]]
        vcoor[:,4] = [Ys[j4,i4] , Xs[j4,i4]]
        #
        alfa, beta = AlfaBeta( vcoor[0,:], vcoor[1,:] )
        del vcoor
        #
        # weights for interpolation
        w1 = (1. - alfa)*(1. - beta)
        w2 =       alfa *(1. - beta)
        w3 =       alfa * beta
        w4 = (1. - alfa)* beta
        #
    else:
        # source mesh (iquadran in the first place) could not be identified...
        w1 = 1. ; w2 = 0. ; w3 = 0. ; w4 = 0. ; # => nearest-point interpolation...
        if ldebug:
            print('   WARNING / sosie_bilin.WeightBL(): will use nearest point value')
            print('             as source mesh could not be identidied!')
            print('    * Model:jP, iP =', j1, i1, ' GPS => ', Ys[j1,i1], Xs[j1,i1])
            print('    * Sat:  yT, xT =', yT, xT)
    #
    return nmp.array([w1, w2, w3, w4])


class BilinTrack:
    ''' 
    '''    
    def __init__( self, Yt, Xt,  Ys, Xs, src_grid_local_angle=[], \
                  k_ew_per=-1, rd_found_km=100., np_box_r=4, freq_talk=0 ):
        #
        (self.Nt,)  = nmp.shape(Yt)
        self.Yt     = Yt
        self.Xt     = Xt
        self.Ys     = Ys
        self.Xs     = Xs
        self.sangle = src_grid_local_angle
        self.kewp   = k_ew_per
        self.rfound = rd_found_km
        self.nprad  = np_box_r

        self.NP = nmp.zeros((self.Nt,2)  , dtype=nmp.int64) ; # nearest point
        self.SM = nmp.zeros((self.Nt,4,2), dtype=nmp.int64) ; # source mesh
        self.WB = nmp.zeros((self.Nt,4))                    ; # weights

        self.ftalk = 2*self.Nt ; # => will never talk!
        if freq_talk>0: self.ftalk = int(freq_talk) ; # will talk every "ftalk" increments....
        
        print('\n *** Finding nearest points on source (model) grid... (rd_found_km, np_box_r =',rd_found_km, np_box_r,')')
        self.NP = self.nrpt( )
        print('     ***    Done! *** \n')

        print('  *** Determining source meshes...')
        self.SM = self.srcm( )
        print('     ***    Done! *** \n')
        
        print('  *** Computing bilinear weights...')
        self.WB = self.wght( )
        print('     ***    Done! *** \n')
        
    def nrpt( self ):
        #
        xnp = nmp.zeros((self.Nt,2), dtype=nmp.int64)
        #
        [jj,ji] = [self.nprad,self.nprad] ; # stupid first guess here...
        #
        for jt in range(self.Nt):
            ltalk = ((jt+1)%self.ftalk==0)

            if ltalk: print('      +++ Treated point: '+str(jt+1)+'/'+str(self.Nt), \
                            '\n          ==> Sat. coordinates:    ', round(self.Yt[jt],3), round(self.Xt[jt],3))
            
            [jj,ji] = NearestPoint( (self.Yt[jt],self.Xt[jt]), self.Ys, self.Xs, \
                                    rd_found_km=self.rfound, j_prv=jj, i_prv=ji, np_box_r=self.nprad )
            xnp[jt,:] = [jj,ji]

            if ltalk: print('          ==> Model nearest point: ', \
                            round(self.Ys[jj,ji],3),round(self.Xs[jj,ji]%360.,3),' (',jj,ji,')')
        #
        return xnp
    
    def srcm( self ):
        # 
        xsp = nmp.zeros((self.Nt,4,2), dtype=nmp.int64)

        for jt in range(self.Nt):
            [jP,iP] = self.NP[jt,:]            
            angle = 0.
            if nmp.shape(self.sangle) == nmp.shape(self.Ys): angle = self.sangle[jP,iP]
            #
            xsp[jt,:,:] = IDSourceMesh( (self.Yt[jt],self.Xt[jt]), self.Ys, self.Xs, jP, iP, \
                                        k_ew_per=self.kewp, grid_s_angle=angle )
        #
        return xsp

    def wght( self ):
        # 
        x4w = nmp.zeros((self.Nt,4))
        #
        for jt in range(self.Nt):
            x4w[jt,:] = WeightBL( (self.Yt[jt],self.Xt[jt]), self.Ys, self.Xs, self.SM[jt,:,:] )
        #
        return x4w
