#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

import numpy as nmp
from netCDF4 import Dataset, num2date, default_fillvals
from calendar import timegm
from datetime import datetime as dtm
from .config import ldebug, rmissval

cabout_nc = 'Created with Gonzag package => https://github.com/brodeau/gonzag'




def GetTimeInfo( ncfile ):
    '''
    # Inspect time dimension
    # Get number of time-records, first and last date
    # Return them + dates as UNIX epoch time, aka "seconds since 1970-01-01 00:00:00" (integer)
    '''
    if ldebug: print(' *** [GetTimeInfo()] Getting calendar/time info in '+ncfile+' ...')
    id_f = Dataset(ncfile)
    list_dim = id_f.dimensions.keys()
    list_var = id_f.variables.keys()
    for cd in [ 'time', 'time_counter', 'TIME', 'record', 't', 'none' ]:
        if cd in list_dim: break
    if cd == 'none': MsgExit('found no time-record dimension in file '+ncfile)
    if ldebug: print('   => time/record dimension is "'+cd+'"')
    nt = id_f.dimensions[cd].size
    cv = cd ; # ASSUMING VAR == DIM !!! Is it bad???
    if not cv in list_var: MsgExit('name of time variable is different than name of time dimension')
    clndr = id_f.variables[cv]
    dt1 = num2date( clndr[0],    clndr.units, clndr.calendar )
    dt2 = num2date( clndr[nt-1], clndr.units, clndr.calendar )
    id_f.close()
    #
    cfrmt = '%Y-%m-%d %H:%M:%S'
    it1 = timegm( dtm.strptime( dt1.strftime(cfrmt) , cfrmt ).timetuple() )
    it2 = timegm( dtm.strptime( dt2.strftime(cfrmt) , cfrmt ).timetuple() )
    if ldebug: print('   => first and last time records: ',dt1,'--',dt2,' (UNIX epoch: ', it1,'--',it2,')\n')
    #
    return nt, (dt1,dt2), (it1,it2)


def GetTimeEpochVector( ncfile, kt1=0, kt2=0, isubsamp=1, lquiet=False ):
    '''
    # Get the time vector in the netCDF file, (from index kt1 to kt2, if these 2 are != 0!)
    # returns it as:
    # isubsamp: subsampling !!!
    # lquiet: shut the f* up!
    #  => ivt: time as UNIX epoch time, aka "seconds since 1970-01-01 00:00:00" (integer)
    '''
    ltalk = ( ldebug and not lquiet )
    cv_t_test = [ 'time', 'time_counter', 'TIME', 'record', 't', 'none' ]
    id_f = Dataset(ncfile)
    list_var = id_f.variables.keys()
    for cv in cv_t_test:
        if cv in list_var:
            clndr = id_f.variables[cv]
            if ltalk: print(' *** [GetTimeEpochVector()] reading "'+cv+'" in '+ncfile+' and converting it to Epoch time...')
            if kt1>0 and kt2>0:
                if kt1>=kt2: MsgExit('mind the indices when calling GetTimeEpochVector()')
                vdate = num2date( clndr[kt1:kt2+1:isubsamp], clndr.units, clndr.calendar )
                cc = 'read...'
            else:
                vdate = num2date( clndr[::isubsamp],         clndr.units, clndr.calendar )
                cc = 'in TOTAL!'
            break
    id_f.close()
    if cv == 'none': MsgExit('found no time-record variable in file '+ncfile+' (possible fix: "cv_t_test" in "GetTimeEpochVector()")')
    #
    # Create the Unix Epoch time version:
    nt = len(vdate)
    ivt = nmp.zeros(nt, dtype=nmp.int64)
    cfrmt = '%Y-%m-%d %H:%M:%S'
    for jt in range(nt): ivt[jt] = timegm( dtm.strptime( vdate[jt].strftime(cfrmt) , cfrmt ).timetuple() )
    #
    if ltalk: print('   => '+str(nt)+' records '+cc+'\n')
    return ivt


def GetModelCoor( ncfile, what ):
    '''
    #   list_dim = list(id_f.dimensions.keys()) ;  print(" list dim:", list_dim)
    '''
    cv_coor_test = nmp.array([[ 'lat','latitude', 'nav_lat','gphit','LATITUDE', 'none' ],
                              [ 'lon','longitude','nav_lon','glamt','LONGITUDE','none' ]])
    if   what ==  'latitude': ii = 0
    elif what == 'longitude': ii = 1
    else: MsgExit(' "what" argument of "GetModelCoor()" only supports "latitude" and "longitude"')
    #
    id_f = Dataset(ncfile)
    list_var = list(id_f.variables.keys())
    for ncvar in cv_coor_test[ii,:]:
        if ncvar in list_var: break
    if ncvar == 'none': MsgExit('could not find '+what+' array into model file (possible fix: "cv_coor_test" in "GetModelCoor()")')
    #
    nb_dim = len(id_f.variables[ncvar].dimensions)
    if   nb_dim==1: MsgExit('FIX ME! Model '+what+' is 1D')
    elif nb_dim==2: xwhat = id_f.variables[ncvar][:,:]
    elif nb_dim==3: xwhat = id_f.variables[ncvar][0,:,:]
    else: MsgExit('FIX ME! Model '+what+' has a weird number of dimensions')
    id_f.close()
    if ldebug: print(' *** [GetModelCoor()] Read model '+what+' (variable is "'+ncvar+'", with '+str(nb_dim)+' dimensions!',nmp.shape(xwhat),'\n')
    #
    return xwhat


def GetModelLSM( ncfile, what ):
    '''
    # Returns the land-sea mask on the source/moded domain: "1" => ocean point, "0" => land point
    # => 2D array [integer]
    '''
    print('\n *** what we use to define model land-sea mask:\n    => "'+what+'" in "'+ncfile+'"\n')
    l_fill_val = (what[:10]=='_FillValue')
    ncvar = what
    if l_fill_val: ncvar = what[11:]
    #
    id_f = Dataset(ncfile)
    ndim = len(id_f.variables[ncvar].dimensions)
    if l_fill_val:
        # Mask is constructed out of variable and its missing value
        if ndim!=3: MsgExit(ncvar+' is expected to have 3 dimensions')
        xmsk = 1 - id_f.variables[ncvar][0,:,:].mask
    else:
        # Mask is read in mask file...
        if   ndim==2: xmsk = id_f.variables[ncvar][:,:]
        elif ndim==3: xmsk = id_f.variables[ncvar][0,:,:]
        elif ndim==4: xmsk = id_f.variables[ncvar][0,0,:,:]
        else: MsgExit('FIX ME! Mask '+ncvar+' has a weird number of dimensions:'+str(ndim))
    #
    id_f.close()
    return xmsk.astype(int)


def GetModel2DVar( ncfile, ncvar, kt=0 ):
    '''
    #   Fetches the 2D field "ncvar" at time record kt into "ncfile"
    '''
    if ldebug: print(' *** [GetModel2DVar()] Reading model "'+ncvar+'" at record kt='+str(kt)+' in '+ncfile)
    id_f = Dataset(ncfile)    
    x2d = id_f.variables[ncvar][kt,:,:]
    id_f.close()
    if ldebug: print('')
    return x2d




def GetSatCoord( ncfile, what,  kt1=0, kt2=0 ):
    '''
    # Get latitude (what=='latitude') OR longitude (what=='longitude') vector
    # in the netCDF file, (from index kt1 to kt2, if these 2 are != 0!)
    '''
    cv_coor_test = nmp.array([[ 'lat','latitude', 'LATITUDE',  'none' ],
                              [ 'lon','longitude','LONGITUDE', 'none' ]])
    if   what ==  'latitude': ii = 0
    elif what == 'longitude': ii = 1
    else: MsgExit('"what" argument of "GetSatCoord()" only supports "latitude" and "longitude"')
    #
    id_f = Dataset(ncfile)
    list_var = list(id_f.variables.keys())
    for ncvar in cv_coor_test[ii,:]:
        if ncvar in list_var: break
    if ncvar == 'none': MsgExit('could not find '+what+' array into satellite file (possible fix: "cv_coor_test" in "GetSatCoord()")')
    #
    if ldebug: print(' *** [GetSatCoord()] reading "'+ncvar+'" in '+ncfile+' ...')
    nb_dim = len(id_f.variables[ncvar].dimensions)
    if nb_dim==1:
        if kt1>0 and kt2>0:
            if kt1>=kt2: MsgExit('mind the indices when calling GetSatCoord()')
            vwhat = id_f.variables[ncvar][kt1:kt2+1]
            cc = 'read...'
        else:
            vwhat = id_f.variables[ncvar][:]
            cc = 'in TOTAL!'
    else:
        MsgExit('FIX ME! Satellite '+what+' has a weird number of dimensions (we expect only 1: the time-record!)')
    id_f.close()
    if ldebug: print('   => '+str(vwhat.size)+' records '+cc+'\n')
    #
    return vwhat


def GetSatSSH( ncfile, ncvar,  kt1=0, kt2=0, ikeep=[] ):
    '''
    # Get vector time-series of 'ncvar' in file 'ncfile'!
    #  - from index kt1 to kt2, if these 2 are != 0!
    #  - if (ikeep != []) => only retains parts of the data for which indices are provide into ikeep
    #          (ikeep is an array obtained as the result of a "numpy.where()"
    '''
    if ldebug: print(' *** [GetSatSSH()] Reading satellite "'+ncvar+' in '+ncfile)
    id_f = Dataset(ncfile)
    if kt1>0 and kt2>0:
        if kt1>=kt2: MsgExit('mind the indices when calling GetSatCoord()')
        vssh = id_f.variables[ncvar][kt1:kt2+1]
    else:
        vssh = id_f.variables[ncvar][:]
    id_f.close()
    if len(ikeep) > 0:
        # Keep specified part of the data
        vssh = vssh[ikeep]
    #
    if nmp.ma.is_masked(vssh): vssh[nmp.where( nmp.ma.getmask(vssh) )] = rmissval
    if ldebug: print('')
    return vssh




### OUTPUT :


def Save2Dfield( ncfile, XFLD, xlon=[], xlat=[], name='field', unit='', long_name='', mask=[], clon='nav_lon', clat='nav_lat' ):
    #LOLO IMPROVE !
    (nj,ni) = nmp.shape(XFLD)
    f_o = Dataset(ncfile, 'w', format='NETCDF4')
    f_o.createDimension('y', nj)
    f_o.createDimension('x', ni)
    if (xlon != []) and (xlat != []):
        if (xlon.shape == (nj,ni)) and (xlon.shape == xlat.shape):
            id_lon  = f_o.createVariable(clon ,'f4',('y','x',), zlib=True, complevel=5)
            id_lat  = f_o.createVariable(clat ,'f4',('y','x',), zlib=True, complevel=5)
            id_lon[:,:] = xlon[:,:]
            id_lat[:,:] = xlat[:,:]
    id_fld  = f_o.createVariable(name ,'f4',('y','x',), zlib=True, complevel=5)
    if long_name != '': id_fld.long_name = long_name
    if unit      != '': id_fld.units     = unit
    if nmp.shape(mask) != (0,):
        xtmp = nmp.zeros((nj,ni))
        xtmp[:,:] = XFLD[:,:]
        idx_land = nmp.where( mask < 0.5)
        xtmp[idx_land] = nmp.nan
        id_fld[:,:] = xtmp[:,:]
        del xtmp
    else:
        id_fld[:,:] = XFLD[:,:]
    f_o.about = cabout_nc
    f_o.close()
    return


def SaveTimeSeries( ivt, xd, vvar, ncfile, time_units='unknown', vunits=[], vlnm=[], missing_val=-9999. ):
    '''
    #  * ivt: time vector of length Nt, unit: UNIX Epoch time            [integer]
    #         => aka "seconds since 1970-01-01 00:00:00"
    #  *  xd: 2D numpy array that contains Nf time series of length Nt   [real]
    #          => hence of shape (Nf,Nt)
    #  * vvar: vector of length Nf of the Nf variable names                         [string]
    #  * vunits, vlnm: vectors of length Nf of the Nf variable units and long names [string]
    #  * missing_val: value for missing values...                        [real]
    '''
    (Nf,Nt) = xd.shape
    if len(ivt) != Nt: MsgExit('SaveTimeSeries() => disagreement in the number of records between "ivt" and "xd"')
    if len(vvar)!= Nf: MsgExit('SaveTimeSeries() => disagreement in the number of fields between "vvar" and "xd"')
    l_f_units = (nmp.shape(vunits)==(Nf,)) ; l_f_lnm = (nmp.shape(vlnm)==(Nf,))
    #
    print('\n *** About to write file "'+ncfile+'"...')
    f_o = Dataset(ncfile, 'w', format='NETCDF4')
    f_o.createDimension('time', None)
    id_t = f_o.createVariable('time','f8',('time',))
    id_t.calendar = 'gregorian' ; id_t.units = time_units
    #
    id_d = []
    for jf in range(Nf):
        id_d.append( f_o.createVariable(vvar[jf],'f4',('time',), fill_value=missing_val, zlib=True, complevel=5) )
        if l_f_units: id_d[jf].units   = vunits[jf]
        if l_f_lnm:   id_d[jf].long_name = vlnm[jf]
    #
    print('   ==> writing "time"')
    id_t[:] = ivt.astype(nmp.float64)
    for jf in range(Nf):
        print('   ==> writing "'+vvar[jf]+'"')
        id_d[jf][:] = xd[jf,:]
    f_o.about = cabout_nc
    f_o.close()
    print(' *** "'+ncfile+'" successfully written!\n')
    return 0
