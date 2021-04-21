#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#   ///// https://github.com/brodeau/gonzag \\\\\
#
#       L. Brodeau, 2021
#
############################################################################

import numpy as nmp
import xarray as xr
from netCDF4 import num2date, default_fillvals
from calendar import timegm
from calendar import timegm
from datetime import datetime as dtm
from .config import ldebug, ivrb, rmissval
from .utils  import MsgExit

cabout_nc = 'Created with Gonzag package => https://github.com/brodeau/gonzag'


def ToEpochTime( X, units, calendar ):
    '''
    # INPUT:
    #   * X: time vector or scalar, provided as something like: "seconds/hours/days since <DATE0>"
    #   * units, calendar: units and calendar type [char]
    #
    # OUTPUT:
    #   * vet: time vector converted to UNIX Epoch time,
    #          aka "seconds since 1970-01-01 00:00:00"
    #          => as FLOAT !! not INTEGER !!!
    '''
    cfrmt = '%Y-%m-%d %H:%M:%S'
    #
    lvect = not (nmp.shape(X)==()) ; # True => X is not a scalar
    #
    if lvect:
        nt = len(X)
        t0 = X[0]
    else:
        t0 = X
    if ivrb>0: print(' *** [ToEpochTime()]: original t0 as "'+units+'" => ', t0)
    t0d = num2date( t0, units, calendar )
    if ivrb>0: print(' *** [ToEpochTime()]: intitial date in datetime format => ', t0d)
    #
    rdec = t0d.microsecond*1.E-6 ; # to convert micro seconds part to decimal part of time in seconds...
    t0E = float( timegm( dtm.strptime( t0d.strftime(cfrmt) , cfrmt ).timetuple() ) + rdec ) ; # t0 as "float" Epoch UNIX time!
    #
    if lvect:
        if   units[0:13] == 'seconds since':
            r2s = 1.
        elif units[0:11] == 'hours since':
            r2s = 3600.
        elif units[0:10] == 'days since':
            r2s = 86400.
        else:
            MsgExit('[ToEpochTime()] => unsupported time unit: '+units)
        #
        vet     = nmp.zeros(nt)
        vet[0]  = t0E
        vet[1:] = t0E + (X[1:] - X[0])*r2s
        #
    else:
        vet = t0E
    #
    return vet



def GetTimeInfo( ncfile ):
    '''
    # Inspect time dimension
    # Get number of time-records + first and last date
    # Return them + dates as UNIX Epoch time, aka "seconds since 1970-01-01 00:00:00" (float)
    '''
    if ldebug: print(' *** [GetTimeInfo()] Getting calendar/time info in '+ncfile+' ...')
    id_f = xr.open_zarr(ncfile,decode_cf=False)
    for cd in [ 'time', 'time_counter', 'TIME', 'record', 't', 'none' ]:
        if cd in id_f.coords: break
    if cd == 'none': MsgExit('found no time-record dimension in file '+ncfile)
    if ldebug: print('   => time/record dimension is "'+cd+'"')
    nt = id_f[cd].size
    clndr = id_f[cd]

    dt1 = num2date( clndr[0], clndr.units, clndr.calendar ) ; dt2 = num2date( clndr[nt-1], clndr.units, clndr.calendar )
    rt1 = ToEpochTime( clndr[0],    clndr.attrs['units'], clndr.attrs['calendar'] )
    rt2 = ToEpochTime( clndr[nt-1], clndr.attrs['units'], clndr.attrs['calendar'] )
    id_f.close()
    #
    if ldebug: print('   => first and last time records: ',dt1,'--',dt2,' (UNIX Epoch: ', rt1,'--',rt2,')\n')
    #
    return nt, (rt1,rt2)




def GetTimeEpochVector( ncfile, kt1=0, kt2=0, isubsamp=1, lquiet=False ):
    '''
    # Get the time vector in the netCDF file and returns it as UNIX Epoch time,
    # aka "seconds since 1970-01-01 00:00:00" (float!)
    #
    # INPUT:
    #  * kt1, kt2 : read from index kt1 to kt2 (if these 2 are != 0)
    #  * isubsamp: subsampling !!!
    #  * lquiet: shut the f* up!
    #
    # OUTPUT:
    #  * rvte: time vector in (float) UNIX Epoch time
    '''
    ltalk = ( ldebug and not lquiet )
    cv_t_test = [ 'time', 'time_counter', 'TIME', 'record', 't', 'none' ]
    id_f = xr.open_zarr(ncfile,decode_cf=False)
    for cv in cv_t_test:
        if cv in id_f.coords:
            clndr = id_f[cv]
            cunt  = clndr.attrs['units']
            ccal  = clndr.attrs['calendar']
            if ivrb>0 and ltalk: print(' *** [GetTimeEpochVector()] reading "'+cv+'" in '+ncfile+' and converting it to Epoch time...')
            if kt1>0 and kt2>0:
                if kt1>=kt2: MsgExit('mind the indices when calling GetTimeEpochVector()')
                vdate = clndr[kt1:kt2+1:isubsamp]
                cc = 'read...'
            else:
                vdate = clndr[::isubsamp]
                cc = 'in TOTAL!'
            break
    id_f.close()
    if cv == 'none': MsgExit('found no time-record variable in file '+ncfile+' (possible fix: "cv_t_test" in "GetTimeEpochVector()")')
    #
    rvte = ToEpochTime( vdate, cunt, ccal ) ; # convert to Unix Epoch time
    #
    if ivrb>0 and ltalk: print('   => '+str(len(rvte))+' records '+cc+'\n')
    return rvte



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
    id_f = xr.open_zarr(ncfile)
    for ncvar in cv_coor_test[ii,:]:
        if ncvar in id_f.coords: break
    if ncvar == 'none': MsgExit('could not find '+what+' array into model file (possible fix: "cv_coor_test" in "GetModelCoor()")')
   #
    nb_dim = len(id_f[ncvar].dims)
    if   nb_dim==1: xwhat = id_f[ncvar][:]
    elif nb_dim==2: xwhat = id_f[ncvar][:,:]
    elif nb_dim==3: xwhat = id_f[ncvar][0,:,:]
    else: MsgExit('FIX ME! Model '+what+' has a weird number of dimensions')
    id_f.close()
    if ldebug: print(' *** [GetModelCoor()] Read model '+what+' (variable is "'+ncvar+'", with '+str(nb_dim)+' dimensions!',nmp.shape(xwhat),'\n')
    #
    return xwhat.values


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
    id_f = xr.open_zarr(ncfile)
    ndim = len(id_f[ncvar].dims)
    if l_fill_val:
        # Mask is constructed out of variable and its missing value
        if not ndim in [3,4]: MsgExit(ncvar+' is expected to have 3 or 4 dimensions')
        if ndim==3: xmsk = 1 - nmp.isnan(id_f[ncvar][0,:,:])
        if ndim==4: xmsk = 1 - nmp.isnan(id_f[ncvar][0,0,:,:])
    else:
        # Mask is read in mask file...
        if   ndim==2: xmsk = id_f[ncvar][:,:]
        elif ndim==3: xmsk = id_f[ncvar][0,:,:]
        elif ndim==4: xmsk = id_f[ncvar][0,0,:,:]
        else: MsgExit('FIX ME! Mask '+ncvar+' has a weird number of dimensions:'+str(ndims))
    #
    id_f.close()
    return xmsk.astype(int)


def GetModel2DVar( ncfile, ncvar, kt=0 ):
    '''
    #   Fetches the 2D field "ncvar" at time record kt into "ncfile"
    '''
    if ldebug: print(' *** [GetModel2DVar()] Reading model "'+ncvar+'" at record kt='+str(kt)+' in '+ncfile)
    id_f = xr.open_zarr(ncfile)
    nb_dim = len(id_f[ncvar].dims)
    if nb_dim==3:
        x2d = id_f[ncvar][kt,:,:]
    elif nb_dim==4:
        x2d = id_f[ncvar][kt,0,:,:] ; # taking surface field!
    else: MsgExit('FIX ME! Model "'+ncvar+'" has a weird number of dimensions: '+str(nb_dim))
    id_f.close()
    if ldebug: print('')
    return x2d




def GetSatCoor( ncfile, what,  kt1=0, kt2=0 ):
    '''
    # Get latitude (what=='latitude') OR longitude (what=='longitude') vector
    # in the netCDF file, (from index kt1 to kt2, if these 2 are != 0!)
    '''
    cv_coor_test = nmp.array([[ 'lat','latitude', 'LATITUDE',  'none' ],
                              [ 'lon','longitude','LONGITUDE', 'none' ]])
    if   what ==  'latitude': ii = 0
    elif what == 'longitude': ii = 1
    else: MsgExit('"what" argument of "GetSatCoor()" only supports "latitude" and "longitude"')
    #
    id_f = xr.open_zarr(ncfile)
    for ncvar in cv_coor_test[ii,:]:
        if ncvar in id_f.coords: break
    if ncvar == 'none': MsgExit('could not find '+what+' array into satellite file (possible fix: "cv_coor_test" in "GetSatCoor()")')
    #
    if ldebug: print(' *** [GetSatCoor()] reading "'+ncvar+'" in '+ncfile+' ...')
    nb_dim = len(id_f[ncvar].dims)
    if nb_dim==1:
        if kt1>0 and kt2>0:
            if kt1>=kt2: MsgExit('mind the indices when calling GetSatCoor()')
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
    return vwhat.values


def GetSatSSH( ncfile, ncvar,  kt1=0, kt2=0, ikeep=[] ):
    '''
    # Get vector time-series of 'ncvar' in file 'ncfile'!
    #  - from index kt1 to kt2, if these 2 are != 0!
    #  - if (ikeep != []) => only retains parts of the data for which indices are provide into ikeep
    #          (ikeep is an array obtained as the result of a "numpy.where()"
    '''
    if ldebug: print(' *** [GetSatSSH()] Reading satellite "'+ncvar+' in '+ncfile)
    id_f = xr.open_zarr(ncfile)
    if kt1>0 and kt2>0:
        if kt1>=kt2: MsgExit('mind the indices when calling GetSatSSH()')
        vssh = id_f[ncvar][kt1:kt2+1]
    else:
        vssh = id_f[ncvar][:]
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
    #Mask
    (nj,ni) = nmp.shape(XFLD)
    if nmp.shape(mask) != (0,):
        xtmp = nmp.zeros((nj,ni))
        xtmp[:,:] = XFLD[:,:]
        idx_land = nmp.where( mask < 0.5)
        xtmp[idx_land] = nmp.nan
    else:
        xtmp=XFLD

    #Turn fields into dataset
    foo=xr.DataArray(xtmp,dims=['y','x'])
    foo.name=name
    if unit      != '': foo.attrs["units"] = unit
    if long_name != '': foo.attrs["long_name"] = long_name

    #Save to netcdf
    ds=foo.to_dataset()
    ds.attrs=dict(about=cabout_nc)
    ds.to_netcdf(ncfile,'w', format='NETCDF4')

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

    if Nf == 1:
      foo=xr.DataArray(xd,dims=['time'],coords=[ivt.astype(nmp.float64)])
      foo.name=vvar
      foo.attrs["units"] = vunits
      foo.attrs["long_name"] = vlnm
      ds=foo.to_dataset(name = vvar)
      footime=xr.DataArray(ivt, dims=['time'],coords=[ivt.astype(nmp.float64)])
      footime.attrs["units"] =time_units
      footime.attrs["calendar"] = 'gregorian'
      ds['time']=footime
      ds.attrs=dict(about=cabout_nc)
      ds.to_netcdf(ncfile,'w', format='NETCDF4')
    else:
      foo0=xr.DataArray(xd[0],dims=['time'],coords=[ivt.astype(nmp.float64)])
      foo0.name=vvar[0]
      foo0.attrs["units"] = vunits[0]
      foo0.attrs["long_name"] = vlnm[0]
      ds=foo0.to_dataset(name = vvar[0])
      footime=xr.DataArray(ivt, dims=['time'],coords=[ivt.astype(nmp.float64)])
      footime.attrs["units"] =time_units
      footime.attrs["calendar"] = 'gregorian'
      ds['time']=footime
      for jf in range(1,Nf):
        foo=xr.DataArray(xd[jf],dims=['time'],coords=[ivt.astype(nmp.float64)])
        foo.name=vvar[jf]
        foo.attrs["units"] = vunits[jf]
        foo.attrs["long_name"] = vlnm[jf]
        ds[vvar[jf]]=foo
      ds.attrs=dict(about=cabout_nc)
      ds.to_netcdf(ncfile,'w', format='NETCDF4')


    print(' *** "'+ncfile+'" successfully written!\n')
    return 0
