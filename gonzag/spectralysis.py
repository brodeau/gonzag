#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#   ///// https://github.com/brodeau/gonzag \\\\\
#
#       L. Brodeau, 2021
#
############################################################################

import numpy as nmp
from .config import ivrb
from .utils  import MsgExit

l_tapper      = True ; # apply tappering !
l_detrend_lin = True ; # apply a linear detrending on data segment before computing spectrum...
#l_rm_i_noise  = False


# Because track files downloaded are usually full of inconsistencies / time and distance done
# between 2 consecutive points, we need to spot these inconsitencies with the 2 following criteria:
# ex: normally SARAL takes one measure every 1 s, => ds~=1 s, time during which it should move by
#     roughly dx ~= 7 km on the ground...
#     However, the netCDF files downloaded are full of points with dt=2s and the corresponding dx is
#     still dx ~= 7 km, this is clearly a bug !!!

#rcut_by_time = 1.2 # specify in seconds the time gap between two obs to detect and discontinuity and therefore
##                  # should be considered as a cut!

#rcut_by_dist = 7.8 # same as rcut_by_time, but in terms of distance (in km) between two consecutive points!
##                  # time criterion "rcut_by_time" would have been sufficient if SARAL data was not bugged!!!
##                  # => in SARAL data, spotted two consecutive points with the usual dt and a huge dL !!!
##                  #   => like if the satellite undergone an extremely short huge acceleration !!!
##                  #   => ex: 3rd of August 2016 around 07:53:43 !!!

#nlen_valid_seg_default = 120  # specify the minimum number of points a segment should contain to be considered and used!



#=============== en of configurable part ================================


def FindUnbrokenSegments( VTe, Vd, VM, rcut_time=1.2, rcut_dist=7.8 ):
    '''
    #   * Vd: full series of cumulated distances in km
    '''
    Nr = len(VTe)
    # Will extract the N valid data segments:
    nb_seg   = 0
    idx_strt = [] ; # index of first valid point of the segment
    idx_stop = [] ; # index of last  valid point of the segment
    Vmsk = VM.mask
    jr=0
    while jr < Nr:
        # Ignoring masked values and zeros...
        if (not Vmsk[jr]):
            vm = VM[jr]
            if (vm!=0.) and (vm<100.):
                nb_seg = nb_seg + 1
                if ivrb>1:
                    print('\n --- found seg #'+str(nb_seg)+' !')
                    print(' => starting at jt='+str(jr))
                idx_strt.append(jr)            
                while (not Vmsk[jr+1]) and (VM[jr+1]!=0.) and (VTe[jr+1]-VTe[jr] < rcut_time) and (Vd[jr+1]-Vd[jr] < rcut_dist) and (vm<100.) :
                    jr = jr+1
                    if jr==Nr-1: break
                idx_stop.append(jr)
                if ivrb>1: print(' => and stoping at jt='+str(jr))
        jr = jr+1
    if len(idx_strt) != nb_seg: MsgExit('[FindUnbrokenSegments()] => len(idx_strt) != nb_seg')
    return nmp.array(idx_strt, dtype=nmp.int64), nmp.array(idx_stop, dtype=nmp.int64)


def SegmentSelection(IS_start, IS_stop, np_valid_seg=100):
    '''
    #  
    #  1/ ID the continuous segments that have at least "np_valid_seg" points
    #  2/ Decide of a constant length for segments to stick with! => "Nsl"
    #  3/ ID the segments of length "Nsl"
    #  4/ Segments with more than "Nsl" points are trimmed at both extremities to only preserve
    #     the "Nsl" points at their center:

    #
    # Returns:
    #   * NbSeg: number of selected segments
    #   * Nsl:   trimmed length for all the  selected segments
    #   * IDEDSeg: array containing the start and end index position of each selected segment
    '''
    #
    NbP_seg = IS_stop[:] - IS_start[:] + 1 ; # number of points in each segment
    if ivrb>0: print('\n *** Longest segments has '+str(max(NbP_seg))+' points!\n')
    #
    # 1/ ID the continuous segments that have at least "np_valid_seg" points
    (idx_ok,) = nmp.where(NbP_seg >= np_valid_seg)
    NbSeg = len(idx_ok) ; # number of selected segments we are going to work with
    if NbSeg==0: MsgExit('could not find any valid segment with np_valid_seg = '+str(np_valid_seg))
    #
    # 2/ Decide of a constant length for segments to stick with => Nsl!
    rN = nmp.mean(NbP_seg[idx_ok])
    if ivrb>0: print(' *** [SegmentSelection()]: Mean segment-length for the '+str(NbSeg)+' segments with at least '+str(np_valid_seg)+' points:', round(rN,1))
    Nsl = int(rN/10.)*10 ; # We want a multiple of 10 just below this:
    if ivrb>1: print('  ==> Nsl = '+str(Nsl))
    #
    # 3/ ID the segments of length "Nsl"
    (idx_ok,) = nmp.where(NbP_seg >= Nsl)
    NbSeg = len(idx_ok)
    print(' *** [SegmentSelection()]: Will use '+str(NbSeg)+' segments with a fixed length of '+str(Nsl)+' points!')
    print('     ==> '+str(NbSeg)+' selected segments out of the '+str(len(IS_start))+' available (requested minimum length is '+str(np_valid_seg)+' points)\n')
    #
    # 4/ Segments with more than "Nsl" points to be trimmed at both extremities to only preserve
    #    the "Nsl" points at their center:
    #
    IDEDSeg = nmp.zeros( (NbSeg,2) , dtype=nmp.int64)
    for jp in range(NbSeg):
        js  = idx_ok[jp]        
        it1 = IS_start[js]
        it2 = IS_stop[js]
        nbp = it2-it1+1    
    
        if ivrb>0:
            print(' ###################################')
            print('  *** Seg #'+'%2.2i'%(jp+1)+' of '+cn_box+':')
            print('  ***   => originally '+str(nbp)+' points in this segment (from '+str(it1)+' to '+str(it2)+')')
    
        # nb of points in excess / Nsl:
        nxcs = nbp - Nsl
        jmp_strt = nxcs//2
        jmp_stop = nxcs//2 + nxcs%2
        it1 = it1+jmp_strt
        it2 = it2-jmp_stop
        if ivrb>0:
            print('  ***   => we only retain '+str(it2-it1+1)+' points from '+str(it1)+' to '+str(it2)+'!')
            print(' ###################################\n')
        #
        IDEDSeg[jp,:] = nmp.array([it1,it2])
    #
    del NbP_seg
    #
    return NbSeg, Nsl, IDEDSeg


def Process4FFT( IDseg, Vd, VS, VM ):
    '''
    # Process time-series on segments so that they are ready to undergo the Fast Fourrier Transform
    #
    # Input:
    #   * IDseg: array containing the start and end index position of each segment
    #   * Vd: full series of cumulated distances in km
    #   * VS: full series of SSH for satellite
    #   * VM: full series of SSH for model
    #
    # Output:
    #   * vs_s: processed satellite SSH for all segments [2D array of shape (<number of seg.>,<seg. length>)]
    #   * vs_m: processed model     SSH           "                      "
    #   * dx_sample: typical distance between to points in km
    #
    '''
    from scipy.signal import detrend, tukey
    #
    #
    NbS = IDseg.shape[0] ; # number of segments
    if ivrb>1: print(' *** [Process4FFT()] NbS = ', NbS)

    Nsl = IDseg[0,1]-IDseg[0,0]+1 ; # length of segments
    if ivrb>1: print(' *** [Process4FFT()] Nsl = ', Nsl)
    
    vs_s = nmp.zeros((NbS,Nsl))
    vs_m = nmp.zeros((NbS,Nsl))
    
    for js in range(NbS):
            
        it1 = IDseg[js,0]
        it2 = IDseg[js,1]

        if js==0:
            # Checking the typical distance (in km) between two measures:
            dmean = nmp.mean(Vd[it1+1:it2+1]-Vd[it1:it2])
            if ivrb>0: print(' *** [Process4FFT()]: Mean distance between two consecutive points is '+str(dmean)+' km\n')
            # Sample spacing in [km] (inverse of the sampling rate):
            dx_sample = round(dmean,3)
            if ivrb>0: print('     => will use a spatial sample spacing of '+str(dx_sample)+' km\n')
                    
        vs_s[js,:] = VS[it1:it2+1]
        vs_m[js,:] = VM[it1:it2+1]
        
        # Linear detrending
        if l_detrend_lin:
            if js==0: print(' *** [Process4FFT()]: applying linear detrending...')
            vs_s[js,:] = detrend(vs_s[js,:],type='linear')
            vs_m[js,:] = detrend(vs_m[js,:],type='linear')
        
        # Centering about 0:
        if js==0: print(' *** [Process4FFT()]: centering about 0...')
        vs_s[js,:] = vs_s[js,:] - nmp.mean(vs_s[js,:])
        vs_m[js,:] = vs_m[js,:] - nmp.mean(vs_m[js,:])
        
        # Tappering:
        if l_tapper:
            wdw =  tukey(Nsl,0.5)
            if js==0: print(' *** [Process4FFT()]: applying "Tukey" tappering...')
            vs_s[js,:] = vs_s[js,:]*wdw
            vs_m[js,:] = vs_m[js,:]*wdw
        if js==NbS-1: print('')
    return vs_s, vs_m, dx_sample


def ApplyFFT( IDseg, XS, XM, dx_sample ):
    '''
    # Process time-series on segments so that they are ready to undergo the Fast Fourrier Transform
    #
    # Input:
    #   * IDseg: array containing the start and end index position of each segment
    #   * XS: processed satellite SSH for all segments [2D array of shape (<number of seg.>,<seg. length>)]
    #   * XM: processed model     SSH           "                      "
    #   * dx_sample: typical distance between to points in km
    #
    # Output:
    #   * VK:   wave numbers associated with spectra [1D array of length <seg. length>]
    #   * PS_s: power spectrum of satellite SSH for all segments [2D array of shape (<number of seg.>,<seg. length>)]
    #   * PS_m: power spectrum of model     SSH           "                      "
    #
    '''
    (NbS,Nsl) = XS.shape
    VK   = nmp.zeros(     Nsl )
    PS_s = nmp.zeros((NbS,Nsl))
    PS_m = nmp.zeros((NbS,Nsl))

    print(' *** [ApplyFFT()]: Applying FFT with a dx_sample of ',dx_sample,' km')
    for js in range(NbS):
        if ivrb>0: print(' *** [ApplyFFT()]: applying FFT to segment #',js+1)
        
        it1 = IDseg[js,0]
        it2 = IDseg[js,1]

        # Wave numbers:
        if js==0:
            VK[:] = nmp.fft.fftfreq(Nsl, dx_sample)
            idx   = nmp.argsort(VK)
            VK[:] = VK[idx]
        
        # Power Spectrum:
        vYf_s = 2.*(dx_sample/float(Nsl)) * nmp.abs(nmp.fft.fft(XS[js,:]))**2
        vYf_m = 2.*(dx_sample/float(Nsl)) * nmp.abs(nmp.fft.fft(XM[js,:]))**2

        PS_s[js,:] = vYf_s[idx]
        PS_m[js,:] = vYf_m[idx]
        print('')
    return VK, PS_s, PS_m
