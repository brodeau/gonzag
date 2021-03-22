#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-

############################################################################
#
#       L. Brodeau, 2021
#
############################################################################

from sys import exit
import numpy as nmp

import gonzag as gz
from gonzag.config import ldebug


## Global config:
xlon = nmp.array( [ [ 0.5, 180., 359.5 ] , [ 0.5, 180., 359.5 ] ] )
print(' ===>> ', gz.IsGlobalLongitudeWise2( xlon, resd=1.),'\n' )


## Regional domain with Greenwhich meridian inside:
xlon = nmp.array( [ [ 340., 358., 5. ] , [ 341., 359., 6. ] ] )
print(' ===>> ', gz.IsGlobalLongitudeWise2( xlon, resd=1.) )
print('   ( must return: False, False, -20.0, 6.0 )\n')
#exit(0)

## Regional domain that has the -180:+180 tansitiion inside:
xlon = nmp.array( [ [ 170., 178., 201. ] , [ 169., 181., 200. ] ] )
print(' ===>> ', gz.IsGlobalLongitudeWise2( xlon, resd=1.) )
print('   ( must return: False, True, 169, 201 )\n')


## Normal gentle regional location near 360, no trick:
xlon = nmp.array( [ [ 299., 324., 350. ] , [ 300., 325., 352. ] ] )
print(' ===>> ', gz.IsGlobalLongitudeWise2( xlon, resd=1.) )
print('   ( must return: False, True, 299., 352. )\n')


## Normal gentle regional location near 0, no trick:
xlon = nmp.array( [ [ 30., 40., 60. ] , [ 29., 39., 59. ] ] )
print(' ===>> ', gz.IsGlobalLongitudeWise2( xlon, resd=1.) )
print('   ( must return: False, True, 29., 60. )\n')
