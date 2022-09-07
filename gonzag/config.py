#!/usr/bin/env python3
# -*- Mode: Python; coding: utf-8; indent-tabs-mode: nil; tab-width: 4 -*-
#
#   ///// https://github.com/brodeau/gonzag \\\\\
#
#       L. Brodeau, 2021
#
############################################################################

ldebug = False

ivrb   = 1 ; # level of verbosity... [0-2]

IsZarr = False ; # A VIRER quand test sur extensions fichier prete!

nb_talk = 10 ; # how many times do you want to see a progression message in long loops ?

l_plot_meshes = True ; # if ldebug: then will generate a plot the sources meshes at the `nb_talk` frequency

deg2km = 111.11 ; # Converts degrees to km

R_v  = 6371.0   ; # Volumetric mean radius (km)
R_eq = 6378.137 ; # radius of earth at equator
R_pl = 6356.752 ; # radius of earth at the poles

rmissval = -9999. ; # Flag missing values in NetCDF files...

rfactor = 0.75 ; # this is multiplied to the model grid resolution to obtain the "found distance criterion"
#                #  => the larger the easier to find a nearest point (but maybe not the absolute nearest then...)

search_box_w_km = 500. ; # width (in km) of the small zoom-box on the source (model) domain
#                        # in which NearestPoint() will initially attempt to locate  the nearest point,
#                        # before falling back on the whole source (model) domain if unsuccessful.
#                        # => see "SearchBoxSize()" in "utils.py" for more !


l_save_track_on_model_grid = True



