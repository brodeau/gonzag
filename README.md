# gonzag

Interpolation of 2D GCM gridded data onto 1D satellite tracks


## What I do

I interpolate, in space and time, the SSH (or any other 2D field) from a
non-regular structured gridded domain of an OGCM onto a satellite track.  The
satellite track is provided as time-series of the form time(t), longitude(t),
latitude(s), as usually found in along-track satellite product.

Both SSH 2D+time data and along track data are provided as netCDF file and must
share a common time-period.

2D space interpolation is done via the bilinear method.

As an output gonzag ...


## Why

Ocean model horizontal grids are notorious for being twisted/distorded, hence, in gonzag,
special effort is put on thourough research of nearest point and source mesh
grid prior to bilinear interpolations.


## Getting started

Best way it to perform the small tests...


eNATL60, zoom over the Faroe Islands:


    ./alongtrack_sat_vs_nemo.py -s dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc4 -n adt_unfiltered \
                                -m sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc -v sossheig \
                                -l mesh_mask_eNATL60_Faroe.nc -p -1

* `-p -1` : it's a rectangular extraction, so a regional region, so no East-West perdiodicity
* `-l mesh_mask_eNATL60_Faroe.nc` : we get the model land-sea mask in this file



