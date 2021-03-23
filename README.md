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

Download the test data:



Best way is to learn from the examples! So let's perform some of the small following tests.



#### Light example,  global ORCA1 ssh interpolated to SARAL-AltiKa track

    ./alongtrack_sat_vs_nemo.py -s dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4 -n adt_unfiltered \
                                -m ssh_ORCA1_20170101_20171231_grid_T.nc4 -v ssh \
                                -l 0 \
                                -p 2

* `-s dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc4`: the file containing the 1D (time,lat,lon) satellite track
* `-n adt_unfiltered`: name of variable of interest in satellite track file (won't be used for any calculations, will just be saved in the output file together with model-interpolated tracks
* `-m sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc` the file containing the 2D+t model variable to interpolate, with 2D latitude and longitude arrays
* `-v sossheig` name of variable of interest in model file
* `-l dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc` : we get the model land-sea mask in this file
* `-k tmask`: name of land-sea mask variable is `tmask`
* `-p 2` : ORCA1 is a global NEMO ORCA-type grid so there is an East-West perdiodicity of 2 overlaping points !




#### eNATL60 zoom over the Faroe Islands interpolated to SARAL-AltiKa track

	./alongtrack_sat_vs_nemo.py -s dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc4 -n adt_unfiltered \
							    -m sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc -v sossheig \
	                            -l dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc -k tmask \
	                            -p -1

* `-s dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc4`: the file containing the 1D (time,lat,lon) satellite track
* `-n adt_unfiltered`: name of variable of interest in satellite track file (won't be used for any calculations, will just be saved in the output file together with model-interpolated tracks
* `-m sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc` the file containing the 2D+t model variable to interpolate, with 2D latitude and longitude arrays
* `-v sossheig` name of variable of interest in model file
* `-l dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc` : we get the model land-sea mask in this file
* `-k tmask`: name of land-sea mask variable is `tmask`
* `-p -1` : it's a rectangular extraction, so a regional region, so no East-West perdiodicity



