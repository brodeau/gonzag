# gonzag

A Python3 package for the interpolation of 2D (O)GCM gridded dataset onto 1D satellite tracks.


## What I do

I interpolate, in space and time, the SSH (or any other 2D field) from a
non-regular structured gridded domain of an OGCM onto a satellite track.  The
satellite track is provided as time-series of the form time(t), longitude(t),
latitude(t), as usually found in along-track satellite data product.

Both SSH `2D+time` data and along track data are provided as netCDF files and must
share a common time-period.

2D space interpolation is performed through the bilinear method.

As an output gonzag ...
<br>

## Why

Ocean model horizontal grids are notorious for being twisted/distorted, hence, in gonzag,
special effort is put on thorough research of nearest point and source mesh
grid prior to bilinear interpolations.
<br>

## Dependencies

The following Python3 modules/packages are needed:
* `argparse`
* `numpy`
* `netCDF4`
* `time`
* `datetime`
* `calendar`
* `shapely`
* `matplotlib`
<br>

## How

![**plot**](https://github.com/brodeau/gonzag/blob/main/doc/figs/mesh_jt10500.png) <br>
*Figure 1: Example of a plot produced in debug-mode: `gonzag` let you know where each
nearest- and surrounding- points, as well as bilinear weights associated to each
of the 4 surrounding points, were taken to perform the bilinear interpolation.*

<br>


## Getting started

First download the archive containing the input files for the test suite:
https://drive.google.com/u/1/uc?id=1M5SBsMbUV29-3rGuh16X112hdI5NJ9qq&export=download

Best way is to learn from the examples! So let's perform some of the small following tests.

#### Light example,  global ORCA1 SSH interpolated to SARAL-AltiKa track

    ./alongtrack_sat_vs_nemo.py -s dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4 -n adt_unfiltered \
                                -m ssh_ORCA1_20170101_20171231_grid_T.nc4 -v ssh \
                                -l 0 \
                                -p 2

* `-s dt_global_alg_sla_vxxc_20170402_SARAL-Altika.nc4`: the file containing the 1D (time,lat,lon) satellite track
* `-n adt_unfiltered`: name of variable of interest in satellite track file (won't be used for any calculations, will just be saved in the output file together with model-interpolated tracks
* `-m ssh_ORCA1_20170101_20171231_grid_T.nc4` the file containing the 2D+t model variable to interpolate, with 2D latitude and longitude arrays
* `-v ssh` : name of variable of interest in model file
* `-l 0` : we get the model land-sea mask from the NetCDF `_FillValue` argument of the field `ssh`
* `-p 2` : ORCA1 is a global NEMO ORCA-type gridded domain, so there is an East-West periodicity of 2 overlapping points !

Check out the `xnp_msk.nc` file generated to see nearest-point satellite track on the 2D model grid.

In the output file `results.nc`, you will find time-series of model `ssh` interpolated on the satellite track (with bilinear and nearest point interpolation), as well as the original satellite field `adt_unfiltered` for the relevant time slice and region.

![**plot**](https://github.com/brodeau/gonzag/blob/main/doc/figs/track_ex_ORCA1.svg) <br>
*Figure 2: nearest-points of the satellite track located on the ORCA1 gridded domain, as computed in this example...*

<br>

#### Heavier example, SSH in eNATL60 zoom over the Faroe Islands interpolated to SARAL-AltiKa track

	./alongtrack_sat_vs_nemo.py -s dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc4 -n adt_unfiltered \
	                            -m sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc -v sossheig \
	                            -l dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc -k tmask \
	                            -p -1

* `-s dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc4`: the file containing the 1D (time,lat,lon) satellite track
* `-n adt_unfiltered`: name of variable of interest in satellite track file (won't be used for any calculations, will just be saved in the output file together with model-interpolated tracks
* `-m sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc` the file containing the 2D+t model variable to interpolate, with 2D latitude and longitude arrays
* `-v sossheig` name of variable of interest in model file
* `-l dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc` : we get the model land-sea mask in this file
* `-k tmask`: name of land-sea mask  in file `dt_global_alg_sla_vxxc_JFM_2017_SARAL-Altika.nc` variable is `tmask`
* `-p -1` : it's a rectangular extraction, so a regional region, so NO East-West periodicity!

Check out the `xnp_msk.nc` file generated to see nearest-point satellite track on the 2D model grid.

In the output file `results.nc`, you will find time-series of model `sossheig` interpolated on the satellite track (with bilinear and nearest point interpolation), as well as the original satellite field `adt_unfiltered` for the relevant time slice and region.


![**plot**](https://github.com/brodeau/gonzag/blob/main/doc/figs/track_ex_eNATL60-Faroe.svg) <br>
*Figure 3: nearest-points of the satellite track located on the eNATL60 Faroe zoom gridded domain, as computed in this example...*

<br>



## Production

Set `ldebug = False` in `gonzag/config.py` !

