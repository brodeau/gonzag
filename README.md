# gonzag

A Python3 package for the interpolation of 2D (O)GCM gridded (aka rasterized) data onto 1D satellite tracks.


## What I do

I interpolate, both in space and time, the *sea surface height* (or any other 2D field) provided onto the
horizontal grid (regular or irregular)  onto a 1-dimensional satellite track. The
satellite track is provided as a time-series of the form `time(t), longitude(t),
latitude(t)`, as usually found in along-track satellite data products.

The `2D+time` data of the model and the along-track data of the satellite must be provided in two separate netCDF files and must share a common period of time.

2D space interpolation is performed through the bilinear method.

<br>

## Why

Ocean model horizontal grids are notorious for being twisted/distorted, hence, in gonzag,
special effort is put on thorough research of nearest point and surrounding source mesh
grid points prior to bilinear interpolations.
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

<!-- ![**plot**](https://github.com/brodeau/gonzag/blob/main/doc/figs/mesh_jt10500.png)  -->

<p align="center">
  <img width="540" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/mesh_jt10500.png">
</p>

_Figure 1: Example of a plot produced in debug-mode: gonzag lets you know where each nearest- and surrounding- points, as well as the bilinear weights associated to each of the 4 surrounding points, were taken to perform the bilinear interpolation._

<br>

## Getting started

First download the archive containing the input files for the test suite:

https://drive.google.com/u/1/uc?id=1XZsEp41RckN9ulBDYZ4VFlh8ycvZLbfL&export=download


### You're more into `jupyter-notebooks`

Try the notebooks into the `notebook/` directory.
Currently the [Med-Sea example](https://github.com/brodeau/gonzag/blob/main/notebook/alongtrack_sat_vs_nemo.ipynb).


### You're more into writing scripts

`alongtrack_sat_vs_nemo.py` is a script that gets the job done!

```$ ./alongtrack_sat_vs_nemo.py -h```

    usage: alongtrack_sat_vs_nemo.py [-h] -s FSAT -n VSAT -m FMOD -v VMOD [-l FMSK] [-k VMSK] [-D] 

    I interpolate, in space and time, the SSH (or any other 2D field) from a non-regular structured gridded domain of an 
    OGCM onto a satellite track. 

    optional arguments: 
      -h, --help       show this help message and exit 
      -l FMSK, --fmsk FMSK  land-sea mask on model grid: "-l <ncfile>" or "-l 0" to use "_FillValue" attribute of field of interest 
      -k VMSK, --vmsk VMSK  name of land-sea mask in <ncfile> if "-l <ncfile>" used 
      -D, --distgrid     take into account possible strong local distortion of the model grid 
    
    required arguments: 
      -s FSAT, --fsat FSAT  satellite alongtrack data NetCDF file to read from 
      -n VSAT, --vsat VSAT  name of field of interest in satellite file (default="ssh") 
      -m FMOD, --fmod FMOD  model NetCDF file to read from 
      -v VMOD, --vmod VMOD  name of field of interest in model file (default="ssh")



Best way is to learn from the examples! So let's perform some of the small following tests.



#### Example 1 [lightweight]: monthly global ORCA1 SSH interpolated to SARAL-AltiKa track

Here, monthly-averaged SSH from a global NEMO-ORCA1 simulation is interpolated on a 1-day-long SARAL-AltiKa track (1st of August 2017).

    ./alongtrack_sat_vs_nemo.py -s SARAL_20170801-20170801.nc -n sla_unfiltered \
                                -m ssh_ORCA1_20170101_20171231_grid_T.nc -v ssh -l 0 -D

* `-s SARAL_20170801-20170801.nc`: the file containing the 1D (time,lat,lon) satellite track
* `-n sla_unfiltered`: name of variable of interest in satellite track file (won't be used for any calculations, will just be saved in the output file together with model-interpolated tracks
* `-m ssh_ORCA1_20170101_20171231_grid_T.nc` the file containing the 2D+t model variable to interpolate, with 2D latitude and longitude arrays
* `-v ssh` : name of variable of interest in model file
* `-l 0` : we get the model land-sea mask from the NetCDF `_FillValue` argument of the field `ssh`
* `-D` : as a global ORCA-type of grid, ORCA1 grid gets pretty distorted in the North, this option forces the computation of the grid distortion (local rotation), so that extra thoroughness is brought while building bilinear mapping in highly distorted regions.
<!-- * `-p 2` : ORCA1 is a global NEMO ORCA-type gridded domain, so there is an East-West periodicity of 2 overlapping points ! -->

Check out the `xnp_msk.nc` file generated to see nearest-point satellite track on the 2D model grid.

In the output file `results.nc`, you will find time-series of model `ssh` interpolated on the satellite track (with bilinear and nearest point interpolation), as well as the original satellite field `sla_unfiltered` for the relevant time slice and region.

<!-- ![**plot**](https://github.com/brodeau/gonzag/blob/main/doc/figs/track_ex_ORCA1.svg) <br> -->
<p align="center">
  <img width="400" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/track_ex_ORCA1.svg">
</p>

_Figure 2: Nearest-points of the satellite track located on the ORCA1 gridded domain, as computed in this example..._

<br>

<p align="center">
  <img width="400" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/ORCA1_MEAN_NEMO--Altimetry_pow-spectrum.svg">
</p>

_Figure 3: Spectral-analysis comparison of the virtual/model-based and satellite tracks. The tragedy of model data at coarse resolution (1 degree + monthly average): no eddies, no energy!_

<br>


#### Example 2: regional hourly SSH of GLORYS12 interpolated to SARAL-AltiKa track

Here, hourly SSH in a South-Western Pacific zoom of the GLORYS12 global reanalysis of Mercator Ocean at 12<sup>th</sup> of a degree (based on NEMO-ORCA12, available on CMEMS data-center) is interpolated on the SARAL-AltiKa track in March 2019.

    ./alongtrack_sat_vs_nemo.py -s SARAL_20190301-20190331.nc -n sla_unfiltered \
                                -m zos_GLORYS12V1_SouthWestPac_20190301-20190331_hourly.nc -v zos -l 0

Note: the `-D` option is not used here, because CMEMS provides NEMO-ORCA-based products on a regular `lat-lon` grid; besides, on the original ORCA12 grid, it would not have been necessary either as ORCA grids are not distorted in southern mid-latitudes.

<p align="center">
  <img width="400" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/track_ex_SouthWestPac_G12.svg">
</p>

_Figure 4: Nearest-points of the satellite track located on the GLORYS12 gridded domain, as computed in this example..._

<br>

<p align="center">
  <img width="400" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/SouthWestPac_G12_MEAN_NEMO--Altimetry_pow-spectrum.svg">
</p>

_Figure 5: Spectral-analysis comparison of the virtual/model-based and satellite tracks. You can see that despite the relative high resolution of the model (12<sup>th</sup> of a degree + hourly) the absence of tidal motion, and hence internal tides, in GLORYS12, yields a critical lack of energy at smaller scales!_

<br>



#### Example 3: regional hourly SSH of eNATL60 interpolated to Sentinel-3A track

Here, hourly SSH in a West Med zoom of the [eNATL60 simulation](https://github.com/ocean-next/eNATL60) with tidal motion at 60<sup>th</sup> of a degree is interpolated on the Sentinel-3A track in February 2017.

    ./alongtrack_sat_vs_nemo.py -s SENTINEL3A_20170130-20170303.nc -n sla_unfiltered \
                                -m sossheig_box_WestMed_eNATL60-BLBT02_20170201-20170228.nc -v sossheig -l 0

Note: the `-D` option is not used here as well, because the eNATL60 is only weakly distorted over the Med Sea.

<p align="center">
  <img width="400" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/track_ex_eNATL60-WestMed.svg">
</p>

_Figure 6: Nearest-points of the satellite track located on the eNATL60 zoom, as computed in this example..._

<br>

<p align="center">
  <img width="400" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/WestMed_MEAN_NEMO--Altimetry_pow-spectrum.svg">
</p>

_Figure 7: Spectral-analysis comparison of the virtual/model-based and satellite tracks._

<br>


#### Example 4 [heavy duty]: regional hourly SSH of eNATL60 interpolated to SARAL track

Here, hourly SSH in a [eNATL60](https://github.com/ocean-next/eNATL60) zoom over the Faeroe Islands  at 60<sup>th</sup> of a degree is interpolated on the SARAL track in JFM 2017.

    ./alongtrack_sat_vs_nemo.py -s SARAL_20170101-20170331.nc -n sla_unfiltered \
                                -m sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc -v sossheig \
                                -l sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc -k tmask -D

* `-l sossheig_box_Faroe_eNATL60-BLBT02_20170101-20170331.nc` : land-sea mask field is in this file
* `-k tmask` : name of land-sea mask field is "tmask"
* `-D` : high-latitude in an ORCA-based NEMO type of grid means substantial grid distortion


<p align="center">
  <img width="400" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/track_ex_eNATL60-Faroe.svg">
</p>

_Figure 8: Nearest-points of the satellite track located on the eNATL60 zoom, as computed in this example..._

<br>

<p align="center">
  <img width="400" src="https://github.com/brodeau/gonzag/blob/main/doc/figs/Faroe_MEAN_NEMO--Altimetry_pow-spectrum.svg">
</p>

_Figure 9: Spectral-analysis comparison of the virtual/model-based and satellite tracks. Model with super high resolution and tiny box: large eddies are not represented, energy at smaller scales in excellent agreement with satellite observations!_

<br>


## Production

Set `ldebug = False` in `gonzag/config.py` !

## FAQ

Spectral analysis and associated figures are done by means of script `plot_spectra_SSH_sat_track.py` of package [`climporn`](https://github.com/brodeau/climporn), which directly takes the `results.nc` generated by `alongtrack_sat_vs_nemo.py` as an input...

