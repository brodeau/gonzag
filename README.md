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

