# ocean-ic

Create MOM and NEMO 'cold start' (temperature and salinity) initial conditions from global reanalysis. 

# Use

Download a reanalysis dataset, e.g.:

https://www.nodc.noaa.gov/OC5/woa13/woa13data.html
http://www.esrl.noaa.gov/psd/data/gridded/data.nodc.woa98.html
http://www.esrl.noaa.gov/psd/data/gridded/data.godas.html
https://reanalyses.org/ocean/overview-current-reanalyses

The horizontal and vertical model grid definitions as well as the land-sea mask are also needed.

Example command using GODAS
```
$ ./makeic.py --temp_var pottmp ocean_hgrid.nc ocean_vgrid.nc ocean_mask.nc pottmp.2016.nc pottmp.2016.nc
```

# How it works

1. The reanalysis/obs dataset is regridded in the vertical to have the same depth and levels as the model grid. Linear interpolation is used for this. If the model is deeper than the obs then the deepest value is extended.

2. Since the obs dataset is often limited it is extended to cover the whole globe. This is done based on nearest neighbours.

3. The obs dataset is then regridded onto the model grid using bilinear interpolation. 

4. The model land sea mask is applied and restart written out.

