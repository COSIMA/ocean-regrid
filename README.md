# ocean-ic

Create MOM and NEMO 'cold start' (temperature and salinity) initial conditions from global reanalysis.

# Use

Download a reanalysis dataset, many can be found here:

https://reanalyses.org/ocean/overview-current-reanalyses

The horizontal and vertical model grid definitions as well as the land-sea mask are also needed.

Example command creating a MOM initial condition from GODAS reanalysis.
```
$ ./makeic.py --temp_var pottmp --salt_var salt --ocean_mask ocean_mask.nc ocean_hgrid.nc ocean_vgrid.nc pottmp.2016.nc salt.2016.nc
```

Creating NEMO initial condition from GODAS:
```
$ ./makeic.py --temp_var pottmp --salt_var salt --model NEMO coordinates.nc data_1m_potential_temperature_nomask.nc pottmp.2016.nc salt.2016.nc
```

# How it works

1. The reanalysis/obs dataset is regridded in the vertical to have the same depth and levels as the model grid. Linear interpolation is used for this. If the model is deeper than the obs then the deepest value is extended.

2. Since the obs dataset is often limited latitudinally it is extended to cover the whole globe. This is done based on nearest neighbours.

3. The obs dataset is then regridded onto the model grid using bilinear interpolation.

4. The model land sea mask is applied and initial condition written out.

# Limitations

To use the created IC the ocean model needs to be able to start from rest.

