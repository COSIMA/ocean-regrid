# ocean-ic

Create MOM and NEMO 'cold start' (temperature and salinity) initial conditions from either GODAS or ORAS4 reanalysis.

# Use

Download a reanalysis dataset, many can be found here:

https://reanalyses.org/ocean/overview-current-reanalyses

The horizontal and vertical model grid definitions as well as the land-sea mask are also needed, in the case or ORAS4 this is a separate file, for GODAS it is contained within the data file.

Example command creating a MOM initial condition from GODAS reanalysis.
```
$ ./makeic.py --obs_name GODAS --model_mask ocean_mask.nc ocean_hgrid.nc ocean_vgrid.nc pottmp.2016.nc salt.2016.nc
```

Creating NEMO initial condition from GODAS:
```
$ ./makeic.py --model_name NEMO --obs_name GODAS coordinates.nc data_1m_potential_temperature_nomask.nc pottmp.2016.nc salt.2016.nc
```

Creating MOM initial conditions from ORAS4:
```
$ ./makeic.py --model_name MOM --obs_grid coords_T.nc --model_mask ocean_mask.nc ocean_hgrid.nc ocean_vgrid.nc thetao_oras4_1m_2014_grid_T.nc so_oras4_1m_2014_grid_T.nc --output mom_oras4_ic.nc
```

# How it works

1. The reanalysis/obs dataset is regridded in the vertical to have the same depth and levels as the model grid. Linear interpolation is used for this. If the model is deeper than the obs then the deepest value is extended.

2. In the case of GODAS since the obs dataset is limited latitudinally it is extended to cover the whole globe. This is done based on nearest neighbours.

3. The obs dataset is then regridded onto the model grid using bilinear interpolation.

4. The model land sea mask is applied and initial condition written out.

# Limitations

* Because the IC only includes salt and temperature the ocean model needs to be able to start from rest.
* When using GODAS reanalysis the values at high latitudes are unphysical due to limited observations.
* WARNING: presently the final interpolation assumes that the obs/reanalysis grid is a regular lat-lon grid, this is not the case for ORAS4.

# Example output

![MOM IC based on GODAS reanalysis](https://raw.github.com/nicjhan/ocean-ic/master/examples/MOM_IC_GODAS.png)

Notice the unusual lines at high latitude. This occurs because the GODAS temperature dataset is limited to 65 deg N, beyond that nearest neighbour values are being used.
