# ocean-regrid

Regrid GODAS or ORAS4 reanalysis to MOM or NEMO grids. In general this means regridding from regular lat-lon (GODAS) or tripolar (ORAS4) to tripolar (MOM and NEMO) grids.

# Dependencies

This tool is written in Python and depends on many different Python packages. See section 'Install' below for instructions on how to download all of the Python dependencies. It also depends on
 [ESMF_RegridWeightGen](https://www.earthsystemcog.org/projects/regridweightgen/) program to perform regridding between non-rectilinear grids.

# Install

ESMF releases can be found here: http://www.earthsystemmodeling.org/download/data/releases.shtml

# Use

Download a reanalysis dataset, many can be found here:

https://reanalyses.org/ocean/overview-current-reanalyses

The horizontal and vertical model grid definitions as well as the land-sea mask are also needed, in the case or ORAS4 this is a separate file, for GODAS it is contained within the data file.

Example command regridding GODAS reanalysis to MOM:
```
$ ./regrid.py --obs_name GODAS --model_mask ocean_mask.nc ocean_hgrid.nc ocean_vgrid.nc pottmp.2016.nc salt.2016.nc
```

Creating NEMO initial condition from GODAS:
```
$ ./regrid.py --model_name NEMO --obs_name GODAS coordinates.nc data_1m_potential_temperature_nomask.nc pottmp.2016.nc salt.2016.nc
```

Creating MOM initial conditions from ORAS4:
```
$ ./regrid.py ORAS4 coords_T.nc coords_T.nc thetao_oras4_1m_2014_grid_T.nc thetao MOM ocean_hgrid.nc ocean_vgrid.nc ocean_out.nc temp --dest_mask ocean_mask.nc
```

Creating NEMO initial condisionf from ORAS4:
```
$ ./regrid.py --model_name MOM --obs_grid coords_T.nc --model_mask ocean_mask.nc ocean_hgrid.nc ocean_vgrid.nc thetao_oras4_1m_2014_grid_T.nc so_oras4_1m_2014_grid_T.nc --output mom_oras4_ic.nc
```

# Testing

# How it works

1. The reanalysis/obs dataset is regridded in the vertical to have the same depth and levels as the model grid. Linear interpolation is used for this. If the model is deeper than the obs then the deepest value is extended.

2. In the case of GODAS since the obs dataset is limited latitudinally it is extended to cover the whole globe. This is done based on nearest neighbours.

3. The obs dataset is then regridded onto the model grid using weights calculated with ESMF_RegridWeightGen. Various regridding schemes are supported includeing distance weighted nearest neighbour, bilinear and conservative.

4. The model land sea mask is applied and initial condition written out.

# Limitations

* When using GODAS reanalysis the values at high latitudes are unphysical due to limited observations.

# Example output

![MOM IC based on GODAS reanalysis](https://raw.github.com/nicjhan/ocean-ic/master/examples/MOM_IC_GODAS.png)

