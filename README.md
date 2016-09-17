# ocean-regrid

Regrid global ocean data in 3d. Suppurts GODAS, ORAS4 reanalysis grids and MOM, NEMO model grids. Handles missing data and grids with mismatched domains.

# Dependencies

This tool is written in Python and depends on many different Python packages. See section 'Install' below for instructions on how to download all of the Python dependencies. It also depends on
 [ESMF_RegridWeightGen](https://www.earthsystemcog.org/projects/regridweightgen/) program to perform regridding between non-rectilinear grids.

# Install

ESMF releases can be found here: http://www.earthsystemmodeling.org/download/data/releases.shtml

# Example Use

The horizontal and vertical model grid definitions as well as the land-sea mask are needed.

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

1. The source dataset is regridded in the vertical to have the same depth and levels as the destination grid. Linear interpolation is used for this. If the destination is deeper than the source then the deepest value is extended.

2. If the source dataset is limited latitudinally it is extended to cover the whole globe. This is done based on nearest neighbours.

3. The source dataset is then regridded using weights calculated with ESMF_RegridWeightGen. Various regridding schemes are supported includeing distance weighted nearest neighbour, bilinear and conservative.

4. The destination land sea mask is applied.

# Example output

![MOM IC based on GODAS reanalysis](https://raw.github.com/nicjhan/ocean-ic/master/examples/MOM_IC_GODAS.png)

