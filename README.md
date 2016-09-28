# ocean-regrid

Regrid global ocean data in 3d. Suppurts GODAS, ORAS4 reanalysis grids and MOM, NEMO model grids. Handles missing data and grids with mismatched domains.

# Build status

[![Build Status](https://travis-ci.org/nicjhan/ocean-regrid.svg?branch=master)](https://travis-ci.org/nicjhan/ocean-regrid)

# Install

## Python dependencies

1. Download and install [Anaconda](https://www.continuum.io/downloads) for your platform.
2. Install the [git](https://git-scm.com/) revision control system if you don't already have it.
3. Download ocean-regrid:
```{bash}
$ git clone --recursive https://github.com/nicjhan/ocean-regrid.git
$ cd ocean-regrid
```
4. Setup the Anaconda environment. This will download all the necessary Python packages.
```{bash}
$ conda env create -f regrid.yml
$ source activate regrid
```

## ESMF dependencies

Install ESMF_RegridWeightGen. ESMF releases can be found [here](http://www.earthsystemmodeling.org/download/data/releases.shtml).

There is a bash script contrib/build_esmf.sh which the testing system uses to build ESMF. This may be useful in addition to the ESMF installation docs.

# Example Use

Regrid ORAS4 reanalysis to MOM 0.25 degree tripolar grid:
```
$ cd test
$ wget http://s3-ap-southeast-2.amazonaws.com/dp-drop/ocean-regrid/test/test_data.tar.gz
$ tar zxvf test_data.tar.gz
$ cd test_data/input
$ ../../../regrid.py ORAS4 coords_T.nc coords_T.nc thetao_oras4_1m_2014_grid_T.nc thetao \
    MOM ocean_hgrid.nc ocean_vgrid.nc mom_oras4_temp.nc temp --dest_mask ocean_mask.nc
$ ncview mom_oras4_temp.nc
```

OR:

```
$ python -m pytest
$ ncview test/test_data/output/mom_oras4_temp.nc
```

# How it works

1. The source dataset is regridded in the vertical to have the same depth and levels as the destination grid. Linear interpolation is used for this. If the destination is deeper than the source then the deepest value is extended.

2. If the source dataset is limited latitudinally it is extended to cover the whole globe. This is done based on nearest neighbours.

3. The source dataset is then regridded using weights calculated with ESMF_RegridWeightGen. Various regridding schemes are supported includeing distance weighted nearest neighbour, bilinear and conservative.

4. The destination land sea mask is applied.

