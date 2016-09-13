#!/usr/bin/env python

from __future__ import print_function

import sys, os
import time
import argparse
import numpy as np
import numba 
import subprocess as sp
import netCDF4 as nc
from scipy import interp
from scipy import ndimage as nd

from mom_grid import MomGrid
from nemo_grid import NemoGrid
from regular_grid import RegularGrid
from tripolar_grid import TripolarGrid
from godas_grid import GodasGrid
from oras_grid import OrasGrid

from file_util import write_mom_ic, write_nemo_ic
from util import normalise_lons

"""
Create ocean model IC based on reanalysis data.
"""

GODAS_BERING_STRAIGHT = [416, 184]

def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()

def regrid_columns(data, src_grid, dest_grid):
    """
    Regrid vertical columns of data from src to dest.
    """

    assert len(data.shape) == 3
    assert data.shape[0] == len(src_grid.z)

    # This gets modified.
    data = np.ma.copy(data)

    # Create masked array of the correct shape.
    tmp = np.zeros((len(dest_grid.z), data.shape[1], data.shape[2]))
    new_data = np.ma.array(tmp, mask=np.ones_like(tmp))

    # Iterate through columns and regrid each.
    for lat in range(data.shape[1]):
        for lon in range(data.shape[2]):
            if np.all(src_grid.mask[:, lat, lon]):
                continue
            if np.any(src_grid.mask[:, lat, lon]):
                # Masked values expected to be at depth. Find these and fill
                # with nearest neighbour.
                for d in range(data.shape[0]-2, 0, -1):
                    if (not src_grid.mask[d, lat, lon]) and \
                            src_grid.mask[d+1, lat, lon]:
                        data[d:, lat, lon] = data[d, lat, lon]

            # 1d linear interpolation/extrapolation
            new_data[:, lat, lon] = interp(dest_grid.z, src_grid.z, data[:, lat, lon])

    return new_data


def extend_obs(var, obs_grid, global_grid, arctic_filler=None):
    """
    Use nearest neighbour to extend obs over the whole globe, including land.
    """

    # Create masked array of the correct shape.
    tmp = np.zeros((global_grid.num_levels, global_grid.num_lat_points,
                    global_grid.num_lon_points))
    new_data = np.ma.array(tmp, mask=global_grid.mask)

    # Drop obs data into new grid at correct location
    lat_min_idx = find_nearest_index(global_grid.y_t[:, 0], np.min(obs_grid.y_t[:]))
    if np.max(global_grid.y_t[:]) <= np.max(obs_grid.y_t[:]):
        new_data[:, lat_min_idx:, :] = var[:]
    else:
        lat_max_idx = find_nearest_index(global_grid.y_t[:, 0], np.max(obs_grid.y_t[:]))
        new_data[:, lat_min_idx:lat_max_idx+1, :] = var[:]

    # Fill in missing values on each level with nearest neighbour
    for l in range(var.shape[0]):
        ind = nd.distance_transform_edt(new_data[l, :, :].mask,
                                        return_distances=False,
                                        return_indices=True)
        tmp = new_data[l, :, :]
        tmp = tmp[tuple(ind)]
        new_data[l, :, :] = tmp[:, :]

    return new_data


def fill_arctic(obs_data, global_data, obs_grid, global_grid):
    """
    In the case of GODAS, there is no data past 65N and the nearest neighbour
    approach is not ideal because the low salinity of the Baltic gets
    propogated into the Arctic. So instead fill the Arctic with a
    'representative value' taken from the Bering Straight.
    """

    assert 'GODAS' in obs_grid.description

    GODAS_ARCTIC_REPRESENTATIVE_VALUE = GODAS_BERING_STRAIGHT
    filler = obs_data[:, GODAS_ARCTIC_REPRESENTATIVE_VALUE[0],
                         GODAS_ARCTIC_REPRESENTATIVE_VALUE[1]]

    arctic_idx = find_nearest_index(global_grid.y_t[:, 0],
                                    np.max(obs_grid.y_t[:, 0]))
    arctic_idx -= 1
    sh = global_data[:, arctic_idx:, :].shape

    filler = np.stack([filler[:]] * sh[1], axis=1)
    filler = np.stack([filler[:]] * sh[2], axis=2)
    global_data[:, arctic_idx:, :] = filler[:]

    return global_data

@numba.jit
def apply_weights(src, dest_shape, n_s, n_b, row, col, s):
    """
    Apply ESMF regirdding weights.
    """

    dest = np.ndarray(dest_shape).flatten()
    dest[:] = 0.0
    src = src.flatten()

    for i in xrange(1, n_s):
        dest[row[i]-1] = dest[row[i]-1] + s[i]*src[col[i]-1]

    return dest.reshape(dest_shape)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('model_hgrid', help='Model horizontal grid spec file.')
    parser.add_argument('model_vgrid', help='Model vertical grid spec file.')
    parser.add_argument('temp_obs_file', help='Temp from GODAS obs reanalysis.')
    parser.add_argument('salt_obs_file', help='Salt from GODAS obs reanalysis.')
    parser.add_argument('--model_mask', default=None, help='Model land-sea mask file.')
    parser.add_argument('--model_name', default='MOM',
                        help='Which model to create IC for, can be MOM or NEMO.')
    parser.add_argument('--obs_name', default='ORAS4',
                        help='Which obs reanalysis to use GODAS or ORAS4.')
    parser.add_argument('--obs_grid', default=None,
                        help='Observations grid definition file, must be provided for ORAS4')
    parser.add_argument('--output', default='ocean_ic.nc',
                        help='Name of the output file.')
    parser.add_argument('--time_index', default=0,
                        help='The time index of the data to use.')
    parser.add_argument('--regrid_tool', default='scipy',
                        help='The regridding method to use, scipy or esmf.')
    parser.add_argument('--regrid_weights', default=None,
                        help='The regridding weights file to use.')

    args = parser.parse_args()

    assert args.model_name == 'MOM' or args.model_name == 'NEMO'
    assert args.obs_name == 'GODAS' or args.obs_name == 'ORAS4'

    # FIXME: if using esmf check that we have ESMF installed.

    if args.obs_name == 'ORAS4' and args.obs_grid is None:
        print('\n Error: --obs_grid must be used for ORAS4 reanalysis\n', file=sys.stderr)
        parser.print_help()
        return 1

    # Destination grid
    if args.model_name == 'MOM':
        title = 'MOM tripolar t-cell grid'
        model_grid = MomGrid(args.model_hgrid, args.model_vgrid, args.model_mask, title)
    else:
        title = 'NEMO tripolar t-cell grid'
        model_grid = NemoGrid(args.model_hgrid, args.model_vgrid, args.model_mask, title)

    # Source grid
    if args.obs_name == 'ORAS4':
        obs_grid = OrasGrid(args.obs_grid, description='ORAS4')
    else:
        obs_grid = GodasGrid(args.temp_obs_file, description='GODAS')

    # Source-like grids but extended to the whole globe, including maximum
    # depth. The reanalysis grids have limited domain and/or depth.
    if args.obs_name == 'ORAS4':
        global_grid = TripolarGrid(obs_grid, model_grid.z, description='ORAS4')
    else:
        num_lat_points = int(180.0 / obs_grid.dy)
        num_lon_points = int(360.0 / obs_grid.dx)
        description = 'GODAS Equidistant Lat Lon Grid'
        global_grid = RegularGrid(num_lon_points, num_lat_points, model_grid.z,
                                  description=description)

    # Read in temperature and salinity data.
    if args.obs_name == 'ORAS4':
        temp_var = 'thetao'
        salt_var = 'so'
    else:
        temp_var = 'pottmp'
        salt_var = 'salt'

    with nc.Dataset(args.temp_obs_file) as obs:
        otemp = obs.variables[temp_var][args.time_index, :]
        temp_units = obs.variables[temp_var].units
    with nc.Dataset(args.salt_obs_file) as obs:
        osalt = obs.variables[salt_var][args.time_index, :]
        salt_units = obs.variables[salt_var].units
        # Convert salt from kg/kg to g/kg
        if salt_units == 'kg/kg':
            osalt *= 1000.0

    # Regrid obs columns onto model vertical grid.
    print('Vertical regridding/extrapolation ...')
    otemp = regrid_columns(otemp, obs_grid, global_grid)
    osalt = regrid_columns(osalt, obs_grid, global_grid)

    # Now extend obs to cover whole globe
    print('Extending obs to global domain ...')
    gtemp = extend_obs(otemp, obs_grid, global_grid)
    gsalt = extend_obs(osalt, obs_grid, global_grid)

    if args.obs_name == 'GODAS':
        print('Filling Arctic with representational value ...')
        gtemp = fill_arctic(otemp, gtemp, obs_grid, global_grid)
        gsalt = fill_arctic(osalt, gsalt, obs_grid, global_grid)

    # Destination arrays    
    mtemp = np.ndarray((model_grid.num_levels, model_grid.num_lat_points,
                      model_grid.num_lon_points))
    msalt = np.ndarray((model_grid.num_levels, model_grid.num_lat_points,
                      model_grid.num_lon_points))

    # Now apply regridding. We have several options here.
    if args.regrid_tool == 'scipy':
        from mpl_toolkits import basemap

        # Move lons and data onto range 0-360
        # FIXME: assert that lats don't need to be shifted as well.
        mlons, _ = normalise_lons(model_grid.x_t)
        glons, gtemp = normalise_lons(global_grid.x_t, gtemp)
        _, gsalt = normalise_lons(global_grid.x_t, gsalt)

        # Bilinear interpolation over all levels.
        # FIXME: print warning for case where src grid is not regular.
        print('Regridding to model grid')
        for src, dest in [(gtemp, mtemp), (gsalt, msalt)]:
            for l in range(gtemp.shape[0]):
                print('.', end='')
                sys.stdout.flush()
                dest[l, :, :] = basemap.interp(src[l,:,:],
                                               glons[150, :],
                                               global_grid.y_t[:, 150],
                                               mlons, model_grid.y_t,
                                               order=1)
        print('')

    else:
        assert args.regrid_tool == 'esmf'

        # Write the source and destination grids out in SCRIP format.
        global_grid.write_scrip('global_grid_scrip.nc')
        model_grid.write_scrip('model_grid_scrip.nc')

        if not args.regrid_weights:    
            args.regrid_weights = 'regrid_weights.nc'
            # Create regrid weights files using ESMF_RegridWeightGen
            ret = sp.call(['ESMF_RegridWeightGen',
                           '-s', 'global_grid_scrip.nc',
                           '-d', 'model_grid_scrip.nc',
                           '-m', 'bilinear', '-w', args.regrid_weights])
            assert(ret == 0)
            assert(os.path.exists(args.regrid_weights))

        # Creating the remapping weights files is a computationally intensive
        # task. For simplicity call external tool for this.
        with nc.Dataset(args.regrid_weights) as wf:
            n_s = wf.dimensions['n_s'].size
            n_b = wf.dimensions['n_b'].size
            row = wf.variables['row'][:]
            col = wf.variables['col'][:]
            s = wf.variables['S'][:]

        for src, dest in [(gtemp, mtemp), (gsalt, msalt)]:
            for l in range(gtemp.shape[0]):
                start_time = time.time()    
                dest[l, :, :] = apply_weights(src[l, :, :], dest.shape[1:],
                                              n_s, n_b, row, col, s)

    print('Writing out')
    if args.model_name == 'MOM':
        # Apply ocean mask.
        if model_grid.mask is not None:
            mask = np.stack([model_grid.mask] * model_grid.num_levels)
            # Shift mask
            _, mask = normalise_lons(model_grid.x_t, mask)
            temp = np.ma.array(mtemp, mask=mask)
            salt = np.ma.array(msalt, mask=mask)

        write_mom_ic(model_grid, temp, salt, args.output, ''.join(sys.argv))
    else:
        write_nemo_ic(model_grid, mtemp, msalt, args.output, ''.join(sys.argv))

if __name__ == '__main__':
    sys.exit(main())
