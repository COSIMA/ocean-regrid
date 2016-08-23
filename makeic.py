#!/usr/bin/env python

from __future__ import print_function

import sys, os
import argparse
import numpy as np
import netCDF4 as nc
from scipy import interp
from scipy import ndimage as nd
from mpl_toolkits import basemap

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
            if all(src_grid.mask[:, lat, lon]):
                continue
            if any(src_grid.mask[:, lat, lon]):
                # Masked values expected to be at depth. Find these and fill
                # with nearest neighbour.
                for d in range(data.shape[0]-2, 0, -1):
                    if (not src_grid.mask[d, lat, lon]) and \
                            src_grid.mask[d+1, lat, lon]:
                        data[d:, lat, lon] = data[d, lat, lon]

            # 1d linear interpolation/extrapolation
            new_data[:, lat, lon] = interp(dest_grid.z, src_grid.z, data[:, lat, lon])

    return new_data


def extend_obs(var, obs_grid, global_grid):
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

    args = parser.parse_args()

    assert args.model_name == 'MOM' or args.model_name == 'NEMO'
    assert args.obs_name == 'GODAS' or args.obs_name == 'ORAS4'

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
        obs_grid = OrasGrid(args.obs_grid)
    else:
        obs_grid = GodasGrid(args.temp_obs_file)

    # Source-like grids but extended to the whole globe, including maximum
    # depth. The reanalysis grids have limited domain and/or depth.
    if args.obs_name == 'ORAS4':
        global_grid = TripolarGrid(obs_grid, model_grid.z)
    else:
        num_lat_points = int(180.0 / obs_grid.dy)
        num_lon_points = int(360.0 / obs_grid.dx)
        description = '{}x{} Equidistant Lat Lon Grid'
        global_grid = RegularGrid(num_lon_points, num_lat_points, model_grid.z,
                                  description=description)

    # Read in temperature and salinity data.
    if args.obs_name == 'ORAS4':
        temp_var = 'thetao'
        salt_var = 'so'
    else:
        temp_var = 'temp'
        salt_var = 'salt'

    with nc.Dataset(args.temp_obs_file) as obs:
        temp = obs.variables[temp_var][args.time_index, :]
        temp_units = obs.variables[temp_var].units
    with nc.Dataset(args.salt_obs_file) as obs:
        salt = obs.variables[salt_var][args.time_index, :]
        salt_units = obs.variables[salt_var].units
        # Convert salt from kg/kg to g/kg
        if salt_units == 'kg/kg':
            salt *= 1000.0

    # Regrid obs columns onto model vertical grid.
    print('Vertical regridding/extrapolation ...')
    temp = regrid_columns(temp, obs_grid, global_grid)
    salt = regrid_columns(salt, obs_grid, global_grid)

    # Now extend obs to cover whole globe
    print('Extending obs to global domain ...')
    gtemp = extend_obs(temp, obs_grid, global_grid)
    gsalt = extend_obs(salt, obs_grid, global_grid)

    # Move lons and data onto range 0-360
    # FIXME: assert that lats don't need to be shifted as well.
    mlons, _ = normalise_lons(model_grid.x_t)
    glons, gtemp = normalise_lons(global_grid.x_t, gtemp)
    _, gsalt = normalise_lons(global_grid.x_t, gsalt)

    # Bilinear interpolation over all levels.
    # FIXME: should use sosie or ESMF for this in the case where obs grid is
    # not rectilinear.
    print('Regridding to model grid')
    mtemp = np.ndarray((model_grid.num_levels, model_grid.num_lat_points,
                      model_grid.num_lon_points))
    msalt = np.ndarray((model_grid.num_levels, model_grid.num_lat_points,
                      model_grid.num_lon_points))
    for src, dest in [(gtemp, mtemp), (gsalt, msalt)]:
        for l in range(gtemp.shape[0]):
            print('.', end='')
            sys.stdout.flush()
            dest[l, :, :] = basemap.interp(src[l,:,:],
                                           glons[0, :],
                                           global_grid.y_t[:, 0],
                                           mlons, model_grid.y_t,
                                           order=1)
    print('')

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
        write_nemo_ic(model_grid, temp, salt, args.output, ''.join(sys.argv))

if __name__ == '__main__':
    sys.exit(main())
