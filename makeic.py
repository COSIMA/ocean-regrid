#!/usr/bin/env python

from __future__ import print_function

import sys, os
import argparse
import numpy as np
import netCDF4 as nc

from scipy import interp
from scipy import ndimage as nd
import mpl_toolkits
from mpl_toolkits.basemap import Basemap

from mom_grid import MomGrid
from grid import Grid
from latlon_grid import LatLonGrid

"""
"""

def regrid_columns(data, src_z, dest_z, plot_results=False):
    """
    Regrid vertical columns of data from src_z to dest_z.
    """

    assert len(data.shape) == 3
    assert data.shape[0] == len(src_z)

    # This gets modified.
    data = np.ma.copy(data)

    # Create masked array of the correct shape.
    tmp = np.zeros((len(dest_z), data.shape[1], data.shape[2]))
    new_data = np.ma.array(tmp, mask=np.ones_like(tmp))

    # Iterate through columns and regrid each.
    for lat in range(data.shape[1]):
        for lon in range(data.shape[2]):
            if all(data[:, lat, lon].mask):
                continue
            if any(data[:, lat, lon].mask):
                # Masked values expected to be at depth. Find these and fill
                # with nearest neighbour.
                for d in range(data.shape[0]-2, 0, -1):
                    if (not data.mask[d, lat, lon]) and \
                            data.mask[d+1, lat, lon]:
                        data[d:, lat, lon] = data[d, lat, lon]

            # 1d linear interpolation/extrapolation
            new_data[:, lat, lon] = interp(dest_z, src_z, data[:, lat, lon])

    return new_data


def extend_obs(obs_grid, global_grid, var):
    """
    Use nearest neighbour to extend obs over the whole globe, including land.
    """

    # Create masked array of the correct shape.
    tmp = np.zeros((global_grid.num_levels, global_grid.num_lat_points,
                    global_grid.num_lon_points))
    new_data = np.ma.array(tmp, mask=np.ones_like(tmp))

    # Drop obs data into new grid at correct location
    def find_nearest_index(array, value):
        return (np.abs(array - value)).argmin()
    lat_min_idx = find_nearest_index(global_grid.y_t[:, 0], obs_grid.y_t[0, 0])
    lat_max_idx = find_nearest_index(global_grid.y_t[:, 0], obs_grid.y_t[-1, 0])
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
    parser.add_argument('ocean_hgrid', help='Ocean horizontal grid spec file.')
    parser.add_argument('ocean_vgrid', help='Ocean vertical grid spec file.')
    parser.add_argument('ocean_mask', help='Ocean land-sea mask file.')
    parser.add_argument('temp_obs_file', help='Temp from GODAS obs reanalysis.')
    parser.add_argument('salt_obs_file', help='Salt from GODAS obs reanalysis.')
    parser.add_argument('--output', default='ocean_ic.nc',
                        help='Name of the output file.')
    parser.add_argument('--temp_var', default='temp',
                        help='Name of the obs temperature variable')
    parser.add_argument('--salt_var', default='salt',
                        help='Name of the obs salt variable')
    parser.add_argument('--time_index', default=0,
                        help='The time index of the data to use.')
    parser.add_argument('--x_var', default='lon',
                        help='Name of the obs longitude variable')
    parser.add_argument('--y_var', default='lat',
                        help='Name of the obs latitude variable')
    parser.add_argument('--z_var', default='level',
                        help='Name of the obs vertical variable')
    parser.add_argument('--use_esmf', default=False,
                        help='Do regridding with ESMF, otherwise use basemap.')
    args = parser.parse_args()

    # Destination grid
    title = 'MOM tripolar t-cell grid'
    mom_grid = MomGrid(args.ocean_hgrid, args.ocean_vgrid, args.ocean_mask, title)

    # Read in obs file.
    with nc.Dataset(args.temp_obs_file) as obs:
        temp = obs.variables[args.temp_var][args.time_index, :]
        lons = obs.variables[args.x_var][:]
        lats = obs.variables[args.y_var][:]
        z = obs.variables[args.z_var][:]

    # Source grid.
    title = '{}x{} Cylindrical Equidistant Projection Grid'.format(len(lons), len(lats))
    obs_grid = Grid(lons, lats, z, temp.mask[0, :, :], title)

    # Regrid obs columns onto model vertical grid.
    temp = regrid_columns(temp, z, mom_grid.z, plot_results=True)

    # Source-like grid but extended to the whole globe.
    num_lat_points = 180.0 / abs(lats[1] - lats[0])
    num_lon_points = 360.0 / abs(lons[1] - lons[0])
    title = '{}x{} Equidistant Lat Lon Grid'
    mask = np.zeros((num_lat_points, num_lon_points))
    global_grid = LatLonGrid(num_lon_points, num_lat_points, mom_grid.z,
                              mask, title)
    # Now extend obs to cover whole globe
    gtemp = extend_obs(obs_grid, global_grid, temp)

    # Move lons from -280 to 80 to 0 to 360 for interpolation step.
    # Use basemap shift grid.
    x_t = np.copy(mom_grid.x_t)
    x_t[mom_grid.x_t < 0] = mom_grid.x_t[mom_grid.x_t < 0] + 360

    # Bilinear interpolation
    temp_bi = mpl_toolkits.basemap.interp(gtemp[0,:,:],
                                          global_grid.x_t[0, :],
                                          global_grid.y_t[:, 0],
                                          x_t, mom_grid.y_t, order=1)
    # Apply ocean mask.
    temp = np.ma.array(temp_bi, mask=mom_grid.mask)

    import pdb
    pdb.set_trace()


if __name__ == '__main__':
    sys.exit(main())
