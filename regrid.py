#!/usr/bin/env python

from __future__ import print_function

import sys, os
import tempfile
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

from file_util import create_mom_output, write_mom_output_at_time
from file_util import create_nemo_output, write_nemo_output_at_time
from util import normalise_lons, get_time_origin

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
    new_data = np.ma.array(tmp, mask=np.ones_like(tmp), copy=True)

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


def extend_to_global(var, src_grid, global_grid, arctic_filler=None):
    """
    Use nearest neighbour to extend obs/source over the whole globe, including land.
    """

    # Create masked array of the correct shape.
    tmp = np.zeros((global_grid.num_levels, global_grid.num_lat_points,
                    global_grid.num_lon_points))
    new_data = np.ma.array(tmp, mask=global_grid.mask, copy=True)
    # Mask everything by default. The code below fills masked values with
    # nearest neighbour.
    new_data.mask[:] = True

    # Drop obs data into new grid at correct location
    lat_min_idx = find_nearest_index(global_grid.y_t[:, 0], np.min(src_grid.y_t[:]))
    if np.max(global_grid.y_t[:]) <= np.max(src_grid.y_t[:]):
        new_data[:, lat_min_idx:, :] = var[:]
    else:
        lat_max_idx = find_nearest_index(global_grid.y_t[:, 0], np.max(src_grid.y_t[:]))
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

def fill_arctic(src_data, global_data, src_grid, global_grid):
    """
    In the case of GODAS, there is no data past 65N and the nearest neighbour
    approach is not ideal because the low salinity of the Baltic gets
    propogated into the Arctic. So instead fill the Arctic with a
    'representative value' taken from the Bering Straight.
    """

    assert 'GODAS' in src_grid.description

    GODAS_ARCTIC_REPRESENTATIVE_VALUE = GODAS_BERING_STRAIGHT
    filler = src_data[:, GODAS_ARCTIC_REPRESENTATIVE_VALUE[0],
                         GODAS_ARCTIC_REPRESENTATIVE_VALUE[1]]

    arctic_idx = find_nearest_index(global_grid.y_t[:, 0],
                                    np.max(src_grid.y_t[:, 0]))
    arctic_idx -= 1
    sh = global_data[:, arctic_idx:, :].shape

    filler = np.stack([filler[:]] * sh[1], axis=1)
    filler = np.stack([filler[:]] * sh[2], axis=2)
    global_data[:, arctic_idx:, :] = filler[:]

    return global_data


def extend_src_data(src_data, src_grid, global_src_grid):
    """
    Extend data to go to full depth and cover global domain.
    """

    # Extend src to go to full depth
    print('Vertical regridding/extrapolation ...')
    src_data = regrid_columns(src_data, src_grid, global_src_grid)

    # Now extend src to cover whole globe
    print('Extending obs to global domain ...')
    global_src_data = extend_to_global(src_data, src_grid, global_src_grid)

    # Possibly fill in Arctic
    if 'GODAS' in src_grid.description:
        print('Filling Arctic with representational value ...')
        global_src_data = fill_arctic(src_data, global_src_data, src_grid,
                                      global_src_grid)
    return global_src_data


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


def regrid(regrid_weights, src_data, dest_grid): 
    """
    Regrid a single time index of data.
    """

    print('Horizontal regridding ...')
    # Destination arrays
    dest_data = np.ndarray((dest_grid.num_levels, dest_grid.num_lat_points,
                            dest_grid.num_lon_points))

    with nc.Dataset(regrid_weights) as wf:
        n_s = wf.dimensions['n_s'].size
        n_b = wf.dimensions['n_b'].size
        row = wf.variables['row'][:]
        col = wf.variables['col'][:]
        s = wf.variables['S'][:]

    for l in range(src_data.shape[0]):
        dest_data[l, :, :] = apply_weights(src_data[l, :, :], dest_data.shape[1:],
                                           n_s, n_b, row, col, s)

    return dest_data

def check_dependencies(use_mpi):

    ret = sp.call(['which', 'ESMF_RegridWeightGen'])
    if ret:
        print('\n Error: regrid.py program depends on on ESMF_RegridWeightGen which is not installed.\n',
               file=sys.stderr)
        return False

    if use_mpi:
        ret = sp.call(['which', 'mpirun'])
        if ret:
            print('\n Error: mpirun must be installed when the --use_mpi flag is used.\n',
                   file=sys.stderr)
            return False

    return True

def do_regridding(src_name, src_hgrid, src_vgrid, src_data_file, src_var,
                  dest_name, dest_hgrid, dest_vgrid, dest_data_file, dest_var,
                  dest_mask=None, month=None, regrid_weights=None, use_mpi=False):

    if not check_dependencies(use_mpi):
        return None

    # Destination grid
    if dest_name == 'MOM':
        if dest_mask is None:
            print('\n Error: please provide a --dest_mask when regridding to MOM.\n',
                  file=sys.stderr)
            return None
        title = 'MOM tripolar t-cell grid'
        dest_grid = MomGrid(dest_hgrid, dest_vgrid, dest_mask, title)
    else:
        title = 'NEMO tripolar t-cell grid'
        dest_grid = NemoGrid(dest_hgrid, dest_vgrid, dest_mask, title)

    # Source grid
    if src_name == 'ORAS4':
        src_grid = OrasGrid(src_hgrid, description='ORAS4')
    else:
        src_grid = GodasGrid(src_hgrid, description='GODAS')

    # An extra, source-like grids but extended to the whole globe, including
    # maximum depth. The reanalysis grids have limited domain and/or depth.
    if src_name == 'ORAS4':
        global_src_grid = TripolarGrid(src_grid, dest_grid.z,
                                       description='ORAS4')
    else:
        num_lat_points = int(180.0 / src_grid.dy)
        num_lon_points = int(360.0 / src_grid.dx)
        description = 'GODAS Equidistant Lat Lon Grid'
        global_src_grid = RegularGrid(num_lon_points, num_lat_points,
                                      dest_grid.z, description=description)

    # Write the source and destination grids out in SCRIP format. We override
    # the src mask because we want to cover everything.
    _, global_src_grid_scrip = tempfile.mkstemp(suffix='.nc')
    unmask_all = np.zeros_like(global_src_grid.mask, dtype=int)
    global_src_grid.write_scrip(global_src_grid_scrip, mask=unmask_all)
    _, dest_grid_scrip = tempfile.mkstemp(suffix='.nc')
    dest_grid.write_scrip(dest_grid_scrip)

    # Creating the remapping weights files is a computationally intensive
    # task. For simplicity call an external tool for this.
    if regrid_weights is None or not os.path.exists(regrid_weights):
        if regrid_weights is None:
            regrid_weights = 'regrid_weights.nc'
        mpi = []
        if use_mpi:
            mpi = ['mpirun', '-n', '8']

        try:
            sp.check_output(mpi + ['ESMF_RegridWeightGen',
                                   '-s', global_src_grid_scrip,
                                   '-d', dest_grid_scrip,
                                   '-m', 'bilinear', '-w', regrid_weights])
        except sp.CalledProcessError, e:
            print("Error: ESMF_RegridWeightGen failed return code {}".format(e.returncode),
                  file=sys.stderr)
            print(e.output, file=sys.stderr)
            log = 'PET0.RegridWeightGen.Log'
            if os.path.exists(log):
                print('Contents of {}:'.format(log), file=sys.stderr)
                with open(log) as f:
                    print(f.read(), file=sys.stderr)
            return None

        assert(os.path.exists(regrid_weights))

    # Create output file
    time_origin = get_time_origin(src_data_file)
    if not os.path.exists(dest_data_file):
        if dest_name == 'MOM':
            create_mom_output(dest_grid, dest_data_file, time_origin,
                              ''.join(sys.argv))
        else:
            create_nemo_output(dest_grid, dest_data_file, time_origin,
                               ''.join(sys.argv))

    # Do regridding on each time point.
    f = nc.Dataset(src_data_file)
    src_var = f.variables[src_var]

    if month is not None:
        time_idxs = [month - 1]
    else:
        time_idxs = range(src_var.shape[0])
    time_points = f.variables['time'][time_idxs]

    for t_idx, t_pt in zip(time_idxs, time_points):
        src_data = src_var[t_idx, :]
        src_data = extend_src_data(src_data, src_grid, global_src_grid)
        dest_data = regrid(regrid_weights, src_data, dest_grid)

        # Write out
        try:
            units = src_var.units
        except AttributeError:
            units = ''
        try:
            long_name = src_var.long_name
        except AttributeError:
            long_name = ''

        if dest_name == 'MOM':
            # Apply ocean mask.
            if dest_grid.mask is not None:
                mask = np.stack([dest_grid.mask] * dest_grid.num_levels)
                dest_data = np.ma.array(dest_data, mask=mask)
            write_mom_output_at_time(dest_data_file, dest_var, long_name,
                                     units, dest_data, t_idx, t_pt)
        else:
            write_nemo_output_at_time(dest_data_file, dest_var, long_name,
                                      units, dest_data, t_idx, t_pt)

    f.close()
    return regrid_weights

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('src_name', help="""
                        Name of src data/grid, must be GODAS or ORAS4""")
    parser.add_argument('src_hgrid', help='Input horizontal grid spec file.')
    parser.add_argument('src_vgrid', help='Input vertical grid spec file.')
    parser.add_argument('src_data_file', help='File containing reanalysis dataset.')
    parser.add_argument('src_var', help='Name of input variable to regrid.')

    parser.add_argument('dest_name', help="""
                        Name of dest data/grid, must be MOM or NEMO""")
    parser.add_argument('dest_hgrid', help='Output horizontal grid spec file.')
    parser.add_argument('dest_vgrid', help='Output vertical grid spec file.')
    parser.add_argument('dest_data_file', help='Name of the destination/output file.')
    parser.add_argument('dest_var', help='Name of the destination/output variable.')
    parser.add_argument('--dest_mask', default=None, help='Destination land-sea mask file.')
    parser.add_argument('--month', default=None, help='Regrid a single month. Default is all.')

    parser.add_argument('--regrid_weights', default=None,
                        help="""
                        The name of the regridding weights file. Will be created if it doesn't exist
                        """)
    parser.add_argument('--use_mpi', action='store_true', default=False,
                        help="""Use MPI to when calculating the regridding weights.
                               This will speed up the calculation considerably.""")
    parser.add_argument('--append', default=False, action='store_true',
                        help='Append to destination file.')
    args = parser.parse_args()

    assert args.dest_name == 'MOM' or args.dest_name == 'NEMO'
    assert args.src_name == 'GODAS' or args.src_name == 'ORAS4'

    if os.path.exists(args.dest_data_file) and not args.append:
        print("Output file {} already exists, ".format(args.dest_data_file) + \
              "please move use the --append option.", file=sys.stderr)
        return 1

    ret = do_regridding(args.src_name, args.src_hgrid, args.src_vgrid,
                        args.src_data_file, args.src_var,
                        args.dest_name, args.dest_hgrid, args.dest_vgrid,
                        args.dest_data_file, args.dest_var,
                        args.dest_mask, args.month, args.regrid_weights,
                        args.use_mpi)
    return ret is None

if __name__ == '__main__':
    sys.exit(main())
