#!/usr/bin/env python

from __future__ import print_function

import sys, os
import tempfile
import argparse
import numpy as np
import numba
import copy
import subprocess as sp
import netCDF4 as nc
from numpy import interp
from scipy import ndimage as nd

from .mom_grid import MomGrid
from .mom1_grid import Mom1Grid
from .nemo_grid import NemoGrid
from .regular_grid import RegularGrid
from .tripolar_grid import TripolarGrid
from .godas_grid import GodasGrid
from .oras_grid import OrasGrid
from .woa_grid import WoaGrid


from regridder import util

"""
Create ocean model IC based on reanalysis data.
"""

GODAS_BERING_STRAIGHT = [416, 184]

def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()

def regrid_columns(src_data, src_grid, dest_grid, temp_or_salt):
    """
    Regrid vertical columns of data from src to dest.

    1. Fill in source bathymetry with nearest neighbours.
    2. interpolate and extrapolate down so that dest columns are filled to the
    bottom.
    """

    assert temp_or_salt == 'temp' or temp_or_salt == 'salt'
    assert len(src_data.shape) == 3
    assert src_data.shape[0] == len(src_grid.z)

    # These gets modified.
    src_data_a = np.copy(src_data)
    src_data_b = np.copy(src_data)

    # Create masked array of the correct shape.
    tmp = np.zeros((len(dest_grid.z), src_data.shape[1], src_data.shape[2]))
    new_data = np.ma.array(tmp, mask=np.ones_like(tmp), copy=True)

    # We need to fill in missing values in the source bathymetry. Do this in
    # steps:
    #   1a. Fill with horizontal nearest neighbour so that deep point
    #   have valued based on neighbours at the same depth.
    #   1b. Fill in the land with horizontal nearest neighbour then go down the
    #   water column and for every masked grid cell choose between the
    #   horizontal nearest neighbour (from 1a) and the cell above. The choice
    #   is based on which one will lead to a more stable water column.
    #   2. Regrid, extrapolate to the bottom of the dest columns.

    # 1a. At every level fill everything with nearest neighbour. This
    # effectively removes bathymetry with the new (previously masked) deep
    # points having values based on neighbours at the same depth.
    for lev in range(src_data.shape[0]):
        ind = nd.distance_transform_edt(src_grid.mask[lev, :, :],
                                        return_distances=False,
                                        return_indices=True)
        tmp = src_data_a[lev, :, :]
        tmp = tmp[tuple(ind)]
        src_data_a[lev, :, :] = tmp[:]

    # 1b. First fill in the top level.
    ind = nd.distance_transform_edt(src_grid.mask[0, :, :],
                                   return_distances=False,
                                   return_indices=True)
    tmp = src_data_b[0, :, :]
    tmp = tmp[tuple(ind)]
    src_data_b[0, :, :] = tmp[:]

    # Now fill in all missing values down the columns.
    for lev in range(1, src_data.shape[0]):
        cmask = np.where(src_grid.mask[lev, :, :])
        if temp_or_salt == 'temp':
            best = np.minimum(src_data_b[lev-1, cmask[0], cmask[1]],
                              src_data_a[lev, cmask[0], cmask[1]])
            src_data_b[lev, cmask[0], cmask[1]] = best[:]
        else:
            # salt
            best = np.maximum(src_data_b[lev-1, cmask[0], cmask[1]],
                              src_data_a[lev, cmask[0], cmask[1]])
            src_data_b[lev, cmask[0], cmask[1]] = best[:]


    # Step 2. Iterate through columns and regrid each to the bottom of the
    # destination column.
    for lat in range(src_data.shape[1]):
        for lon in range(src_data.shape[2]):
            # 1d linear interpolation/extrapolation
            new_data[:, lat, lon] = interp(dest_grid.z, src_grid.z,
                                            src_data_b[:, lat, lon])

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


def extend_src_data(src_data, src_grid, global_src_grid, temp_or_salt):
    """
    Extend data to go to full depth and cover global domain.
    """

    # Extend src to go to full depth
    print('Vertical regridding/extrapolation ...')
    src_data = regrid_columns(src_data, src_grid, global_src_grid,
                              temp_or_salt)

    # Now extend src to cover whole globe
    print('Extending obs to global domain ...')
    global_src_data = extend_to_global(src_data, src_grid, global_src_grid)

    # Possibly fill in Arctic
    if 'GODAS' in src_grid.description:
        print('Filling Arctic with representational value ...')
        global_src_data = fill_arctic(src_data, global_src_data, src_grid,
                                      global_src_grid)

    return global_src_data

def apply_weights(src, dest_shape, n_s, n_b, row, col, s):
    """
    Apply ESMF regirdding weights.
    """

    dest = np.ndarray(dest_shape).flatten()
    dest[:] = 0.0
    src = src.flatten()

    for i in range(1, n_s):
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

def smooth_all(data):

    sigma = (2, 5, 5)

    new_data = np.copy(data)
    new_data[:, :, :] = nd.filters.gaussian_filter(data[:, :, :], sigma)
    return new_data

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

def is_var_temp_or_salt(src_var, dest_var):

    for v in [src_var.lower(), dest_var.lower()]:
        if v == 'salt' or v == 'vosaline' or v == 'practical_salinity':
            return 'salt'
        if v == 'temp' or v == 'votemper' or v == 'pottmp' or v == 'potential_temperature' or v == 'conservative_temperature':
            return 'temp'

def check_src_data_ranges(src_data, temp_or_salt):

    if temp_or_salt == 'temp':
        assert np.max(src_data) < 320
        assert np.min(src_data) >= -10

    if temp_or_salt == 'salt':
        assert np.max(src_data) >= 0
        assert np.max(src_data) < 60


def do_regridding(src_name, src_hgrids, src_vgrid, src_data_file, src_var,
                  dest_name, dest_hgrid, dest_vgrid, dest_data_file, dest_var,
                  dest_mask=None, month=None, regrid_weights=None, use_mpi=False,
                  write_ic=False):

    if not check_dependencies(use_mpi):
        return None

    filenames = list(src_hgrids) + [src_vgrid, src_data_file, dest_hgrid, dest_vgrid]
    if dest_mask is not None:
        filenames.append(dest_mask)

    if check_files(filenames):
        return None;

    temp_or_salt = is_var_temp_or_salt(src_var, dest_var)

    if dest_name == 'MOM':
        title = 'MOM tripolar 0.25 degree t-cell grid'
        dest_grid = MomGrid(dest_hgrid, dest_vgrid, dest_mask, title)
    elif dest_name == 'MOM1':
        title = 'MOM tripolar 1 degree t-cell grid'
        dest_grid = Mom1Grid(dest_hgrid, dest_vgrid, dest_mask, title)
    else:
        title = 'NEMO tripolar t-cell grid'
        dest_grid = NemoGrid(dest_hgrid, dest_vgrid, dest_mask, title)

    # Source grid
    if src_name == 'ORAS4':
        assert len(src_hgrids) >= 1 and len(src_hgrids) <= 3
        src_grid = OrasGrid(src_hgrids, description='ORAS4')
    elif src_name == 'GODAS':
        assert len(src_hgrids) == 1
        src_grid = GodasGrid(src_hgrids[0], description='GODAS')
    elif src_name == 'WOA':
        assert len(src_hgrids) == 1
        src_grid = WoaGrid(src_hgrids[0], description='WOA')
    else:
        print('\n Error: invalid source name: {}.\n'.format(src_name),
            file=sys.stderr)
        return None
        
    # The source grids need to be extended to the whole globe, including
    # maximum depth. The reanalysis grids have limited domain and/or depth.
    if src_name == 'ORAS4':
        global_src_grid = TripolarGrid(src_grid, dest_grid.z,
                                       description='ORAS4')
    elif src_name == 'GODAS': 
        num_lat_points = int(180.0 / src_grid.dy)
        num_lon_points = int(360.0 / src_grid.dx)
        description = 'GODAS Equidistant Lat Lon Grid'
        global_src_grid = RegularGrid(num_lon_points, num_lat_points,
                                      dest_grid.z, description=description)
    elif src_name == 'WOA':
        # WOA is global, just needs to have the depth fixed.
        global_src_grid = copy.deepcopy(src_grid)
        global_src_grid.z = dest_grid.z
        global_src_grid.mask = np.zeros((global_src_grid.num_levels,
                                         global_src_grid.num_lat_points,
                                         global_src_grid.num_lon_points), dtype='int')


    # Write the source and destination grids out in SCRIP format. We override
    # the masks, we want to cover everything.
    _, global_src_grid_scrip = tempfile.mkstemp(suffix='.nc')
    global_src_grid.write_scrip(global_src_grid_scrip,
            mask=np.zeros_like(global_src_grid.mask, dtype=int))
    _, dest_grid_scrip = tempfile.mkstemp(suffix='.nc')
    dest_grid.write_scrip(dest_grid_scrip,
            mask=np.zeros_like(dest_grid.mask, dtype=int))

    print('global_src_grid_scrip {}'.format(global_src_grid_scrip))
    print('dest_grid_scrip {}'.format(dest_grid_scrip))

    # Creating the remapping weights files is a computationally intensive
    # task. For simplicity call an external tool for this.
    if regrid_weights is None or not os.path.exists(regrid_weights):
        if regrid_weights is None:
            _, regrid_weights = tempfile.mkstemp(suffix='.nc')
        mpi = []
        if use_mpi:
            mpi = ['mpirun', '-n', '8']

        try:
            sp.check_output(mpi + ['ESMF_RegridWeightGen',
                                   '-s', global_src_grid_scrip,
                                   '-d', dest_grid_scrip,
                                   '-m', 'bilinear', '-w', regrid_weights])
        except sp.CalledProcessError as e:
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
    time_origin = util.get_time_origin(src_data_file)
    if not os.path.exists(dest_data_file):
        if 'MOM' in dest_name:
            util.create_mom_output(dest_grid, dest_data_file, time_origin,
                              ''.join(sys.argv))
        else:
            util.create_nemo_output(dest_grid, dest_data_file, time_origin,
                               ''.join(sys.argv))

    # Do regridding on each time point.
    f = nc.Dataset(src_data_file)
    src_var = f.variables[src_var]
    if src_name == 'ORAS4':
        # FIXME: ORAS4 hack to deal with duplicate rows/columns
        src_data = src_grid.fix_data_shape(src_var[:])
        # Also add mask to data
        src_data = src_grid.apply_grid_mask(src_data)
    else:
        src_data = src_var[:]

        if src_name == 'GODAS':
            # Give the grid a new mask, this is because there are tiny differences
            # in the mask for each time point of GODAS data.
            new_mask = np.sum(src_data.mask[:, :, :, :], axis=0)
            new_mask[np.where(new_mask > 1)] = 1
            src_grid.set_mask(new_mask)

    check_src_data_ranges(src_data, temp_or_salt)

    if month is not None:
        time_idxs = [month - 1]
    else:
        time_idxs = range(src_var.shape[0])
    time_points = f.variables['time'][time_idxs]
    for t_idx, t_pt in zip(time_idxs, time_points):
        ext_src_data = extend_src_data(src_data[t_idx, :], src_grid, global_src_grid,
                                        temp_or_salt)
        dest_data = regrid(regrid_weights, ext_src_data, dest_grid)

        # FIXME: issue with regridding in bottom left corner of grid. This is
        # masked in any case.
        dest_data[np.where(dest_data <= np.min(ext_src_data))] = np.min(ext_src_data)

        # FIXME: run a smoother to remove sharp edges associated with missing data.
        if (src_name == 'GODAS' or src_name == 'ORAS4') and \
            dest_name == 'MOM' and write_ic:
            dest_data = smooth_all(dest_data)

        # Write out
        try:
            units = src_var.units
        except AttributeError:
            units = ''
        try:
            long_name = src_var.long_name
        except AttributeError:
            long_name = ''

        # Input file has units in hours, convert to days.
        if 'hours since' in f.variables['time'].units:
            t_pt = int(t_pt / 24.)

        if 'MOM' in dest_name:
            # Apply ocean mask.
            if dest_grid.mask is not None:
                mask = np.stack([dest_grid.mask] * dest_grid.num_levels)
                dest_data = np.ma.array(dest_data, mask=mask)
            util.write_mom_output_at_time(dest_data_file, dest_var, long_name,
                                     units, dest_data, t_idx, t_pt, write_ic)
        else:
            util.write_nemo_output_at_time(dest_data_file, dest_var, long_name,
                                      units, dest_data, t_idx, t_pt, write_ic)

    f.close()
    for f in [global_src_grid_scrip, global_src_grid_scrip + '_test',
              dest_grid_scrip, dest_grid_scrip + '_test']:
        try:
            os.remove(f)
            pass
        except OSError:
            pass

    return regrid_weights

def check_files(filenames):

    for filename in filenames:
        if not os.path.exists(filename):
            print("File {} not found, please check it's location.".format(filename),
                    file=sys.stderr)
            return 1

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

    assert args.dest_name == 'MOM' or args.dest_name == 'MOM1' or \
        args.dest_name == 'NEMO'
    assert args.src_name == 'GODAS' or args.src_name == 'ORAS4'

    if os.path.exists(args.dest_data_file) and not args.append:
        print("Output file {} already exists, ".format(args.dest_data_file) + \
              "please move use the --append option.", file=sys.stderr)
        return 1

    ret = do_regridding(args.src_name, (args.src_hgrid,), args.src_vgrid,
                        args.src_data_file, args.src_var,
                        args.dest_name, args.dest_hgrid, args.dest_vgrid,
                        args.dest_data_file, args.dest_var,
                        args.dest_mask, args.month, args.regrid_weights,
                        args.use_mpi)
    return ret is None

if __name__ == '__main__':
    sys.exit(main())
