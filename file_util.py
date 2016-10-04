
import sys
import netCDF4 as nc
import numpy as np
from mom_grid import MomGrid

def row_idx_largest_lon(lons):
    """
    The row index with the largest lon.
    """
    r, _  = np.unravel_index(np.argmax(lons), lons.shape)

    return r

def col_idx_largest_lat(lats):
    """
    The col index with the largest lat.
    """
    _, c  = np.unravel_index(np.argmax(lats), lats.shape)

    return c


def create_mom_output(ocean_grid, filename, start_date, history):

    f = nc.Dataset(filename, 'w')

    f.createDimension('GRID_X_T', ocean_grid.num_lon_points)
    f.createDimension('GRID_Y_T', ocean_grid.num_lat_points)
    f.createDimension('ZT', ocean_grid.num_levels)
    f.createDimension('time')

    lons = f.createVariable('GRID_X_T', 'f8', ('GRID_X_T'))
    lons.long_name = 'Nominal Longitude of T-cell center'
    lons.units = 'degree_east'
    lons.modulo = 360.
    lons.point_spacing = 'even'
    lons.axis = 'X'
    # MOM needs this to be a single dimension
    row = row_idx_largest_lon(ocean_grid.x_t[:])
    lons[:] = ocean_grid.x_t[row, :]

    lats = f.createVariable('GRID_Y_T', 'f8', ('GRID_Y_T'))
    lats.long_name = 'Nominal Latitude of T-cell center'
    lats.units = 'degree_north'
    lats.point_spacing = 'uneven'
    lats.axis = 'Y'
    # MOM needs this to be a single dimension
    col = col_idx_largest_lat(ocean_grid.y_t[:])
    lats[:] = ocean_grid.y_t[:, col]

    zt = f.createVariable('ZT', 'f8', ('ZT'))
    zt.long_name = 'zt'
    zt.units = 'meters'
    zt.positive = 'downdown'
    zt.point_spacing = 'uneven'
    zt.axis = 'Z'
    zt[:] = ocean_grid.z[:]

    time = f.createVariable('time', 'f8', ('time'))
    time.long_name = 'time'
    time.units = "days since {}-{}-{} 00:00:00".format(str(start_date.year).zfill(4),
                                                       str(start_date.month).zfill(2),
                                                       str(start_date.day).zfill(2))
    time.cartesian_axis = "T"
    time.calendar_type = "GREGORIAN"
    time.calendar = "GREGORIAN"

    f.close()

def write_mom_output_at_time(filename, var_name, var_longname, var_units,
                             var_data, time_idx, time_pt, write_ic=False):

    with nc.Dataset(filename, 'r+') as f:
        if not f.variables.has_key(var_name):
            var = f.createVariable(var_name, 'f8',
                                   ('time', 'ZT', 'GRID_Y_T', 'GRID_X_T'),
                                   fill_value=-1.e+34)
            var.missing_value = -1.e+34
            var.long_name = var_longname
            var.units = var_units

        var = f.variables[var_name]

        if write_ic:
            var[0, :] = var_data[:]
            f.variables['time'][0] = time_pt
        else:
            var[time_idx, :] = var_data[:]
            f.variables['time'][time_idx] = time_pt


def create_nemo_output(ocean_grid, filename, start_date, history):

    f = nc.Dataset(filename, 'w')

    f.createDimension('y', ocean_grid.num_lat_points)
    f.createDimension('x', ocean_grid.num_lon_points)
    f.createDimension('z', ocean_grid.num_levels)
    f.createDimension('time_counter')

    lats = f.createVariable('nav_lat', 'f8', ('y', 'x'))
    lats[:] = ocean_grid.y_t[:]

    lons = f.createVariable('nav_lon', 'f8', ('y', 'x'))
    lons[:] = ocean_grid.x_t[:]

    depth = f.createVariable('depth', 'f8', ('z'))
    depth[:] = ocean_grid.z[:]

    time = f.createVariable('time_counter', 'f8', ('time_counter'))
    time.long_name = 'time'
    time.units = "days since {}-{}-{} 00:00:00".format(str(start_date.year).zfill(4),
                                                       str(start_date.month).zfill(2),
                                                       str(start_date.day).zfill(2))
    time.cartesian_axis = "T"

    f.close()

def write_nemo_output_at_time(filename, var_name, var_longname, var_units,
                              var_data, time_idx, time_pt, write_ic=False):

    with nc.Dataset(filename, 'r+') as f:
        if not f.variables.has_key(var_name):
            var = f.createVariable(var_name, 'f8', ('time_counter', 'z', 'y', 'x'))
            var.long_name = var_longname
            var.units = var_units

        var = f.variables[var_name]
        if write_ic:
            var[0, :] = var_data[:]
            f.variables['time_counter'][0] = time_pt
        else:
            var[time_idx, :] = var_data[:]
            f.variables['time_counter'][time_idx] = time_pt
