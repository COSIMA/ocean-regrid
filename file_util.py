
import netCDF4 as nc
from mom_grid import MomGrid

def create_mom_output(ocean_grid, filename,
                      var_name, var_longname, var_units, history):

    f = nc.Dataset(filename, 'w')

    f.createDimension('GRID_Y_T', ocean_grid.num_lat_points)
    f.createDimension('GRID_X_T', ocean_grid.num_lon_points)
    f.createDimension('ZT', ocean_grid.num_levels)
    f.createDimension('time')

    lats = f.createVariable('GRID_Y_T', 'f8', ('GRID_Y_T'))
    lats.long_name = 'Nominal Latitude of T-cell center'
    lats.units = 'degree_north'
    lats.point_spacing = 'uneven'
    lats.axis = 'Y'
    # FIXME
    lats[:] = ocean_grid.y_t[:, 0]

    lons = f.createVariable('GRID_X_T', 'f8', ('GRID_X_T'))
    lons.long_name = 'Nominal Longitude of T-cell center'
    lons.units = 'degree_east'
    lons.modulo = 360.
    lons.point_spacing = 'even'
    lons.axis = 'X'
    # FIXME
    lons[:] = ocean_grid.x_t[0, :]

    zt = f.createVariable('ZT', 'f8', ('ZT'))
    zt.long_name = 'zt'
    zt.units = 'meters'
    zt.positive = 'downdown'
    zt.point_spacing = 'uneven'
    zt.axis = 'Z'
    # FIXME
    zt[:] = ocean_grid.z[:]

    time = f.createVariable('time', 'f8', ('time'))
    time.long_name = 'time'
    time.units = "days since 0001-01-01 00:00:00"
    time.cartesian_axis = "T"
    time.calendar_type = "GREGORIAN"
    time.calendar = "GREGORIAN"

    temp = f.createVariable(var_name, 'f8', ('time', 'ZT', 'GRID_Y_T', 'GRID_X_T'), fill_value=-1.e+34)
    temp.missing_value = -1.e+34
    temp.long_name = var_longname
    temp.units = var_units
    temp.history = history

    f.close()

def write_mom_output_at_time(filename, var_name, var_data, time_idx):

    with nc.Dataset(filename, 'r+') as f:
        var = f.variables[var_name]
        var[time_idx, :] = var_data[:]


def create_nemo_output(ocean_grid, filename,
                            var_name, var_longname, var_units, history):

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
    time.units = "days since 0001-01-01 00:00:00"
    time.cartesian_axis = "T"

    f.createVariable(var_name, 'f8', ('time_counter', 'z', 'y', 'x'))

    f.close()

def write_nemo_output_at_time(filename, var_name, var_data, time_idx):

    with nc.Dataset(filename, 'r+') as f:
        var = f.variables[var_name]
        var[time_idx, :] = var_data[:]
