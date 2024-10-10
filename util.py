
import numpy as np
import re
import datetime as dt
import netCDF4 as nc

def normalise_lons(lons, data=None):
    """
    Normalise longitudes to 0-360 deg. Perform the same transformation on data.
    """

    # Remove -ves
    new_lons = np.copy(lons)
    new_lons[lons < 0] = lons[lons < 0] + 360
    lons_copy = np.copy(new_lons)

    if data is not None:
        new_data = np.copy(data)
    else:
        new_data = None

    # Use changes in 2nd derivative to find jumps. Then offset (the +2) to get
    # element directly after jump. It just works OK.
    jumps = list(np.where(np.diff(new_lons[0, :], 2))[0][::2] + 2)

    # Beginning and sizes of continuous segments of lons.
    segs = [0] + jumps
    sizes = np.diff(segs + [len(new_lons[0, :])])

    # Sort according to value of lon at segment begin index.
    segs = zip(new_lons[0, segs], segs, sizes)
    src_segs = sorted(segs, key=lambda x : x[0])

    dest_idx = 0
    for i, (_, src_idx, size) in enumerate(src_segs):
        new_lons[:, dest_idx:dest_idx+size] = lons_copy[:, src_idx:src_idx+size]

        if new_data is not None:
            if len(data.shape) == 3:
                new_data[:, :, dest_idx:dest_idx+size] = data[:, :, src_idx:src_idx+size]
            else:
                new_data[:, dest_idx:dest_idx+size] = data[:, src_idx:src_idx+size]

        dest_idx += size

        if i+1 == len(segs):
            break

    return new_lons, new_data

def get_time_origin(filename):
    """
    Parse time.units to find the start/origin date of the file. Return a
    datetime.date object.
    """

    date_search_strings = ['\d{4}-\d{2}-\d{2}','\d{4}-\d{1}-\d{2}',
                            '\d{4}-\d{2}-\d{1}','\d{4}-\d{1}-\d{1}']

    with nc.Dataset(filename) as f:
        time_var = f.variables['time']
        assert 'months since' in time_var.units or \
               'days since' in time_var.units or \
               'hours since' in time_var.units, \
            "Time units doesn't have expected format: {}".format(time_var.units)
        for ds in date_search_strings:
            m = re.search(ds, time_var.units)
            if m is not None:
                break
        assert m is not None
        date = dt.datetime.strptime(m.group(0), '%Y-%m-%d')

    return dt.date(date.year, date.month, date.day)

def col_idx_largest_lat(lats):
    """
    The col index with the largest lat.
    """
    _, c  = np.unravel_index(np.argmax(lats), lats.shape)

    return c


def create_mom_output(ocean_grid, filename, start_date, history):

    f = nc.Dataset(filename, 'w')

    f.createDimension('lon', ocean_grid.num_lon_points)
    f.createDimension('lat', ocean_grid.num_lat_points)
    f.createDimension('depth', ocean_grid.num_levels)
    f.createDimension('time', 1)

    lons = f.createVariable('lon', 'f8', ('lon'))
    lons.long_name = 'Nominal Longitude of T-cell center'
    lons.units = 'degree_east'
    lons.modulo = 360.
    lons.point_spacing = 'even'
    lons.axis = 'X'
    # MOM needs this to be a single dimension
    lons[:] = ocean_grid.x_t[ocean_grid.x_t.shape[0] // 2, :]

    lats = f.createVariable('lat', 'f8', ('lat'))
    lats.long_name = 'Nominal Latitude of T-cell center'
    lats.units = 'degree_north'
    lats.point_spacing = 'uneven'
    lats.axis = 'Y'
    # MOM needs this to be a single dimension
    col = col_idx_largest_lat(ocean_grid.y_t[:])
    lats[:] = ocean_grid.y_t[:, col]

    depth = f.createVariable('depth', 'f8', ('depth'))
    depth.long_name = 'depth'
    depth.units = 'meters'
    depth.positive = 'down'
    depth.point_spacing = 'uneven'
    depth.axis = 'Z'
    depth[:] = ocean_grid.z[:]

    time = f.createVariable('time', 'f8', ('time'))
    time.long_name = 'time'
    time.units = "months since {}-{}-{} 00:00:00".format(str(start_date.year).zfill(4),
                                                       str(start_date.month).zfill(2),
                                                       str(start_date.day).zfill(2))
    time.cartesian_axis = "T"
    time.calendar_type = "noleap"
    time.calendar = "noleap"

    f.close()

def write_mom_output_at_time(filename, var_name, var_longname, var_units,
                             var_data, time_idx, time_pt, write_ic=False):

    with nc.Dataset(filename, 'r+') as f:
        if not var_name in f.variables:
            var = f.createVariable(var_name, 'f8',
                                   ('time', 'depth', 'lat', 'lon'),
                                   fill_value=-1.e+34, zlib=True, complevel=5, shuffle=True)
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
