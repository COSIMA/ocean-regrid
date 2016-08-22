
import netCDF4 as nc
from mom_grid import MomGrid

def write_mom_ic(ocean_grid, temp_data, salt_data, filename, history):

    f = nc.Dataset(filename, 'w')

    f.createDimension('GRID_Y_T', ocean_grid.num_lat_points)
    f.createDimension('GRID_X_T', ocean_grid.num_lon_points)
    f.createDimension('ZT', ocean_grid.num_levels)

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

    temp = f.createVariable('temp', 'f8', ('ZT', 'GRID_Y_T', 'GRID_X_T'), fill_value=-1.e+34)
    temp.missing_value = -1.e+34
    temp.long_name = "Temperature"
    temp.units = "deg C"
    temp.history = history
    temp[:] = temp_data[:]

    salt = f.createVariable('salt', 'f8', ('ZT', 'GRID_Y_T', 'GRID_X_T'), fill_value=-1.e+34)
    salt.missing_value = -1.e+34
    salt.long_name = "Salinity"
    salt.units = "psu"
    salt.history = history
    salt[:] = salt_data[:]

    f.close()

def write_nemo_ic(ocean_grid, temp, salt, filename, history):
    pass
