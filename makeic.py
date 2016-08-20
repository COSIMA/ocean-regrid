#!/usr/bin/env python

from __future__ import print_function

import sys, os
import argparse
from scipy.interpolate import UnivariateSpline
from scipy import interp
import numpy as np
import netCDF4 as nc

from mom_grid import MomGrid
from grid import Grid

import ESMF

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

    args = parser.parse_args()

    title = 'MOM tripolar t-cell grid'
    mom_grid = MomGrid(args.ocean_hgrid, args.ocean_vgrid, args.ocean_mask, title)
    mom_scrip_file = 'model_scrip.nc'
    mom_grid.write_scrip(mom_scrip_file, ' '.join(sys.argv))
    # Dest grid in scrip format.
    dest_grid = ESMF.Grid(filename=mom_scrip_file,
                            filetype=ESMF.FileFormat.SCRIP)

    # Read in obs file.
    with nc.Dataset(args.temp_obs_file) as obs:
        temp = obs.variables[args.temp_var][args.time_index, :]
        xt = obs.variables[args.x_var][:]
        yt = obs.variables[args.y_var][:]
        z = obs.variables[args.z_var][:]

    # Source grid in scrip format.
    title = '{}x{} Cylindrical Equidistant Projection Grid'.format(len(xt), len(yt))
    obs_grid = Grid(xt, yt, z, temp.mask[0, :, :], title)
    obs_scrip_file = 'obs_scrip.nc'
    obs_grid.write_scrip(obs_scrip_file, ' '.join(sys.argv))

    src_grid = ESMF.Grid(filename=obs_scrip_file,
                            filetype=ESMF.FileFormat.SCRIP)

    dest_grid = ESMF.Grid(filename=mom_scrip_file,
                            filetype=ESMF.FileFormat.SCRIP)

    #dest_grid = ESMF.Grid(filename=mom_scrip_file,
    #                        filetype=ESMF.FileFormat.SCRIP)

    # Create temp and salinity source fields.
    temp_src = ESMF.Field(src_grid, 'temp_src',
                            staggerloc=ESMF.StaggerLoc.CENTER)
    temp_dest = ESMF.Field(dest_grid, 'temp_dest',
                            staggerloc=ESMF.StaggerLoc.CENTER)

    temp_dest.data[:] = 1e20
    print('HERE 1')
    import pdb
    pdb.set_trace()

    # Create an object to regrid data 
    regrid = ESMF.Regrid(temp_src, temp_dest,
			 #regrid_method=ESMF.RegridMethod.BILINEAR,
			 regrid_method=ESMF.RegridMethod.NEAREST_STOD,
			 unmapped_action=ESMF.UnmappedAction.ERROR)
    print('HERE 2')

    # Regrid obs columns onto model vertical grid.
    temp = regrid_columns(temp, z, mom_grid.z, plot_results=True)
    temp_src.data[:] = temp[0, :, :].transpose()

    print('HERE 3')
    import pdb
    pdb.set_trace()

    # Do the regridding
    temp_dest = regrid(temp_src, temp_dest)

    print('HERE 4')

    import pdb
    pdb.set_trace()

if __name__ == '__main__':
    sys.exit(main())
