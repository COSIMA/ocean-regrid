#!/usr/bin/env python

from __future__ import print_function

import sys, os
import argparse
import netCDF4 as nc

import regrid

"""
Regrid Ocean reanalysis
"""

grid_defs_error = \
"""
Grid definitions directory {} not found.
please download it with:
wget http://s3-ap-southeast-2.amazonaws.com/dp-drop/ocean-nudge/grid_defs.tar.gz
and unzip into the same directory as this executable.
"""

def grid_defs_dir():
    """
    Get path to directory where MOM, NEMO, GODAS and ORAS4 grid definitions are
    found.
    """

    if getattr(sys, 'frozen', False):
        basedir = sys._MEIPASS
    else:
        basedir = os.path.dirname(os.path.realpath(__file__))

    return os.path.join(basedir, 'grid_defs')

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('reanalysis_name', help="""
                        Name of src data/grid, must be GODAS or ORAS4""")
    parser.add_argument('reanalysis_data_file',
                        help='Reanalysis data file to regrid.')
    parser.add_argument('reanalysis_var', help="""
        Reanalysis variable to regrid. Must be 'salt' or 'temp'""")
    parser.add_argument('model_name', help="""
                        Name of model, must be MOM or NEMO""")
    parser.add_argument('output_file', help='Name of the destination/output file.')
    parser.add_argument('--regrid_weights', help='Name file to store regridding weights.')
    args = parser.parse_args()

    assert args.model_name == 'MOM' or args.model_name == 'NEMO'
    assert args.reanalysis_name == 'GODAS' or args.reanalysis_name == 'ORAS4'
    assert args.reanalysis_var == 'salt' or args.reanalysis_var == 'temp'

    # Set up the reanalysis model and grid definitions.
    grid_defs = grid_defs_dir()
    if not os.path.exists(grid_defs):
        print(grid_defs_error.format(grid_defs), file=sys.stderr)
        return 1

    if args.reanalysis_name == 'GODAS':
        reanalysis_hgrids = (os.path.join(grid_defs, 'pottmp.2016.nc'),)
        reanalysis_vgrid = os.path.join(grid_defs, 'pottmp.2016.nc')
    else:
        reanalysis_hgrids = (os.path.join(grid_defs, 'coordinates_grid_T.nc'),
                             os.path.join(grid_defs, 'coordinates_grid_U.nc'),
                             os.path.join(grid_defs, 'coordinates_grid_V.nc'))
        reanalysis_vgrid = os.path.join(grid_defs, 'coordinates_grid_T.nc')

    if args.model_name == 'MOM':
        model_hgrid = os.path.join(grid_defs, 'ocean_hgrid.nc')
        model_vgrid = os.path.join(grid_defs, 'ocean_vgrid.nc')
        model_mask = os.path.join(grid_defs, 'ocean_mask.nc')
    else:
        model_hgrid = os.path.join(grid_defs, 'coordinates.nc')
        model_vgrid = os.path.join(grid_defs, 'data_1m_potential_temperature_nomask.nc')
        model_mask = None

    # Read in temperature and salinity data.
    if args.reanalysis_name == 'ORAS4':
        temp_src_var = 'thetao'
        salt_src_var = 'so'
    else:
        temp_src_var = 'pottmp'
        salt_src_var = 'salt'

    if args.model_name == 'MOM':
        temp_dest_var = 'temp'
        salt_dest_var = 'salt'
    else:
        temp_dest_var = 'votemper'
        salt_dest_var = 'vosaline'

    if args.reanalysis_var == 'salt':
        src_var = salt_src_var
        dest_var = salt_dest_var
    else:
        src_var = temp_src_var
        dest_var = temp_dest_var

    # Regrid temp and salt, write out to the same file.
    weights = None
    weights = regrid.do_regridding(args.reanalysis_name, reanalysis_hgrids,
                                   reanalysis_vgrid,
                                   args.reanalysis_data_file, src_var,
                                   args.model_name, model_hgrid,
                                   model_vgrid, args.output_file, dest_var,
                                   dest_mask=model_mask, regrid_weights=weights)
    if weights is None:
        return 1

if __name__ == '__main__':
    sys.exit(main())
