#!/usr/bin/env python

from __future__ import print_function

import sys, os
import subprocess as sp
import argparse

"""
Regrid Ocean reanalysis
"""

grid_defs_error = \
"""
Grid definitions directory {} not found.
please download it with:
wget http://s3-ap-southeast-2.amazonaws.com/dp-drop/ocean-regrid/grid_defs.tar.gz
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
    parser.add_argument('src_name', help="""
                        Name of src data/grid, must be GODAS or ORAS4""")
    parser.add_argument('src_data_file',
                        help='Source data file to regrid.')
    parser.add_argument('src_var', help="""
        Source variable to regrid. Must be 'salt' or 'temp'""")
    parser.add_argument('dest_name', help="""
                        Name of desintation grid, must be MOM or NEMO""")
    parser.add_argument('output_file', help='Name of the destination/output file.')
    args = parser.parse_args()

    assert args.dest_name == 'MOM' or args.dest_name == 'MOM1' or \
        args.dest_name == 'NEMO'
    assert args.src_name == 'GODAS' or args.src_name == 'ORAS4'
    assert args.src_var == 'salt' or args.src_var == 'temp'

    # Set up the src and dest and grid definitions.
    grid_defs = grid_defs_dir()
    if not os.path.exists(grid_defs):
        print(grid_defs_error.format(grid_defs), file=sys.stderr)
        return 1

    if args.src_name == 'GODAS':
        src_hgrids = (os.path.join(grid_defs, 'pottmp.2016.nc'),)
        src_vgrid = os.path.join(grid_defs, 'pottmp.2016.nc')
    else:
        src_hgrids = (os.path.join(grid_defs, 'coordinates_grid_T.nc'),
                      os.path.join(grid_defs, 'coordinates_grid_U.nc'),
                      os.path.join(grid_defs, 'coordinates_grid_V.nc'))
        src_vgrid = os.path.join(grid_defs, 'coordinates_grid_T.nc')

    if args.dest_name == 'MOM':
        dest_hgrid = os.path.join(grid_defs, 'ocean_hgrid.nc')
        dest_vgrid = os.path.join(grid_defs, 'ocean_vgrid.nc')
        dest_mask = os.path.join(grid_defs, 'ocean_mask.nc')
        mm_arg = ['--dest_mask', dest_mask]
    elif args.dest_name == 'MOM1':
        dest_hgrid = os.path.join(grid_defs, 'grid_spec.nc')
        dest_vgrid = os.path.join(grid_defs, 'grid_spec.nc')
        dest_mask = os.path.join(grid_defs, 'grid_spec.nc')
        mm_arg = ['--dest_mask', dest_mask]
    else:
        dest_hgrid = os.path.join(grid_defs, 'coordinates.nc')
        dest_vgrid = os.path.join(grid_defs, 'data_1m_potential_temperature_nomask.nc')
        dest_mask = None
        mm_arg = []

    # Read in temperature and salinity data.
    if args.src_name == 'ORAS4':
        temp_src_var = 'thetao'
        salt_src_var = 'so'
    else:
        temp_src_var = 'pottmp'
        salt_src_var = 'salt'

    if args.dest_name == 'MOM' or args.dest_name == 'MOM1':
        temp_dest_var = 'temp'
        salt_dest_var = 'salt'
    else:
        temp_dest_var = 'votemper'
        salt_dest_var = 'vosaline'

    if args.src_var == 'salt':
        src_var = salt_src_var
        dest_var = salt_dest_var
    else:
        src_var = temp_src_var
        dest_var = temp_dest_var

    # Regrid temp and salt, write out to the same file.
    args = [args.src_name, src_hgrids[0], src_vgrid,
            args.src_data_file, src_var, args.dest_name,
            dest_hgrid, dest_vgrid, args.output_file, dest_var] + mm_arg
    exe = os.path.join(os.path.dirname(__file__), 'regrid.py')
    return sp.call([exe] + args)

if __name__ == '__main__':
    sys.exit(main())
