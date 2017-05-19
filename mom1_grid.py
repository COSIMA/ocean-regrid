#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc

from .grid import Grid

class Mom1Grid(Grid):

    def __init__(self, h_grid_def, v_grid_def, mask_file,
                    description='MOM 1 deg tripolar'):
        """
        MOM 1 degree grid.
        """

        with nc.Dataset(h_grid_def) as f:

            # Select points from double density horizontal grid.
            # t-cells.
            x_t = f.variables['x_T'][:]
            y_t = f.variables['y_T'][:]

            # u-cells.
            x_u = f.variables['x_C'][:]
            y_u = f.variables['y_C'][:]

            self.area_t = f.variables['area_T'][:]
            self.area_u = f.variables['area_C'][:]

            self.clon_t = f.variables['x_vert_T'][:]
            self.clat_t = f.variables['y_vert_T'][:]
            self.clon_u = f.variables['x_vert_C'][:]
            self.clat_u = f.variables['y_vert_C'][:]

        with nc.Dataset(v_grid_def) as f:
            z = f.variables['zt'][:]

        with nc.Dataset(mask_file) as f:
            mask = np.zeros_like(f.variables['wet'], dtype=bool)
            mask[f.variables['wet'][:] == 0.0] = True

        super(Mom1Grid, self).__init__(x_t, y_t, z, mask, description)

    def make_corners(self):
        """
        Corners have already been read.
        """
        pass
