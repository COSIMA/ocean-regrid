#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
import exceptions

from grid import Grid

def find_nearest_index(array, value):
    return (np.abs(array - value)).argmin()

class OrasGrid(Grid):

    def __init__(self, grid_def, description=''):

        with nc.Dataset(grid_def) as f:

            # Select points from double density horizontal grid. Only
            # need t-points.
            try:
                x_t = f.variables['nav_lon'][:]
                y_t = f.variables['nav_lat'][:]
                z = f.variables['deptht'][:]
            except KeyError:
                x_t = f.variables['lon'][:]
                y_t = f.variables['lat'][:]
                z = f.variables['depth'][:]

            mask = np.zeros_like(f.variables['tmask'], dtype=bool)
            mask[f.variables['tmask'][:] == 0.0] = True

        super(OrasGrid, self).__init__(x_t, y_t, z, mask, description)

    def make_corners(self):
        raise exceptions.NotImplementedError
