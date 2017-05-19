#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
from .grid import Grid

class GodasGrid(Grid):

    def __init__(self, grid_def, description=''):

        with nc.Dataset(grid_def) as f:

            # Select points from double density horizontal grid. Only
            # need t-points.
            x_t = f.variables['lon'][:]
            y_t = f.variables['lat'][:]
            try:
                z = f.variables['level'][:]
            except KeyError:
                z = f.variables['depth'][:]

            try:
                mask = f.variables['pottmp'][0, :].mask[:]
            except KeyError:
                mask = f.variables['POT'][0, :].mask[:]

        super(GodasGrid, self).__init__(x_t, y_t, z, mask, description)


    def set_mask(self, new_mask):
        self.mask = new_mask[:]
