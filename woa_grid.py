#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
from .grid import Grid

class WoaGrid(Grid):

    def __init__(self, grid_def, description=''):

        with nc.Dataset(grid_def) as f:

            # Get lon and lat
            x_t = f.variables['lon'][:]
            y_t = f.variables['lat'][:]
            try:
                z = f.variables['depth'][:]
            except KeyError:
                z = f.variables['level'][:]

            try:
                mask = f.variables['practical_salinity'][0, :, :, :].mask[:]
            except KeyError:
                mask = f.variables['s_an'][0, :, :, :].mask[:]

        super(WoaGrid, self).__init__(x_t, y_t, z, mask, description)


    def set_mask(self, new_mask):
        self.mask = new_mask[:]
