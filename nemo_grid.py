#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
from grid import Grid

class NemoGrid(Grid):

    def __init__(self, h_grid_def, v_grid_def, mask_file, description):

        with nc.Dataset(h_grid_def) as f:

            # Select points from double density horizontal grid. Only
            # need t-points.
            x_t = f.variables['nav_lon'][:]
            y_t = f.variables['nav_lat'][:]

        with nc.Dataset(v_grid_def) as f:
            z = f.variables['depth'][:]

        if mask_file is None:
            mask = np.zeros_like(x_t, dtype=bool)
        else:
            with nc.Dataset(mask_file) as f:
                mask = np.zeros_like(f.variables['mask'], dtype=bool)
                mask[f.variables['mask'][:] == 0.0] = True

        super(NemoGrid, self).__init__(x_t, y_t, z, mask, description)

        self.num_lat_points = y_t.shape[0]
        self.num_lon_points = y_t.shape[1]
