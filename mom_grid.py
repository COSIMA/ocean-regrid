#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc

from grid import Grid


class MomGrid(Grid):

    def __init__(self, h_grid_def, v_grid_def, mask_file, description):

        with nc.Dataset(h_grid_def) as f:

            # Select points from double density horizontal grid. Only
            # need t-points.
            x_t = f.variables['x'][1::2,1::2]
            y_t = f.variables['y'][1::2,1::2]
            self.x_vt = f.variables['x'][:]
            self.y_vt = f.variables['y'][:]

        with nc.Dataset(v_grid_def) as f:
            # Only take cell centres.
            z = f.variables['zeta'][1::2]

        with nc.Dataset(mask_file) as f:
            mask = np.zeros_like(f.variables['mask'], dtype=bool)
            mask[f.variables['mask'][:] == 0.0] = True

        super(MomGrid, self).__init__(x_t, y_t, z, mask, description)

        self.num_lat_points = y_t.shape[0]
        self.num_lon_points = y_t.shape[1]

    def make_corners(self):

        # Uses double density grid to figure out corners.
        x = self.x_vt
        y = self.y_vt

        # Corners of t points. Index 0 is bottom left and then
        # anti-clockwise.
        clon = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clon[:] = np.NAN
        clon[:,:,0] = x[0:-1:2,0:-1:2]
        clon[:,:,1] = x[0:-1:2,2::2]
        clon[:,:,2] = x[2::2,2::2]
        clon[:,:,3] = x[2::2,0:-1:2]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clat[:] = np.NAN
        clat[:,:,0] = y[0:-1:2,0:-1:2]
        clat[:,:,1] = y[0:-1:2,2::2]
        clat[:,:,2] = y[2::2,2::2]
        clat[:,:,3] = y[2::2,0:-1:2]
        assert(not np.isnan(np.sum(clat)))

        self.clon_t = clon
        self.clat_t = clat
