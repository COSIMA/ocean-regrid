#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
from grid import Grid

class NemoGrid(Grid):

    def __init__(self, h_grid_def, v_grid_def, mask_file, description):

        with nc.Dataset(h_grid_def) as f:

            # Get t-points.
            x_t = f.variables['glamt'][:]
            y_t = f.variables['gphit'][:]

            # These variables hold the corners
            self.glamf = f.variables['glamf'][:]
            self.gphif = f.variables['gphif'][:]

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

    def make_corners(self):

        # These are the top righ-hand corner of t cells.
        glamf = self.glamf
        gphif = self.gphif

        # Extend south so that Southern most cells can have bottom corners.
        gphif_new = np.ndarray((gphif.shape[0] + 1, gphif.shape[1] + 1))
        gphif_new[1:, 1:] = gphif[:]
        gphif_new[0, 1:] = gphif[0, :] - abs(gphif[1, :] - gphif[0, :])

        glamf_new = np.ndarray((glamf.shape[0] + 1, glamf.shape[1] + 1))
        glamf_new[1:, 1:] = glamf[:]
        glamf_new[0, 1:] = glamf[0, :]

        # Repeat first longitude so that Western most cells have left corners.
        gphif_new[:, 0] = gphif_new[:, -1]
        glamf_new[:, 0] = glamf_new[:, -1]

        gphif = gphif_new
        glamf = glamf_new

        # Corners of t points. Index 0 is bottom left and then
        # anti-clockwise.
        clon = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clon[:] = np.NAN
        clon[:,:,0] = glamf[0:-1,0:-1]
        clon[:,:,1] = glamf[0:-1,1:]
        clon[:,:,2] = glamf[1:,1:]
        clon[:,:,3] = glamf[1:,0:-1]
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.x_t.shape[0], self.x_t.shape[1], 4))
        clat[:] = np.NAN
        clat[:,:,0] = gphif[0:-1,0:-1]
        clat[:,:,1] = gphif[0:-1,1:]
        clat[:,:,2] = gphif[1:,1:]
        clat[:,:,3] = gphif[1:,0:-1]
        assert(not np.isnan(np.sum(clat)))

        self.clon_t = clon
        self.clat_t = clat
