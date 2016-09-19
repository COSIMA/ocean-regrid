#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
from grid import Grid
from util import normalise_lons

class TripolarGrid(Grid):
    """
    This is a global tripolar grid with given depth. To build it an existing
    tripolar grid is used.
    """

    def __init__(self, tripolar_grid, levels, description=''):

        # We may need to extend the input grid Southward.
        if np.min(tripolar_grid.y_t[:]) > -82:

            dy = tripolar_grid.y_t[1, 0] - tripolar_grid.y_t[0, 0]
            new_rows = int(np.rint((abs(-82 - np.min(tripolar_grid.y_t)) / dy)))

            x_t = np.zeros((tripolar_grid.x_t.shape[0] + new_rows, tripolar_grid.x_t.shape[1]))
            y_t = np.zeros((tripolar_grid.x_t.shape[0] + new_rows, tripolar_grid.x_t.shape[1]))

            x_t[new_rows:, :] = tripolar_grid.x_t[:]
            x_t[:new_rows, :] = np.stack([tripolar_grid.x_t[0, :]]*new_rows)

            y_t[new_rows:, :] = tripolar_grid.y_t[:]
            for n in range(new_rows-1, -1, -1):
                y_t[n, :] = y_t[n+1, :] - dy

            tripolar_mask = np.zeros_like(tripolar_grid.mask, dtype=bool)
            tripolar_mask[tripolar_grid.mask[:] == 0.0] = True
            assert len(tripolar_mask.shape) == 3

            # Drop the mask in, with new rows being masked by default.
            mask = np.ndarray((levels.shape[0], y_t.shape[0], y_t.shape[1]), dtype=bool)
            mask[:] = True
            diff = abs(mask.shape[0] - tripolar_grid.mask.shape[0])
            if mask.shape[0] > tripolar_mask.shape[0]:
                mask[:-diff, new_rows:, :] = tripolar_mask[:]
            else:
                mask[:, new_rows:, :] = tripolar_mask[:-diff, :, :]
        else:
            x_t = tripolar_grid.x_t[:]
            y_t = tripolar_grid.x_y[:]
            mask = tripolar_grid.mask[:]

        super(TripolarGrid, self).__init__(x_t, y_t, levels, mask, description)

    def make_corners(self):

        x = self.x_t
        y = self.y_t

        dx_half = np.empty_like(x)
        dy_half = np.empty_like(x)

        dx_half[:,1:] = (x[:, 1:] - x[:, 0:-1]) / 2.0
        dy_half[1:,:] = (y[1:, :] - y[0:-1, :]) / 2.0

        # Need to extend South
        dy_half[0, 1:] = dy_half[1, 1:]
        dx_half[0, 1:] = dx_half[1, 1:]

        # and West
        dy_half[:, 0] = dy_half[:, 1]
        dx_half[:, 0] = dx_half[:, 1]

        clon = np.empty((self.num_lat_points, self.num_lon_points, 4))
        clon[:] = np.NAN

        clon[:,:,0] = x - dx_half
        clon[:,:,1] = x + dx_half
        clon[:,:,2] = x + dx_half
        clon[:,:,3] = x - dx_half
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.num_lat_points, self.num_lon_points, 4))
        clat[:] = np.NAN
        clat[:,:,0] = y - dy_half
        clat[:,:,1] = y - dy_half
        clat[:,:,2] = y + dy_half
        clat[:,:,3] = y + dy_half
        assert(not np.isnan(np.sum(clat)))

        self.clon_t = clon
        self.clat_t = clat
