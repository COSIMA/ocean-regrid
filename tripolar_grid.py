#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc
import exceptions
from grid import Grid

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
            for n in range(new_rows-1, 0, -1):
                y_t[n, :] = y_t[n+1, :] - dy

            # Drop the mask in, with new rows being masked by default.
            mask = np.ndarray((levels.shape[0], y_t.shape[0], y_t.shape[1]))
            mask[:] = True
            mask[levels.shape[0]-tripolar_grid.mask.shape[0]:, new_rows:, :] = tripolar_grid.mask[:]
        else:
            x_t = tripolar_grid.x_t[:]
            y_t = tripolar_grid.x_y[:]
            mask = tripolar_grid.mask[:]

        super(TripolarGrid, self).__init__(x_t, y_t, levels, mask, description)

    def make_corners(self):
        raise exceptions.NotImplementedError
