#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc

from .grid import Grid

class OrasGrid(Grid):

    def __init__(self, grid_defs, description=''):

        # FIXME: only support T grid for now.
        grid_def = None
        for gd in grid_defs:
            if 'coordinates_grid_T.nc' in gd:
                grid_def = gd
        assert grid_def is not None

        with nc.Dataset(grid_def) as f:

            # The oras grid is unusual. It has one duplicate row and two
            # duplicate columns. We cut these off. It also has a fully masked
            # level at the bottom, also cut this off.
            try:
                x_t = f.variables['nav_lon'][:-1, :-2]
                y_t = f.variables['nav_lat'][:-1, :-2]
                z = f.variables['deptht'][:-1]
            except KeyError:
                x_t = f.variables['lon'][:-1, :-2]
                y_t = f.variables['lat'][:-1, :-2]
                z = f.variables['depth'][:-1]

            mask = np.zeros_like(f.variables['tmask'][:-1, :-1, :-2], dtype=bool)
            mask[f.variables['tmask'][:-1, :-1, :-2] == 0.0] = True

        super(OrasGrid, self).__init__(x_t, y_t, z, mask, description)


    def fix_data_shape(self, data_array):
        """
        Since the ORAS grid has had duplicate rows and columns cut off we need
        to do the same for the data. FIXME: better way to do this.
        """

        assert len(data_array.shape) >= 2 and len(data_array.shape) <= 4

        if len(data_array.shape) == 4:
            new_array = data_array[:, :-1, :-1, :-2]
            assert(new_array.shape[2] == self.x_t.shape[0])
            assert(new_array.shape[3] == self.x_t.shape[1])
        elif len(data_array.shape) == 3:
            new_array = data_array[:-1, :-1, :-2]
            assert(new_array.shape[1] == self.x_t.shape[0])
            assert(new_array.shape[2] == self.x_t.shape[1])
        else:
            new_array = data_array[:-1, :-2]
            assert(new_array.shape[0] == self.x_t.shape[0])
            assert(new_array.shape[1] == self.x_t.shape[1])

        return new_array

    def apply_grid_mask(self, data_array):

        assert len(data_array.shape) >= 2 and len(data_array.shape) <= 4

        if len(data_array.shape) == 4:
            new_mask = np.stack([self.mask[:]] * data_array.shape[0], axis=0)
            return np.ma.array(data_array, mask=new_mask)
        elif len(data_array.shape) == 3:
            return np.ma.array(data_array, mask=self.mask)
        else:
            return np.ma.array(data_array, mask=self.mask[0, :, :])

    def make_corners(self):
        raise exceptions.NotImplementedError
