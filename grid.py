#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import netCDF4 as nc

class Grid(object):

    def __init__(self, lons, lats, levels, mask, description):

        self.num_lat_points = len(lats)
        self.num_lon_points = len(lons)
        self.num_levels = len(levels)

        if len(lons.shape) == 1:
            # We expect this to be a regular grid.
            assert np.allclose(np.diff(lons),
                     np.array([lons[1] - lons[0]]*(len(lons)-1)))
            assert np.allclose(np.diff(lats),
                     np.array([lats[1] - lats[0]]*(len(lats)-1)), atol=1e-4)
            # Turn into tiled
            self.x_t = np.tile(lons, (self.num_lat_points, 1))
            self.y_t = np.tile(lats, (self.num_lon_points, 1))
            self.y_t = self.y_t.transpose()
        else:
            self.x_t = lons
            self.y_t = lats

        self.z = levels
        self.description = description
        self.mask = mask

        self.clon_t = None
        self.clat_t = None

    def make_corners(self):

	x = self.x_t
	y = self.y_t

        dx = (x[0][1] - x[0][0]) / 2.0
        dy = (y[1][0] - y[0][0]) / 2.0

        # Set grid corners, we do these one corner at a time. Start at the 
        # bottom left and go anti-clockwise. This is the SCRIP convention.
        clon = np.empty((self.num_lat_points, self.num_lon_points, 4))
        clon[:] = np.NAN
        clon[:,:,0] = x - dx
        clon[:,:,1] = x + dx
        clon[:,:,2] = x + dx
        clon[:,:,3] = x - dx
        assert(not np.isnan(np.sum(clon)))

        clat = np.empty((self.num_lat_points, self.num_lon_points, 4))
        clat[:] = np.NAN
        clat[:,:,0] = y - dy
        clat[:,:,1] = y - dy
        clat[:,:,2] = y + dy
        clat[:,:,3] = y + dy
        assert(not np.isnan(np.sum(clat)))

        # The bottom latitude band should always be Southern extent.
        assert(np.all(clat[0, :, 0] == np.min(y) - dy))
        assert(np.all(clat[0, :, 1] == np.min(y) - dy))

        # The top latitude band should always be Northern extent.
        assert(np.all(clat[-1, :, 2] == np.max(y) + dy))
        assert(np.all(clat[-1, :, 3] == np.max(y) + dy))

        self.clon_t = clon
        self.clat_t = clat

    def write_scrip(self, filename, history):

        self.make_corners()

        f = nc.Dataset(filename, 'w')

        x = self.x_t
        y = self.y_t

        clat = self.clat_t
        clon = self.clon_t
        num_points = self.num_lat_points * self.num_lon_points

        f.createDimension('grid_size', num_points)
        f.createDimension('grid_corners', 4)
        f.createDimension('grid_rank', 2)

        grid_dims = f.createVariable('grid_dims', 'i4', ('grid_rank'))
        # SCRIP likes lon, lat
        grid_dims[:] = [self.num_lon_points, self.num_lat_points]

        center_lat = f.createVariable('grid_center_lat', 'f8', ('grid_size'))
        center_lat.units = 'degrees'
        center_lat[:] = y[:].flatten()

        center_lon = f.createVariable('grid_center_lon', 'f8', ('grid_size'))
        center_lon.units = 'degrees'
        center_lon[:] = x[:].flatten()

        imask = f.createVariable('grid_imask', 'i4', ('grid_size'))
        imask.units = 'unitless'
        # Invert the mask. SCRIP uses zero for points that do not
        # participate.
        imask[:] = np.invert(self.mask[:]).flatten()

        corner_lat = f.createVariable('grid_corner_lat', 'f8',
                                      ('grid_size', 'grid_corners'))
        corner_lat.units = 'degrees'
        corner_lat[:] = clat[:].flatten()

        corner_lon = f.createVariable('grid_corner_lon', 'f8',
                                      ('grid_size', 'grid_corners'))
        corner_lon.units = 'degrees'
        corner_lon[:] = clon[:].flatten()

        f.title = self.description
        f.history = history
        f.close()

