#!/usr/bin/env python

from __future__ import print_function

import numpy as np
from grid import Grid

class RegularGrid(Grid):

    def __init__(self, num_lons, num_lats, levels, mask=None, description=''):

        dx = 360.0 / num_lats
        dy = 180.0 / num_lons
        dx_half = dx / 2
        dy_half = dy / 2

        # Set lats and lons.
        lons = np.linspace(0, 360, num_lons, endpoint=False)
        # lat points exclude the poles.
        lats = np.linspace(-90 + dy_half, 90 - dy_half, num_lats)

        super(RegularGrid, self).__init__(lons, lats, levels, mask, description)
