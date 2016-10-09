
from __future__ import print_function

import pytest
import os
import subprocess as sp
import sh

data_tarball = 'test_data.tar.gz'
data_tarball_url = 'http://s3-ap-southeast-2.amazonaws.com/dp-drop/ocean-regrid/test/test_data.tar.gz'

class TestRegrid():

    @pytest.fixture
    def input_dir(self):
        test_dir = os.path.dirname(os.path.realpath(__file__))
        test_data_dir = os.path.join(test_dir, 'test_data')
        test_data_tarball = os.path.join(test_dir, data_tarball)

        if not os.path.exists(test_data_dir):
            if not os.path.exists(test_data_tarball):
                sh.wget('-P', test_dir, data_tarball_url)
            sh.tar('zxvf', test_data_tarball, '-C', test_dir)

        return os.path.join(test_data_dir, 'input')

    @pytest.fixture
    def output_dir(self):
        test_dir = os.path.dirname(os.path.realpath(__file__))
        test_data_dir = os.path.join(test_dir, 'test_data')

        return os.path.join(test_data_dir, 'output')

    def test_mom(self, input_dir, output_dir):
        """
        Test script makes oasis grids from mom inputs.
        """

        output = os.path.join(output_dir, 'mom_oras4_temp.nc')
        if os.path.exists(output):
            os.remove(output)

        src_name = 'ORAS4'
        src_data_file = os.path.join(input_dir, 'thetao_oras4_1m_2014_grid_T.nc')
        dest_name = 'MOM'
        dest_data_file = output

        args = [src_name, src_data_file, 'temp', dest_name, dest_data_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'regrid_simple.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Check that outputs exist.
        assert(os.path.exists(output))
