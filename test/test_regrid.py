
from __future__ import print_function

import pytest
import os
import subprocess as sp
import sh
import netCDF4 as nc
import numpy as np

data_tarball = 'test_data.tar.gz'
data_tarball_url = 'http://s3-ap-southeast-2.amazonaws.com/dp-drop/ocean-regrid/test/test_data.tar.gz'

def check_output_fields(model_name, output):

    with nc.Dataset(output) as f:
        if 'MOM' in model_name:

            if f.variables.has_key('temp'):
                assert f.variables['temp'].units == 'C' or \
                    'celsius' in f.variables['temp'].units.lower()
                var = f.variables['temp']
                var_name = 'temp'
            else:
                var = f.variables['salt']
                var_name = 'salt'

            x_t = f.variables['GRID_X_T'][:]
            y_t = f.variables['GRID_Y_T'][:]

            # Make sure that lats and lons are nicely spaced apart. Don't want
            # any big jumps.  These are 'representational values' for lat and
            # lon only because they are a single dimension on a irregular grid.
            assert np.max(np.diff(x_t)) <= (360. / float(x_t.shape[0])) * 2
            assert np.max(np.diff(y_t)) <= (180. / float(y_t.shape[0])) * 2

        else:
            assert model_name == 'NEMO'
            var = f.variables['votemper']
            var_name = 'temp'

        assert var_name == 'temp' or var_name == 'salt'

        if var_name == 'temp':
            if var.units == 'C' or 'celsius' in var.units.lower():
                for t in range(var.shape[0]):
                    assert np.max(var[t, :]) < 40.0
                    assert np.min(var[t, :]) > -10.0
            else:
                for t in range(var.shape[0]):
                    assert np.max(var[t, :]) < 320.0
                    assert np.min(var[t, :]) > 260.0
        else:
            for t in range(var.shape[0]):
                assert np.max(var[t, :]) < 50.0
                assert np.min(var[t, :]) >= 0.0


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


    @pytest.mark.mom
    def test_oras4_to_mom(self, input_dir, output_dir):
        """
        Regrid ORAS4 to MOM.
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
        check_output_fields('MOM', output)
        assert(os.path.exists(output))

    @pytest.mark.mom1
    def test_oras4_to_mom1(self, input_dir, output_dir):
        """
        Regrid ORAS4 to MOM 1 degree.
        """

        output = os.path.join(output_dir, 'mom1_oras4_temp.nc')
        if os.path.exists(output):
            os.remove(output)

        src_name = 'ORAS4'
        src_data_file = os.path.join(input_dir, 'so_oras4_1m_2014_grid_T.nc')
        dest_name = 'MOM1'
        dest_data_file = output

        args = [src_name, src_data_file, 'salt', dest_name, dest_data_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'regrid_simple.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Check that outputs exist.
        check_output_fields('MOM', output)
        assert(os.path.exists(output))

    @pytest.mark.godas
    @pytest.mark.nemo
    def test_godas_to_nemo(self, input_dir, output_dir):
        """
        Regrid ORAS4 to MOM.
        """

        output = os.path.join(output_dir, 'nemo_godas_temp.nc')
        if os.path.exists(output):
            os.remove(output)

        src_name = 'GODAS'
        src_data_file = os.path.join(input_dir, 'pottmp.2004.nc')
        dest_name = 'NEMO'
        dest_data_file = output

        args = [src_name, src_data_file, 'temp', dest_name, dest_data_file]

        my_dir = os.path.dirname(os.path.realpath(__file__))
        cmd = [os.path.join(my_dir, '../', 'regrid_simple.py')] + args
        ret = sp.call(cmd)
        assert(ret == 0)

        # Check that outputs exist.
        assert os.path.exists(output)
        check_output_fields('NEMO', output)
