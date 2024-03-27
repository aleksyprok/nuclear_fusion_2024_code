"""
This module contains unit tests for the prepare_profiles module.

The prepare_profiles module provides functions for preparing profiles for a given SPR string.
"""

import os
import netCDF4
from python_scripts import prepare_profiles

def test_get_cdf_filename():
    """
    Test case for the get_cdf_filename function.

    This test case verifies that the get_cdf_filename function returns the expected filename
    for a given SPR string.
    """
    spr_string = 'SPR-045-14'
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(repo_path, "input_data")
    expected_filename = os.path.join(
        repo_path,
        'input_data',
        f'profiles_{spr_string}.CDF'
    )
    acutal_filename = prepare_profiles.get_cdf_filename(spr_string, input_dir)
    assert os.path.abspath(acutal_filename) == expected_filename

def test_read_cdf_file():
    """
    Test case for the read_cdf_file function.

    This test case verifies that the read_cdf_file function reads the CDF file correctly
    and returns the expected data.
    """
    spr_string = 'SPR-045-14'
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(repo_path, "input_data")
    cdf_filename = os.path.join(
        input_dir,
        f'profiles_{spr_string}.CDF'
    )
    psin, ti, te, ne, nd, nt = prepare_profiles.read_cdf_file(cdf_filename)

    # Check that psin[0] is approximately 0 i.e. < 1e-4
    assert abs(psin[0]) < 1e-4
    # Check that ti, te and ne are greater in the core than the edge
    assert ti[0] > ti[-1]
    assert te[0] > te[-1]
    assert ne[0] > ne[-1]
    assert nd[0] > nd[-1]
    assert nt[0] > nt[-1]
    # Check ne > nd > nt
    assert ne[0] > nd[0] > nt[0]

def test_get_gfile():
    """
    Test case for the get_gfile function.

    This test case verifies that the get_gfile function returns the expected gfile path
    for a given gfile filename.
    """
    gfile_filename = 'SPR-045-14.eqdsk'
    current_dir = os.path.dirname(__file__)
    input_data_dir = os.path.join(current_dir, '..', 'input_data')
    gfile_path = os.path.join(input_data_dir, gfile_filename)
    gfile = prepare_profiles.get_gfile(gfile_path)
    assert abs(gfile.rmaxis - 4.211483193) < 1e-6
    assert abs(gfile.zmaxis - 0.0) < 1e-6

def test_calculate_number_of_impurities():
    """
    Test case for the calculate_number_of_impurities function.

    This test case verifies that the calculate_number_of_impurities function returns the expected
    number of impurities for a given CDF file.
    """
    spr_string = 'SPR-045-14'
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(repo_path, "input_data")
    cdf_filename = os.path.join(
        input_dir,
        f'profiles_{spr_string}.CDF'
    )
    num_of_impurities = prepare_profiles.calculate_number_of_impurities(cdf_filename)
    assert num_of_impurities == 3

def test_main():
    """
    Test case for the main function.
    """
    num_markers = 10
    spr_string = 'SPR-045-14'
    current_dir = os.path.dirname(__file__)
    input_data_dir = os.path.join(current_dir, '..', 'input_data')
    output_dir = os.path.join(current_dir, 'output_data')
    prepare_profiles.main(spr_string, input_data_dir, output_dir, num_markers)
    assert os.path.exists(os.path.join(output_dir, f'profile_{spr_string}_ne.dat'))
    assert os.path.exists(os.path.join(output_dir, f'profile_{spr_string}_Te.dat'))
    assert os.path.exists(os.path.join(output_dir, f'profile_{spr_string}_Ti.dat'))
    assert os.path.exists(os.path.join(output_dir, f'ion_info_{spr_string}.dat'))
    assert os.path.exists(os.path.join(output_dir, f'{spr_string}_markers_10.dat'))
