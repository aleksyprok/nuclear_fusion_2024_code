"""
Module for testing the wall module.
"""
import numpy as np
from python_scripts import wall

def test_calc_wall_centroid():
    """
    Test the calc_wall_centroid function.
    """
    r_wall = np.array([0, 1, 1, 0, 0])
    z_wall = np.array([0, 0, 1, 1, 0])
    assert wall.calc_wall_centroid(r_wall, z_wall) == (0.5, 0.5)

def test_get_s_theta_from_rz():
    """
    Test the get_s_theta_from_rz function.
    """
    r_coord = np.array([0.300, 1.001, 0.700, 0.001])
    z_coord = np.array([0.001, 0.200, 0.999, 0.200])
    r_wall = np.array([0, 1, 1, 0, 0])
    z_wall = np.array([0, 0, 1, 1, 0])
    s_theta = np.array([0, 0, 0, 0, 0])
    s_theta = wall.get_s_theta_from_rz(r_coord, z_coord, r_wall, z_wall)
    s_theta_test = np.array([0.3, 1.2, 2.3, 3.8])
    assert np.allclose(s_theta, s_theta_test, atol=1e-8)
