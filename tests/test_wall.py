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