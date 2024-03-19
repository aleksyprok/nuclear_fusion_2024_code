"""
Module to test the pkde module.
"""
import numpy as np
from python_scripts import pkde

def test_extend_coords():
    """
    Test case for the extend_coords function.

    This test case verifies that the extend_coords function returns the expected
    extended coordinates for a given set of coordinates.
    """
    x = np.array([1, 2, 3])
    y = np.array([3, 4, 5])
    x_min = 0
    x_max = 4
    y_min = 2
    y_max = 6
    weights = np.array([1, 2, 3])
    x_extended, y_extended, weights_extended = \
        pkde.extend_coords(x, y, x_min, x_max, y_min, y_max, weights)
    
    # x_test = np.array([-3, -2, -1,
    #                    -3, -2, -1,
    #                    -3, -2, -1,
    #                    +1, +2, +3,
    #                    +1, +2, +3,
    #                    +1, +2, +3,
    #                    +5, +6, +7,
    #                    +5, +6, +7,
    #                    +5, +6, +7])
    # y_test = np.array([-1, +0, 1,
    #                    +3, +4, +5,
    #                    +7, +8, +9,
    #                    -1, +0, +1,
    #                    +3, +4, +5,
    #                    +7, +8, +9,
    #                    -1, +0, +1,
    #                    +3, +4, +5,
    #                    +7, +8, +9])
    # weights_test = np.array([1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3,
    #                          1, 2, 3])
    # Remove indices where x < x_min - delta_x or x > x_max + delta_x
    # or y < y_min - delta_y or y > y_max + delta_y
    # where delta_x = (x_max - x_min) / 2 and delta_y = (y_max - y_min) / 2
    # hence delta_x = 2 and delta_y = 2
    # hence remove elements where x = -3 or 7 and y = -1 or 9
    x_test = np.array([-2, -1,
                       -2, -1,
                       -2,
                       +2, +3,
                       +1, +2, +3,
                       +1, +2,
                       +6,
                       +5, +6,
                       +5, +6])
    y_test = np.array([+0, +1,
                       +4, +5,
                       +8,
                       +0, +1,
                       +3, +4, +5,
                       +7, +8,
                       +0,
                       +3, +4,
                       +7, +8])
    weights_test = np.array([2, 3,
                             2, 3,
                             2,
                             2, 3,
                             1, 2, 3,
                             1, 2,
                             2,
                             1, 2,
                             1, 2])

    assert np.array_equal(x_extended, x_test)
    assert np.array_equal(y_extended, y_test)
    assert np.array_equal(weights_extended, weights_test)
