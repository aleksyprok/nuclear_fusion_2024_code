"""
Module contains the routines for parameterizing the wall.
"""

import numpy as np

def calc_wall_centroid(r_wall, z_wall):
    """
    Calculate the centroid of the wall using formula given here:
    https://en.wikipedia.org/wiki/Centroid#Of_a_polygon
    """

    A = 0.5 * np.sum(r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1])

    CR = np.sum((r_wall[:-1] + r_wall[1:]) * (r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1]))
    CR /= 6 * A
    CZ = np.sum((z_wall[:-1] + z_wall[1:]) * (r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1]))
    CZ /= 6 * A

    return CR, CZ