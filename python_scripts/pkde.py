"""
This module contains routines fro calculating kernel desnity estimates
over 1D and 2D domains. With periodic boundary conditions.
"""

import numpy as np

def extend_coords(x, y, x_min, x_max, y_min, y_max, weights):
    """
    Repeats the coordinates in each dimension by assuming periodicity
    to ensure KDE is well resolved at the boundaries.

    Parameters
    ----------
    x : array
        1D array of x coordinates
    y : array
        1D array of y coordinates
    x_min : float
        Minimum x coordinate
    x_max : float
        Maximum x coordinate
    y_min : float
        Minimum y coordinate
    y_max : float
        Maximum y coordinate
    weights : array
        1D array of weights

    Returns
    -------
    x_extended : array
        1D array of extended x coordinates
    y_extended : array 
        1D array of extended y coordinates
    weights_extended : array
        1D array of extended weights
    """

    period_x = x_max - x_min
    period_y = y_max - y_min

    # Repeat array in each dimension
    x_block =  np.concatenate((x, x, x))
    x_extended = np.concatenate((x_block - period_x, x_block, x_block + period_x))

    y_block = np.concatenate((y - period_y, y, y + period_y))
    y_extended = np.concatenate((y_block, y_block, y_block))

    weights_block = np.concatenate((weights, weights, weights))
    weights_extended = np.concatenate((weights_block, weights_block, weights_block))

    delta_x = period_x / 2
    delta_y = period_y / 2

    reduced_indices = np.where(
        (x_min - delta_x <= x_extended) *
                           (x_extended <= x_max + delta_x) *
        (y_min - delta_y <= y_extended) *
                           (y_extended <= y_max + delta_y))
    x_extended = x_extended[reduced_indices]
    y_extended = y_extended[reduced_indices]
    weights_extended = weights_extended[reduced_indices]

    return x_extended, y_extended, weights_extended
