"""
This module contains routines fro calculating kernel desnity estimates
over 1D and 2D domains. With periodic boundary conditions.
"""

import numpy as np
from scipy import interpolate
from KDEpy import FFTKDE

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
def periodic_kde_1d(x, x_min, x_max, weights, hx,
                    kernel='gaussian',
                    num_grid_points=10**3):
    """
    Calculates a 1D kernel density estimate using FFTs with periodic
    boundary conditions.

    Parameters
    ----------
    x : array
        1D array of x coordinates
    x_min : float
        Minimum x coordinate
    x_max : float
        Maximum x coordinate
    weights : array
        1D array of weights
    hx : float
        Bandwidth in metres
    kernel : str
        KDE kernel
    num_grid_points : int
        Number of grid points

    Returns
    -------
    pdf_fun : function
        KDE estimate.
    """

    # Extending coordiantes
    period_x = x_max - x_min
    x_norm = np.concatenate((x - period_x, x, x + period_x))
    weights_norm = np.concatenate((weights, weights, weights))
    delta_x = period_x / 2
    reduced_indices = np.where(
        (x_min - delta_x <= x_norm) *
                           (x_norm <= x_max + delta_x))
    x_norm = x_norm[reduced_indices]
    weights_norm = weights_norm[reduced_indices]
    # Double number of grid points to account for extention
    num_grid_points *= 2

    # Stretching coordiantes
    x_norm = x_norm / hx

    # Calculating extended and stretched KDE
    fft_kde = FFTKDE(kernel = kernel, bw = 1)
    grid_norm_x, pdf = \
        fft_kde.fit(x_norm, weights = weights_norm).evaluate(num_grid_points)
    pdf *= np.sum(weights_norm) / np.sum(weights) / hx # Renormalise

    # Evaluating extended and stretched KDE on an unstretched grid
    grid_x = hx * grid_norm_x
    pdf_fun = interpolate.UnivariateSpline(grid_x, pdf,
                                           s=0, ext='raise')

    return pdf_fun

def periodic_kde_2d(x, y, x_min, x_max, y_min, y_max, weights,
                    hx, hy,
                    kernel = 'gaussian',
                    num_grid_points = 10**3):
    """
    Calculates a 2D kernel density estimate using FFTs with periodic
    boundary conditions.

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
    hx : float
        Bandwidth in metres
    hy : float
        Bandwidth in metres
    kernel : str
        KDE kernel
    num_grid_points : int
        Number of grid points

    Returns
    -------
    pdf_fun : array
        2D function  which interpolates the FFTKDE
    """

    # Extending coordiantes
    x_norm, y_norm, weights_norm = extend_coords(x, y, x_min, x_max, y_min, \
                                                 y_max, weights)
    # Double number of grid points to account for extention
    num_grid_points *= 2

    # Stretching coordiantes
    x_norm = x_norm / hx
    y_norm = y_norm / hy

    # Calculating extended and stretched KDE
    fft_kde = FFTKDE(kernel=kernel, bw=1)
    coords_norm = np.vstack([x_norm, y_norm]).T
    grid_norm, pdf = \
        fft_kde.fit(coords_norm, weights = weights_norm).evaluate(num_grid_points)
    grid_norm_x, grid_norm_y = np.unique(grid_norm[:, 0]), np.unique(grid_norm[:, 1])
    pdf = pdf.reshape(num_grid_points, num_grid_points) \
        * np.sum(weights_norm) / np.sum(weights) / (hx * hy) # Renormalise

    # Evaluating extended and stretched KDE on an unstretched grid
    # Apply optional scale factor)
    grid_x = grid_norm_x * hx
    grid_y = grid_norm_y * hy
    pdf_fun = interpolate.RectBivariateSpline(grid_x, grid_y, pdf)
    # pdf_fun = interpolate.interp2d(grid_x, grid_y, pdf, \
    #                                kind = 'cubic', bounds_error = True)
    return pdf_fun
