"""
Module to test the pkde module.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from python_scripts import pkde
import warnings
warnings.filterwarnings("error", category=DeprecationWarning)

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

def test_periodic_kde_1d():
    """
    Test case for the periodic_kde_1d function.

    This test case verifies that the periodic_kde_1d function returns the expected
    KDE estimate for a given set of coordinates.
    """
    x = np.array([1, 2, 3])
    x_min = 0
    x_max = 4
    weights = np.array([1, 2, 3])
    hx = 0.3
    kernel = 'gaussian'
    num_grid_points = 10**3
    pdf_fun = pkde.periodic_kde_1d(x, x_min, x_max, weights, hx,
                                   kernel = kernel,
                                   num_grid_points = num_grid_points)

    # Calcualte KDE using scipy, not FFTKDE
    x = np.array([-3, -2, -1, +1, +2, +3, +5, +6, +7])
    sigma_x = np.std(x)
    weights = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3])
    kde_scipy = stats.gaussian_kde(x,
                                   bw_method=hx/sigma_x,
                                   weights=weights)

    # Plot the KDE estimate
    grid_x = np.linspace(x_min, x_max, num_grid_points)
    fig, ax = plt.subplots()
    ax.plot(grid_x, pdf_fun(grid_x))
    ax.plot(grid_x, 3 * kde_scipy(grid_x))
    ax.set_xlabel('x')
    ax.set_ylabel('PDF')
    ax.set_title('KDE estimate')
    ax.legend(['Periodic KDE', 'Scipy KDE'])
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")
    outputpath = os.path.join(output_dir, 'periodic_kde_1d.png')
    fig.savefig(outputpath)

    # Check integral of the PDE is 1 within a 10^-4 tolerance
    assert np.isclose(np.trapz(pdf_fun(grid_x), grid_x),
                      1,
                      atol=1e-4)
    # Check the integral of scipy KDE is 1 within a 10^-4 tolerance
    assert np.isclose(np.trapz(3 * kde_scipy(grid_x), grid_x),
                      1,
                      atol=1e-4)
    # Check the KDE estimate is within 10% of the scipy KDE
    assert np.isclose(np.trapz(pdf_fun(grid_x), grid_x),
                      np.trapz(3 * kde_scipy(grid_x), grid_x),
                      atol=0.1)

def test_periodic_kde_2d():
    """
    Test case for the periodic_kde_2d function.

    This test case verifies that the periodic_kde_2d function returns the expected
    KDE estimate for a given set of coordinates.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_dir = os.path.join(repo_path, "tests", "output_plots")

    x = np.array([1, 2, 3])
    y = np.array([2, 5.6, 2])
    x_min = 0
    x_max = 4
    y_min = 1
    y_max = 5.7
    weights = np.array([1, 2, 3])
    hx = 0.3
    hy = 0.4
    kernel = 'gaussian'
    num_grid_points = 10**3
    pdf_fun = pkde.periodic_kde_2d(x, y, x_min, x_max, y_min, y_max, weights,
                                   hx, hy,
                                   kernel=kernel,
                                   num_grid_points=num_grid_points)
    grid_x = np.linspace(x_min, x_max, num_grid_points)
    grid_y = np.linspace(y_min, y_max, num_grid_points)
    pdf_fft = pdf_fun(grid_x, grid_y)

    # Plot the KDE estimate image
    fig, ax = plt.subplots()
    im = ax.imshow(pdf_fft.T, extent=[x_min, x_max, y_min, y_max],
                   origin='lower')
    ax.set_title('FFT KDE estimate')
    fig.colorbar(im)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.scatter(x, y)
    outputpath = os.path.join(output_dir, 'periodic_kde_2d.png')
    fig.savefig(outputpath)

    # Check integral is equal to 1
    assert np.isclose(np.trapz(np.trapz(pdf_fft, grid_x, axis=0), grid_y, axis =0),
                      1,
                      atol=1e-4)
