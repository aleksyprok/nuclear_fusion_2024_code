"""
This file contains functions used to make the plots that appears in the IAEA FEC 2024
paper.

Note that the code in this file assumes that the LOCUST runs have completed and
the data files are in the output_data/FEC_2024 directory.
"""
import os
import numpy as np
from python_scripts import flux, run

REPOSITORY_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS_DIRECTORY = os.path.join(REPOSITORY_PATH, "output_data",
                                "FEC_2024")
GFILE_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPR-045-16.eqdsk")
WALL_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
NUM_GRID_POINTS = 10**3
H_PHI_ARRAY = np.logspace(-2, 0, 8)
H_THETA_1D_ARRAY = np.logspace(-2, 0, 64)
H_THETA_2D_ARRAY = np.logspace(-2, 0, 8)

def plot_ripple_runs(all_runs):
    """
    This function plots the ripple runs to produce Figure 2 of the paper.
    """

    # First we need to filter the ripple runs
    runs = []
    for run_i in all_runs:
        if run_i.log.analytic_ripple:
            runs.append(run_i)
        if run_i.log.axisymmetric:
            run_axisymmetric = run_i

    # Next sort the runs by rcoil and ncoil
    runs.sort(key=lambda x: (x.log.ncoil, x.log.rcoil))

    rcoils_unique = np.linspace(7, 9, 9)
    ncoils_unique = np.array([12, 16, 18])
    optimum_h_phi = np.zeros(len(runs))
    optimum_h_theta_1d = np.zeros(len(runs))
    optimum_h_theta_2d = np.zeros(len(runs))
    for i, run_i in enumerate(runs):
        run_i.init_wall(WALL_PATH)
        run_i.init_gfile(GFILE_PATH)
        run_i.init_markers()
        run_i.init_flux(num_grid_points=NUM_GRID_POINTS)
        flux.calc_optimum_bandwidth_2d(run_i)
        # Plot amise to check it is correct
        


if __name__ == "__main__":
    RUNS = run.create_runs_list(RUNS_DIRECTORY)
    plot_ripple_runs(RUNS)
