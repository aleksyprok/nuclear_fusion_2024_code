"""
This file contains functions used to make the plots that appears in the IAEA FEC 2024
paper.

Note that the code in this file assumes that the LOCUST runs have completed and
the data files are in the output_data/FEC_2024 directory.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from python_scripts import flux, paper_plots_extra, run

REPOSITORY_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS_DIRECTORY = os.path.join(REPOSITORY_PATH, "output_data",
                                "FEC_2024")
# GFILE_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPR-045-16.eqdsk")
WALL_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
NUM_GRID_POINTS = 10**3
H_PHI_ARRAY = np.logspace(-2, 0, 8)
H_THETA_1D_ARRAY = np.logspace(-2, 0, 64)
H_THETA_2D_ARRAY = np.logspace(-2, 0, 8)
SPECIAL_NODES = (16, 55, 138, 178)

def plot_ripple_runs(all_runs):
    """
    This function plots the ripple runs to produce Figure 2 of the paper.
    """

    def create_csv():
        # First we need to filter the ripple runs
        runs = []
        for run_i in all_runs:
            if run_i.log.analytic_ripple:
                runs.append(run_i)
            if run_i.log.axisymmetric:
                run_axisymmetric = run_i
        # Next sort the runs by rcoil and ncoil
        runs.sort(key=lambda x: (x.log.ncoil, x.log.rcoil))

        runs_metadata = []
        for i, run_i in enumerate(runs):
            output_dir_i = os.path.join(output_dir, f"rcoil_{run_i.log.rcoil}_ncoil_{run_i.log.ncoil}")
            os.makedirs(output_dir_i, exist_ok=True)
            run_i.init_wall(WALL_PATH,
                            special_nodes=SPECIAL_NODES)
            # run_i.init_gfile(GFILE_PATH)
            run_i.init_markers(remap_phi_n=run_i.log.ncoil)
            run_i.init_flux(num_grid_points=NUM_GRID_POINTS,
                            h_phi_array=H_PHI_ARRAY,
                            h_theta_1d_array=H_THETA_1D_ARRAY,
                            h_theta_2d_array=H_THETA_2D_ARRAY)
            flux.calc_optimum_bandwidth_2d(run_i)
            # Plot amise to check it is correct
            paper_plots_extra.plot_amise_2d(run_i, output_dir_i)
            # run_i.flux.h_phi = 1
            # run_i.flux.h_theta_2d = 0.1
            # Calculate the energy flux
            flux.calc_energy_flux_2d(run_i)
            # Plot the energy flux
            paper_plots_extra.plot_energy_flux_2d(run_i, output_dir_i)
            run_i.free_space()
            paper_plots_extra.save_attributes_to_file(run_i, output_dir_i,
                                                      indent=4)
            runs_metadata.append([run_i.log.ncoil, run_i.log.rcoil, np.max(run_i.flux.energy_2d),
                                  run_i.total_energy_flux, run_i.flux.h_phi, run_i.flux.h_theta_2d])

        columns = ['ncoil', 'rcoil', 'max_energy_flux', 'total_energy_flux', 'h_phi', 'h_theta_2d']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'ripple_runs.csv'))

    output_dir = os.path.join(REPOSITORY_PATH, "plots", "ripple_runs")
    make_csv = True
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'ripple_runs.csv'))
    for ncoil in [12, 16, 18]:
        # Make an array called max_energy_flux for each ncoil from df and sort them in ascending order of rcoil
        max_energy_flux = df[df.ncoil == ncoil].sort_values(by='rcoil').max_energy_flux.values
        rcoils = df[df.ncoil == ncoil].sort_values(by='rcoil').rcoil.values
        print(max_energy_flux)
        print(rcoils)
        fig, ax = plt.subplots()
        ax.plot(max_energy_flux)
        plt.show()

if __name__ == "__main__":
    RUNS = run.create_runs_list(RUNS_DIRECTORY)
    plot_ripple_runs(RUNS)
