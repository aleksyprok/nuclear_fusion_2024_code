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
from python_scripts import bootstrap, flux, paper_plots_extra, run

REPOSITORY_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS_DIRECTORY = os.path.join(REPOSITORY_PATH, "output_data",
                                "FEC_2024")
# GFILE_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPR-045-16.eqdsk")
WALL_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
NUM_GRID_POINTS_1D = 10**5
NUM_GRID_POINTS_2D = 10**3
NUM_BOOTSTRAPS = 128
CONFIDENCE_LEVEL = 0.95
H_PHI_ARRAY = np.logspace(-2, 0, 8)
H_THETA_1D_ARRAY = np.logspace(-2, 0, 64)
H_THETA_2D_ARRAY = np.logspace(-2, 0, 8)
SPECIAL_NODES = (16, 55, 138, 178)

def calc_energy_flux(run_i, output_dir_i,
                     previous_run=None):
    """
    This function calculates the energy flux for a given run.
    And plots some of the data along the way.
    """
    os.makedirs(output_dir_i, exist_ok=True)
    run_i.init_wall(WALL_PATH,
                    special_nodes=SPECIAL_NODES)
    run_i.init_markers(remap_phi_n=run_i.log.ncoil)
    run_i.init_flux(num_grid_points_1d=NUM_GRID_POINTS_1D,
                    num_grid_points_2d=NUM_GRID_POINTS_2D,
                    num_bootstraps=NUM_BOOTSTRAPS,
                    h_phi_array=H_PHI_ARRAY,
                    h_theta_1d_array=H_THETA_1D_ARRAY,
                    h_theta_2d_array=H_THETA_2D_ARRAY)
    flux.calc_optimum_bandwidth_2d(run_i)
    # If the optimum bandwidth finder fails then use the value from the previous run.
    if np.isclose(run_i.flux.h_phi, H_PHI_ARRAY[-1], atol=1e-3) and \
        np.isclose(run_i.flux.h_theta_2d, H_THETA_2D_ARRAY[-1], atol=1e-3):
        run_i.flux.h_phi = previous_run.flux.h_phi
        run_i.flux.h_theta_2d = previous_run.flux.h_theta_2d
    paper_plots_extra.plot_amise_2d(run_i, output_dir_i)
    run_i.flux.energy_2d = flux.calc_energy_flux_2d(run_i)
    paper_plots_extra.plot_energy_flux_2d(run_i, output_dir_i)
    be_array_2d = bootstrap.calc_be_2d_array(run_i)
    run_i.flux.conf_band_2d = np.quantile(be_array_2d, CONFIDENCE_LEVEL)
    be_array_total = bootstrap.calc_be_total_array(run_i)
    run_i.flux.conf_band_total = np.quantile(be_array_total, CONFIDENCE_LEVEL)
    run_i.free_space()
    paper_plots_extra.save_attributes_to_file(run_i, output_dir_i,
                                                indent=4)

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
        # Next sort the runs by rcoil and ncoil
        runs.sort(key=lambda x: (x.log.ncoil, x.log.rcoil))

        runs_metadata = []
        for i, run_i in enumerate(runs):
            output_dir_i = os.path.join(output_dir,
                                        f"rcoil_{run_i.log.rcoil}_ncoil_{run_i.log.ncoil}")
            if i==0:
                calc_energy_flux(run_i, output_dir_i)
            else:
                calc_energy_flux(run_i, output_dir_i, previous_run=runs[i-1])
            run_i.flux.max_energy_2d = np.max(run_i.flux.energy_2d)
            runs_metadata.append([run_i.log.ncoil,
                                  run_i.log.rcoil,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_2d])
        columns = ['ncoil', 'rcoil', 'max_energy_flux', 'total_energy_flux', 'conf_band_2d',
                   'conf_band_total', 'h_phi', 'h_theta_2d']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'ripple_runs.csv'))

    output_dir = os.path.join(REPOSITORY_PATH, "plots", "ripple_runs")
    make_csv = True
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'ripple_runs.csv'))
    for run_i in all_runs:
        if run_i.log.axisymmetric:
            run_axisymmetric = run_i
    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    for ncoil in [12, 16, 18]:
        rcoils = df[df.ncoil == ncoil].sort_values(by='rcoil').rcoil.values
        max_energy_flux = df[df.ncoil == ncoil].sort_values(by='rcoil').max_energy_flux.values
        total_energy_flux = df[df.ncoil == ncoil].sort_values(by='rcoil').total_energy_flux.values
        conf_band_2d = df[df.ncoil == ncoil].sort_values(by='rcoil').conf_band_2d.values
        conf_band_total = df[df.ncoil == ncoil].sort_values(by='rcoil').conf_band_total.values
        axs[0].errorbar(rcoils, max_energy_flux, yerr=conf_band_2d)
        total_energy_flux = total_energy_flux / run_axisymmetric.log.pinj * 100
        axs[1].errorbar(rcoils, total_energy_flux, yerr=conf_band_total)
    axs[0].set_ylabel(r'Maximum Alpha Particle Energy Flux [MW/m$^2$]')
    axs[1].set_ylabel(r'Power lost [%]')
    fig.suptitle('TF Ripple Field Results')
    for i in range(2):
        axs[i].set_xlabel('Rcoil')
        axs[i].set_yscale('log')
    plt.show()

if __name__ == "__main__":
    RUNS = run.create_runs_list(RUNS_DIRECTORY)
    plot_ripple_runs(RUNS)
