"""
This file contains functions used to make the plots to compare results from
SPR-045-16 and SPR-046-16b.
"""

import os
import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from python_scripts import my_gfile_reader, paper_plots, run

REPOSITORY_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS_DIRECTORY = os.path.join(REPOSITORY_PATH, "output_data",
                              "spr_045_16_vs_francis_q_profile_processed")
WALL_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
NUM_GRID_POINTS_1D = 10**5
NUM_GRID_POINTS_2D = 10**3
NUM_BOOTSTRAPS = 128
CONFIDENCE_LEVEL = 0.95
H_PHI_ARRAY = np.logspace(-2, 0, 8)
H_THETA_1D_ARRAY = np.logspace(-3, 0, 64)
H_THETA_2D_ARRAY = np.logspace(-2, 0, 8)
SPECIAL_NODES = (16, 55, 138, 178)

def compare_q_profiles():
    """
    Compare the Q-profiles of SPR-045-16 vs SPR046_16b vs 
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    plot_dir = os.path.join(repo_dir, 'plots', 'francis_q_profiles_plots')
    os.makedirs(plot_dir, exist_ok=True)
    eqdsk_path = os.path.join(repo_dir, 'input_data', 'SPR-045-16.eqdsk')
    gfile_spr045_16 = my_gfile_reader.getGfile(eqdsk_path)
    eqdsk_path = os.path.join(repo_dir, 'input_data', 'SPR046_16b.eqdsk')
    gfile_spr046_16b = my_gfile_reader.getGfile(eqdsk_path)
    eqdsk_path = os.path.join(repo_dir, 'input_data', 'SPR-045-14.eqdsk')
    gfile_spr045_14 = my_gfile_reader.getGfile(eqdsk_path)

    fig, ax = plt.subplots()
    ax.plot(gfile_spr045_16.psin, gfile_spr045_16.q,
            label='SPR-045-16')
    ax.plot(gfile_spr046_16b.psin, gfile_spr046_16b.q,
            label='SPR-046-16b')
    ax.plot(gfile_spr045_14.psin, gfile_spr045_14.q,
            label='SPR-045-14')
    ax.legend()
    ax.set_xlabel(r'$\psi_n$')
    ax.set_ylabel(r'$q$')
    fig.savefig(plot_dir + '/q_profile_compare.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

def plot_all_runs():
    """
    This function plots the francis_vs_spr_045_16 runs to produce Figure 2 of the paper.
    """

    def create_csv():
        runs = all_runs

        runs_metadata = []
        for i, run_i in enumerate(runs):
            print("axisymmetric", run_i.log.axisymmetric)
            print("analytic_ripple", run_i.log.analytic_ripple)
            print("eqdsk_fname", run_i.log.eqdsk_fname)
            print("bplasma", run_i.log.bplasma)
            spr_name = run_i.log.eqdsk_fname[:-6]
            print("spr_name", spr_name)
            if run_i.log.axisymmetric:
                remap_phi_n = None
            elif run_i.log.analytic_ripple:
                remap_phi_n = run_i.log.ncoil
            elif run_i.log.bplasma:
                remap_phi_n = run_i.log.bplasma_n

            output_dir_i = os.path.join(output_dir, f"{spr_name}_axisymmetric_"
                                        f"{run_i.log.axisymmetric}_analytic_ripple_"
                                        f"{run_i.log.analytic_ripple}_bplasma_"
                                        f"{run_i.log.bplasma}")
            if i==0:
                paper_plots.calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=remap_phi_n)
            else:
                paper_plots.calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=remap_phi_n,
                                 previous_run=runs[i-1])
            runs_metadata.append([spr_name,
                                  run_i.log.axisymmetric,
                                  run_i.log.analytic_ripple,
                                  run_i.log.bplasma,
                                  run_i.flux.max_energy_1d,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_1d,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_1d,
                                  run_i.flux.h_theta_2d])
        columns = ['spr_name', 'axisymmetric', 'analytic_ripple', 'bplasma',
                   'max_energy_1d', 'max_energy_2d', 'total_energy',
                   'conf_band_1d', 'conf_band_2d', 'conf_band_total',
                   'h_phi', 'h_theta_1d', 'h_theta_2d']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'francis_vs_spr_045_16_runs.csv'))

    all_runs = run.create_runs_list(RUNS_DIRECTORY)
    output_dir = os.path.join(REPOSITORY_PATH, "plots", "francis_vs_spr_045_16_runs")
    make_csv = False
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'francis_vs_spr_045_16_runs.csv'))
    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    spr_names = ["SPR-045-14", "SPR-045-16", "SPR046_16b"]
    for spr_name in spr_names:
        df_filtered = df[(df.spr_name == spr_name)]
        attributes = ["axisymmetric", "analytic_ripple", "bplasma"]
        max_energy_flux = np.zeros(3)
        max_energy_flux_error = np.zeros(3)
        total_energy_flux = np.zeros(3)
        total_energy_flux_error = np.zeros(3)
        for i, attribute in enumerate(attributes):
            df_filtered_2 = df_filtered[df_filtered[attribute]]
            if attribute == "analytic_ripple":
                max_energy_flux[i] = df_filtered_2.max_energy_2d.values[0]
            else:
                max_energy_flux[i] = df_filtered_2.max_energy_1d.values[0]
            max_energy_flux_error[i] = df_filtered_2.conf_band_1d.values[0]
            total_energy_flux[i] = df_filtered_2.total_energy.values[0]
            total_energy_flux_error[i] = df_filtered_2.conf_band_total.values[0]


        axs[0].errorbar(range(3), max_energy_flux,
                        yerr=max_energy_flux_error,
                        linestyle='None',
                        marker='+',
                        label=spr_name)
        axs[1].errorbar(range(3), total_energy_flux,
                        yerr=total_energy_flux_error,
                        linestyle='None',
                        marker='+',
                        label=spr_name)
        xticks = ["Axisymmetric", "Ripple", "RMP"]
        for i in range(2):
            axs[i].legend()
            axs[i].set_xticks(range(3), xticks)
            axs[i].set_yscale('log')
        axs[0].set_ylabel(r'Maximum Alpha Particle Energy Flux [MW/m$^2$]')
        axs[1].set_ylabel(r'Power lost [%]')
    fig.suptitle("Maximum Alpha Particle Energy Flux vs. SPR and 3D Field \n"
                 r"(Ripple configuration: $R_{coil}$ = 7.25, $N_{coil}$ = 16,"
                 "\nRMP configuration: n=3, coil_set=exterior, current=90kAt, plasma_response=True, phase=20)",
                 y=1)
    output_path = os.path.join(output_dir, 'max_and_total_flux_vs_spr_and_3d_field')
    fig.savefig(output_path + ".pdf", bbox_inches='tight')
    fig.savefig(output_path + ".png", bbox_inches='tight',
                dpi=300)
    plt.close(fig)

if __name__ == '__main__':
    # compare_q_profiles()
    start_time = time.time()
    plot_all_runs()
    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.2e} seconds")
    