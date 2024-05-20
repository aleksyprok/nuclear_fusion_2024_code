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
            output_dir_i = os.path.join(output_dir,
                                        f"rcoil_{run_i.log.rcoil}_ncoil_{run_i.log.ncoil}")
            if i==0:
                paper_plots.calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=run_i.log.ncoil)
            else:
                paper_plots.calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=run_i.log.ncoil,
                                 previous_run=runs[i-1])
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
        df.to_csv(os.path.join(output_dir, 'francis_vs_spr_045_16_runs.csv'))

    all_runs = run.create_runs_list(RUNS_DIRECTORY)
    output_dir = os.path.join(REPOSITORY_PATH, "plots", "francis_vs_spr_045_16_runs")
    make_csv = True
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'francis_vs_spr_045_16_runs.csv'))
    # linestyles = ['-', '--', ':']
    # fig, axs = plt.subplots(1, 2)
    # fig_size = fig.get_size_inches()
    # fig_size[0] *= 2
    # fig.set_size_inches(fig_size)
    # for i, ncoil in enumerate([12, 16, 18]):
    #     rcoils = df[df.ncoil == ncoil].sort_values(by='rcoil').rcoil.values
    #     max_energy_flux = df[df.ncoil == ncoil].sort_values(by='rcoil').max_energy_flux.values
    #     total_energy_flux = \
    #         df[df.ncoil == ncoil].sort_values(by='rcoil').total_energy_flux.values / \
    #         run_axisymmetric.log.pinj * 100
    #     conf_band_2d = df[df.ncoil == ncoil].sort_values(by='rcoil').conf_band_2d.values
    #     conf_band_total = \
    #         df[df.ncoil == ncoil].sort_values(by='rcoil').conf_band_total.values / \
    #         run_axisymmetric.log.pinj * 100
    #     axs[0].errorbar(rcoils, max_energy_flux, yerr=conf_band_2d,
    #                     label=r"$N_{coil}$ = "f"{ncoil}",
    #                     linestyle=linestyles[i])
    #     axs[0].axhline(y=run_axisymmetric.flux.max_energy_2d,
    #                    color='k',
    #                    linestyle=':')
    #     axs[0].legend()
    #     axs[1].errorbar(rcoils, total_energy_flux, yerr=conf_band_total,
    #                     label=r"$N_{coil}$ = "f"{ncoil}",
    #                     linestyle=linestyles[i])
    #     axs[1].axhline(y=run_axisymmetric.flux.total_energy / run_axisymmetric.log.pinj * 100,
    #                    color='k',
    #                    linestyle=':')
    #     axs[1].legend()
    # axs[0].set_ylabel(r'Maximum Alpha Particle Energy Flux [MW/m$^2$]')
    # axs[1].set_ylabel(r'Power lost [%]')
    # fig.suptitle('TF francis_vs_spr_045_16 Field Results')
    # for i in range(2):
    #     axs[i].set_xlabel(r'Major radius of TF coil outer limb ($R_{coil}$) [m]')
    #     axs[i].set_yscale('log')
    # output_path = os.path.join(output_dir, 'max_and_total_flux_vs_rcoil')
    # fig.savefig(output_path + ".pdf", bbox_inches='tight')
    # fig.savefig(output_path + ".png", bbox_inches='tight',
    #             dpi=300)
    # plt.close(fig)

if __name__ == '__main__':
    # compare_q_profiles()
    start_time = time.time()
    plot_all_runs()
    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.2e} seconds")
    