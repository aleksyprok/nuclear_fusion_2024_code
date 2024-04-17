"""
This file contains functions used to make the plots that appears in the IAEA FEC 2024
paper.

Note that the code in this file assumes that the LOCUST runs have completed and
the data files are in the output_data/FEC_2024 directory.
"""
import os
import pickle
import time
from matplotlib.colors import Normalize
# pylint: disable=no-name-in-module
from matplotlib.cm import ScalarMappable, viridis
# pylint: enable=no-name-in-module
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import pandas as pd
from python_scripts import bootstrap, flux, paper_plots_extra, paper_plots_3d, run, wall

REPOSITORY_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUNS_DIRECTORY = os.path.join(REPOSITORY_PATH, "output_data",
                                "FEC_2024")
GFILE_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPR-045-16.eqdsk")
WALL_PATH = os.path.join(REPOSITORY_PATH, "input_data", "SPP-001_wall.dat")
NUM_GRID_POINTS_1D = 10**5
NUM_GRID_POINTS_2D = 10**3
NUM_BOOTSTRAPS = 128
CONFIDENCE_LEVEL = 0.95
H_PHI_ARRAY = np.logspace(-2, 0, 8)
H_THETA_1D_ARRAY = np.logspace(-3, 0, 64)
H_THETA_2D_ARRAY = np.logspace(-2, 0, 8)
SPECIAL_NODES = (16, 55, 138, 178)

def calc_energy_flux(run_i, output_dir_i,
                     remap_phi_n=None,
                     previous_run=None):
    """
    This function calculates the energy flux for a given run.
    And plots some of the data along the way.
    """
    os.makedirs(output_dir_i, exist_ok=True)
    run_i.init_wall(WALL_PATH,
                    special_nodes=SPECIAL_NODES)
    run_i.init_markers(remap_phi_n=remap_phi_n)
    run_i.init_flux(num_grid_points_1d=NUM_GRID_POINTS_1D,
                    num_grid_points_2d=NUM_GRID_POINTS_2D,
                    num_bootstraps=NUM_BOOTSTRAPS,
                    h_phi_array=H_PHI_ARRAY,
                    h_theta_1d_array=H_THETA_1D_ARRAY,
                    h_theta_2d_array=H_THETA_2D_ARRAY)
    # run_i.flux.h_phi = 1
    # run_i.flux.h_theta_2d = 0.01
    flux.calc_optimum_bandwidth_1d(run_i)
    flux.calc_optimum_bandwidth_2d(run_i)
    if previous_run is not None:
        # If the optimum bandwidth finder fails then use the value from the previous run.
        if np.isclose(run_i.flux.h_phi, H_PHI_ARRAY[-1], atol=1e-3) and \
            np.isclose(run_i.flux.h_theta_2d, H_THETA_2D_ARRAY[-1], atol=1e-3):
            run_i.flux.h_phi = previous_run.flux.h_phi
            run_i.flux.h_theta_2d = previous_run.flux.h_theta_2d
    paper_plots_extra.plot_amise_1d(run_i, output_dir_i)
    paper_plots_extra.plot_amise_2d(run_i, output_dir_i)
    run_i.flux.energy_1d = flux.calc_energy_flux_1d(run_i)
    run_i.flux.energy_2d = flux.calc_energy_flux_2d(run_i)
    paper_plots_extra.plot_energy_flux_1d(run_i, output_dir_i)
    paper_plots_extra.plot_energy_flux_2d(run_i, output_dir_i)
    be_array_1d = bootstrap.calc_be_1d_array(run_i)
    be_array_2d = bootstrap.calc_be_2d_array(run_i)
    run_i.flux.conf_band_1d = np.quantile(be_array_1d, CONFIDENCE_LEVEL)
    run_i.flux.conf_band_2d = np.quantile(be_array_2d, CONFIDENCE_LEVEL)
    be_array_total = bootstrap.calc_be_total_array(run_i)
    run_i.flux.conf_band_total = np.quantile(be_array_total, CONFIDENCE_LEVEL)
    # run_i.flux.conf_band_2d = 0.1
    # run_i.flux.conf_band_total = 0.1
    run_i.flux.max_energy_1d = np.max(run_i.flux.energy_1d)
    run_i.flux.max_energy_2d = np.max(run_i.flux.energy_2d)
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
                calc_energy_flux(run_i, output_dir_i,
                                 remap_phi_n=run_i.log.ncoil)
            else:
                calc_energy_flux(run_i, output_dir_i,
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
        df.to_csv(os.path.join(output_dir, 'ripple_runs.csv'))

    output_dir = os.path.join(REPOSITORY_PATH, "plots", "ripple_runs")
    make_csv = False
    save_axisymmetric = False
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'ripple_runs.csv'))
    for run_i in all_runs:
        if run_i.log.axisymmetric:
            run_axisymmetric = run_i
    output_dir_i = os.path.join(output_dir, "axisymmetric")
    if save_axisymmetric:
        calc_energy_flux(run_axisymmetric, output_dir_i)
        pickle.dump(run_axisymmetric, open(os.path.join(output_dir, 'run_axisymmetric.pkl'), 'wb'))
    else:
        run_axisymmetric = pickle.load(open(os.path.join(output_dir, 'run_axisymmetric.pkl'), 'rb'))
    linestyles = ['-', '--', ':']
    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    for i, ncoil in enumerate([12, 16, 18]):
        rcoils = df[df.ncoil == ncoil].sort_values(by='rcoil').rcoil.values
        max_energy_flux = df[df.ncoil == ncoil].sort_values(by='rcoil').max_energy_flux.values
        total_energy_flux = \
            df[df.ncoil == ncoil].sort_values(by='rcoil').total_energy_flux.values / \
            run_axisymmetric.log.pinj * 100
        conf_band_2d = df[df.ncoil == ncoil].sort_values(by='rcoil').conf_band_2d.values
        conf_band_total = \
            df[df.ncoil == ncoil].sort_values(by='rcoil').conf_band_total.values / \
            run_axisymmetric.log.pinj * 100
        axs[0].errorbar(rcoils, max_energy_flux, yerr=conf_band_2d,
                        label=r"$N_{coil}$ = "f"{ncoil}",
                        linestyle=linestyles[i])
        axs[0].axhline(y=run_axisymmetric.flux.max_energy_2d,
                       color='k',
                       linestyle=':')
        axs[0].legend()
        axs[1].errorbar(rcoils, total_energy_flux, yerr=conf_band_total,
                        label=r"$N_{coil}$ = "f"{ncoil}",
                        linestyle=linestyles[i])
        axs[1].axhline(y=run_axisymmetric.flux.total_energy / run_axisymmetric.log.pinj * 100,
                       color='k',
                       linestyle=':')
        axs[1].legend()
    axs[0].set_ylabel(r'Maximum Alpha Particle Energy Flux [MW/m$^2$]')
    axs[1].set_ylabel(r'Power lost [%]')
    fig.suptitle('TF Ripple Field Results')
    for i in range(2):
        axs[i].set_xlabel(r'Major radius of TF coil outer limb ($R_{coil}$) [m]')
        axs[i].set_yscale('log')
    output_path = os.path.join(output_dir, 'max_and_total_flux_vs_rcoil')
    fig.savefig(output_path + ".pdf", bbox_inches='tight')
    fig.savefig(output_path + ".png", bbox_inches='tight',
                dpi=300)
    plt.close(fig)

def plot_rmp_runs(all_runs):
    """
    This function plots the RMP runs to produce Figures 3 and 4 of the paper.
    """

    def create_csv():
        # First we need to filter the rmp runs
        runs = []
        for run_i in all_runs:
            if run_i.log.bplasma:
                if run_i.log.bplasma_n > 1:
                    runs.append(run_i)
        # Next sort the runs by rcoil and ncoil
        runs.sort(key=lambda x: x.log.rmp_phase, reverse=False)
        runs.sort(key=lambda x: x.log.rmp_response, reverse=False)
        runs.sort(key=lambda x: x.log.rmp_current, reverse=False)
        runs.sort(key=lambda x: x.log.bplasma_n, reverse=False)
        runs.sort(key=lambda x: x.log.coil_set, reverse=True)
        for run_i in runs:
            print(run_i.log.rmp_phase, run_i.log.rmp_response, run_i.log.rmp_current,
                  run_i.log.bplasma_n, run_i.log.coil_set)

        runs_metadata = []
        for run_i in runs:
            output_dir_i = os.path.join(output_dir,
                                        f"{run_i.log.coil_set}",
                                        f"{run_i.log.bplasma_n}",
                                        f"{run_i.log.rmp_current}_{run_i.log.rmp_response}",
                                        f"{run_i.log.rmp_phase}")
            # calc_energy_flux(run_i, output_dir_i,
            #                  remap_phi_n=run_i.log.bplasma_n)
            calc_energy_flux(run_i, output_dir_i,
                             remap_phi_n=None)
            runs_metadata.append([run_i.log.coil_set,
                                  run_i.log.bplasma_n,
                                  run_i.log.rmp_current,
                                  run_i.log.rmp_response,
                                  run_i.log.rmp_phase,
                                  run_i.flux.max_energy_1d,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_1d,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_1d,
                                  run_i.flux.h_theta_2d])
        columns = ['coil_set', 'bplasma_n', 'rmp_current', 'rmp_response', 'rmp_phase',
                   'max_energy_flux_1d', 'max_energy_flux_2d', 'total_energy_flux',
                   'conf_band_1d', 'conf_band_2d', 'conf_band_total', 'h_phi', 'h_theta_1d',
                   'h_theta_2d']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'rmp_runs.csv'))

    output_dir = os.path.join(REPOSITORY_PATH, "plots", "rmp_runs")
    make_csv = False
    save_axisymmetric = False
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'rmp_runs.csv'))
    for run_i in all_runs:
        if run_i.log.axisymmetric:
            run_axisymmetric = run_i
    output_dir_i = os.path.join(output_dir, "axisymmetric")
    if save_axisymmetric:
        calc_energy_flux(run_axisymmetric, output_dir_i)
        pickle.dump(run_axisymmetric, open(os.path.join(output_dir, 'run_axisymmetric.pkl'), 'wb'))
    else:
        run_axisymmetric = pickle.load(open(os.path.join(output_dir, 'run_axisymmetric.pkl'), 'rb'))
    linestyles = ['--', '-', '--', '-']
    clrs = ['tab:blue', 'tab:blue', 'tab:orange', 'tab:orange']
    for coil_set in ["interior_rmp", "exterior_rmp"]:
        fig, axs = plt.subplots(3, 2)
        fig_size = fig.get_size_inches()
        fig_size[0] *= 2
        fig_size[1] *= 3
        fig.set_size_inches(fig_size)
        for i, bplasma_n in enumerate([2, 3, 4]):
            if coil_set == "interior_rmp":
                if bplasma_n == 2:
                    current0 = 30
                    optimum_phase = 265
                elif bplasma_n == 3:
                    current0 = 50
                    optimum_phase = 173
                elif bplasma_n == 4:
                    current0 = 80
                    optimum_phase = 67
            elif coil_set == "exterior_rmp":
                if bplasma_n == 2:
                    current0 = 50
                    optimum_phase = 61
                elif bplasma_n == 3:
                    current0 = 90
                    optimum_phase = 20
                elif bplasma_n == 4:
                    current0 = 150
                    optimum_phase = 321
            parameters = []
            for rmp_current in [current0, 2 * current0]:
                for rmp_response in [False, True]:
                    parameters.append((rmp_current, rmp_response))
            for j, (rmp_current, rmp_response) in enumerate(parameters):
                df_i = df.loc[(df['coil_set'] == coil_set) & (df['bplasma_n'] == bplasma_n) &
                              (df['rmp_current'] == rmp_current) &
                              (df['rmp_response'] == rmp_response)]
                rmp_phases = df_i['rmp_phase'].values
                max_energy_flux = df_i.max_energy_flux_1d.values
                total_energy_flux = df_i.total_energy_flux.values / run_axisymmetric.log.pinj * 100
                conf_band_1d = df_i.conf_band_1d.values
                conf_band_total = df_i.conf_band_total.values / run_axisymmetric.log.pinj * 100
                axs[i, 0].errorbar(rmp_phases, max_energy_flux, yerr=conf_band_1d,
                                   linestyle=linestyles[j], color=clrs[j])
                axs[i, 0].axhline(y=run_axisymmetric.flux.max_energy_1d,
                                  color='k',
                                  linestyle=':')
                axs[i, 0].axvline(x=optimum_phase, color='k', linestyle=':')
                axs[i, 1].errorbar(rmp_phases, total_energy_flux, yerr=conf_band_total,
                                   linestyle=linestyles[j], color=clrs[j])
                axisymmetric_total_energy_flux = run_axisymmetric.flux.total_energy \
                                               / run_axisymmetric.log.pinj * 100
                axs[i, 1].axhline(y=axisymmetric_total_energy_flux,
                                  color='k',
                                  linestyle=':')
                axs[i, 1].axvline(x=optimum_phase, color='k', linestyle=':')
                axs[i, 0].get_shared_y_axes().joined(axs[i, 0], axs[0, 0])
                axs[i, 1].get_shared_y_axes().joined(axs[i, 1], axs[0, 1])
            blue_line = mlines.Line2D([], [],
                                      color='tab:blue',
                                      label=f'Current: {current0} kAt')
            orange_line = mlines.Line2D([], [],
                                        color='tab:orange',
                                        label=f'Current: {2 * current0} kAt')
            dashed_line = mlines.Line2D([], [],
                                        color='black',
                                        linestyle='-',
                                        label='Plasma response')
            solid_line = mlines.Line2D([], [],
                                       color='black',
                                       linestyle='--',
                                       label='Vacuum')
            leg1 = axs[i, 0].legend(handles=[orange_line, blue_line])
            leg2 = axs[i, 1].legend(handles=[dashed_line, solid_line])
            axs[i, 0].add_artist(leg1)
            axs[i, 1].add_artist(leg2)
            for j in range(2):
                axs[i, j].set_title(f'n = {bplasma_n}')
                axs[i, j].set_yscale('log')
            axs[i, 0].set_ylabel(r'Maximum Alpha Particle Energy Flux [MW m$^{-2}$]')
            output_path = os.path.join(output_dir,
                                       f"max_and_total_flux_vs_phase_{coil_set}")
            if coil_set == "interior_rmp":
                fig.suptitle("In-Vessel ELM Suppression Results",
                             y=0.92)
            elif coil_set == "exterior_rmp":
                fig.suptitle("Out-of-Vessel ELM Suppression Results",
                             y=0.92)
            fig.savefig(output_path + '.pdf', bbox_inches='tight')
            fig.savefig(output_path + '.png', bbox_inches='tight',
                        dpi=300)
            plt.close(fig)

def plot_rmp_distribution(all_runs):
    """
    This function plots the the energy distribution runs to produce Figure 5 of the paper.
    """
    for run_i in all_runs:
        if run_i.log.bplasma and \
           run_i.log.coil_set == "exterior_rmp" and \
           run_i.log.bplasma_n == 3 and \
           run_i.log.rmp_current == 90 and \
           run_i.log.rmp_response and \
           run_i.log.rmp_phase == 20:
            print("Run found!")
            run0 = run_i
    print(run0.log.bplasma_n, run0.log.rmp_current, run0.log.rmp_response,
            run0.log.rmp_phase)
    run0.init_wall(WALL_PATH,
                    special_nodes=SPECIAL_NODES)
    run0.init_markers(remap_phi_n=None)
    run0.init_flux(num_grid_points_1d=NUM_GRID_POINTS_1D,
                   num_grid_points_2d=NUM_GRID_POINTS_2D,
                   num_bootstraps=NUM_BOOTSTRAPS,
                   h_phi_array=H_PHI_ARRAY,
                   h_theta_1d_array=H_THETA_1D_ARRAY,
                   h_theta_2d_array=H_THETA_2D_ARRAY)
    flux.calc_optimum_bandwidth_1d(run0)
    run0.flux.energy_1d = flux.calc_energy_flux_1d(run0)
    # Make copies of run0.flux.energy_1d and run0.flux.s_theta_1d
    # with shorter length using interpolation
    new_length = 1000
    num_edges = new_length - 1
    new_s_theta = np.interp(np.linspace(0, 1, new_length),
                            np.linspace(0, 1, len(run0.flux.s_theta_1d)),
                            run0.flux.s_theta_1d)
    new_energy_flux = np.interp(np.linspace(0, 1, new_length),
                                np.linspace(0, 1, len(run0.flux.energy_1d)),
                                run0.flux.energy_1d)
    x_points, y_points = wall.get_rz_from_s_theta(new_s_theta,
                                                  run0.wall.r,
                                                  run0.wall.z)
    new_energy_flux_mid = 0.5 * (new_energy_flux[1:] + new_energy_flux[:-1])

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)

    norm = Normalize(vmin=0, vmax=np.max(run0.flux.energy_1d))
    scalar_map = ScalarMappable(norm=norm, cmap=viridis,)
    for i in range(num_edges):
        color = scalar_map.to_rgba(new_energy_flux_mid[i])
        axs[0].plot(x_points[i:i+2], y_points[i:i+2],
                    color=color)
    plt.colorbar(scalar_map, ax=axs[0], orientation='vertical')
    axs[0].set_aspect('equal')
    axs[0].set_xlim(0, 7.5)
    axs[0].set_ylim(-10, 10)
    axs[0].set_xlabel('R [m]')
    axs[0].set_ylabel('Z [m]')
    clrs = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    symbols = ['+', 'x', 'o', 's']
    axs[1].plot(run0.flux.s_theta_1d, run0.flux.energy_1d, 'k')
    y_mid = 0.5 * (axs[1].get_ylim()[1] - axs[1].get_ylim()[0])
    for i, s_nod_i in enumerate(run0.wall.special_nodes):
        axs[0].plot(run0.wall.r[s_nod_i], run0.wall.z[s_nod_i],
                    color=clrs[i],
                    marker=symbols[i])
        axs[1].axvline(x=run0.wall.s_nodes[s_nod_i],
                       color=clrs[i])
        axs[1].plot(run0.wall.s_nodes[s_nod_i], y_mid,
                    color=clrs[i],
                    marker=symbols[i])

    axs[1].set_xlabel(r'Poloidal Distance along the Inner Wall ($s_\theta$) [m]')
    fig.suptitle(r'Maximum Alpha Particle Energy Flux [MW m$^{-2}$]',
                 x=0.6)
    output_dir = os.path.join(REPOSITORY_PATH, "plots")
    output_path = os.path.join(output_dir,
                               "energy_flux_distribution_rmp")
    fig.savefig(output_path + '.pdf', bbox_inches='tight')
    fig.savefig(output_path + '.png', bbox_inches='tight',
                dpi=300)
    plt.close(fig)

def plot_rwm_runs(all_runs):
    """
    This function plots the RWM runs to produce Figure 6 of the paper.
    """

    def create_csv():
        # First we need to filter the rwm runs
        runs = []
        for run_i in all_runs:
            if run_i.log.bplasma:
                if run_i.log.bplasma_n == 1:
                    runs.append(run_i)
        # Next sort the runs by bscale
        runs.sort(key=lambda x: x.log.rwm_bscale, reverse=False)

        runs_metadata = []
        for run_i in runs:
            output_dir_i = os.path.join(output_dir,
                                        f"bscale_{run_i.log.rwm_bscale}")
            calc_energy_flux(run_i, output_dir_i)
            runs_metadata.append([run_i.log.rwm_bscale,
                                  run_i.flux.max_energy_2d,
                                  run_i.flux.total_energy,
                                  run_i.flux.conf_band_2d,
                                  run_i.flux.conf_band_total,
                                  run_i.flux.h_phi,
                                  run_i.flux.h_theta_2d])
        columns = ['bscale', 'max_energy_flux', 'total_energy_flux', 'conf_band_2d',
                   'conf_band_total', 'h_phi', 'h_theta_2d']
        df = pd.DataFrame(runs_metadata, columns=columns)
        df.to_csv(os.path.join(output_dir, 'rwm_runs.csv'))

    output_dir = os.path.join(REPOSITORY_PATH, "plots", "rwm_runs")
    make_csv = False
    save_axisymmetric = False
    if make_csv:
        create_csv()
    df = pd.read_csv(os.path.join(output_dir, 'rwm_runs.csv'))
    output_dir_i = os.path.join(output_dir, "axisymmetric")
    if save_axisymmetric:
        for run_i in all_runs:
            if run_i.log.axisymmetric:
                run_axisymmetric = run_i
        calc_energy_flux(run_axisymmetric, output_dir_i)
        pickle.dump(run_axisymmetric,
                    open(os.path.join(output_dir_i, 'axisymmetric_run.pickle'), 'wb'))
    else:
        run_axisymmetric = pickle.load(open(os.path.join(output_dir_i, 
                                                         'axisymmetric_run.pickle'), 'rb'))

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)
    bscales = df.sort_values(by='bscale').bscale.values
    max_energy_flux = df.sort_values(by='bscale').max_energy_flux.values
    total_energy_flux = df.sort_values(by='bscale').total_energy_flux.values \
                      / run_axisymmetric.log.pinj * 100
    conf_band_2d = df.sort_values(by='bscale').conf_band_2d.values
    conf_band_total = df.sort_values(by='bscale').conf_band_total.values \
                    / run_axisymmetric.log.pinj * 100
    axs[0].errorbar(bscales, max_energy_flux, yerr=conf_band_2d)
    axs[1].errorbar(bscales, total_energy_flux, yerr=conf_band_total)
    axs[0].set_ylabel(r'Maximum Alpha Particle Energy Flux [MW m$^{-2}$]')
    axs[1].set_ylabel('Power lost [%]')
    for i in range(2):
        axs[i].set_xlabel('Magnetic field strength of the RWM at the sensors [T]')
        axs[i].set_xscale('log')
        axs[i].set_yscale('log')
    fig.suptitle('RWM field results')

    output_path = os.path.join(output_dir, 'max_and_total_flux_vs_bscale')
    fig.savefig(output_path + '.pdf', bbox_inches='tight')
    fig.savefig(output_path + '.png', bbox_inches='tight',
                dpi=300)
    plt.close(fig)

if __name__ == "__main__":
    start_time = time.time()
    RUNS = run.create_runs_list(RUNS_DIRECTORY)
    paper_plots_3d.coil_plot_3d(gfile_path=GFILE_PATH)
    plot_ripple_runs(RUNS)
    plot_rmp_runs(RUNS)
    plot_rmp_distribution(RUNS)
    plot_rwm_runs(RUNS)
    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.2e} seconds")
