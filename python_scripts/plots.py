"""
This module contains routines for plotting results.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from python_scripts import run

def plot_total_flux_over_runs(runs_dir):
    """
    This function plots the total flux over a range of runs.
    Note this reads in the LOG*.out files not the FINAL_STATE*.dat files
    which give slightly different results.
    """

    # repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # runs_dir = os.path.join(repo_path, "output_data", "FEC_2024_missing_31_45_processed")
    runs = run.create_runs_list(runs_dir)
    for run_i in runs:
        run_i.update_log()

    print("len(runs) = ", len(runs))

    runs_axisymmetric = []
    for run_i in runs:
        if run_i.log.axisymmetric:
            runs_axisymmetric.append(run_i)
    assert len(runs_axisymmetric) == 1
    runs_axisymmetric = runs_axisymmetric[0]
    total_fluxes_ratio_axisymmetric = runs_axisymmetric.log.total_stopped_power \
                                    / runs_axisymmetric.log.pinj * 100
    print(total_fluxes_ratio_axisymmetric)

    # Ripple plot
    ripple_runs = []
    for run_i in runs:
        if run_i.log.analytic_ripple:
            ripple_runs.append(run_i)
    ripple_runs.sort(key=lambda x: x.log.rcoil)
    ripple_runs.sort(key=lambda x: x.log.ncoil)
    print("len(ripple_runs) = ", len(ripple_runs))
    rcoils = np.arange(7, 9.25, 0.25)
    ncoils = np.array([12, 16, 18])
    total_fluxes_ratio = np.zeros((len(rcoils), len(ncoils)))
    for run_n in ripple_runs:
        i = np.where(rcoils == run_n.log.rcoil)[0][0]
        j = np.where(ncoils == run_n.log.ncoil)[0][0]
        total_fluxes_ratio[i, j] = run_n.log.total_stopped_power \
                                 / run_n.log.pinj * 100
    fig, ax = plt.subplots()
    for j, ncoil in enumerate(ncoils):
        ax.plot(rcoils, total_fluxes_ratio[:, j],
                label=str(ncoil))
    ax.set_xlabel("Rcoil [m]")
    ax.set_ylabel("Power lost [%]")
    ax.set_yscale("log")
    ax.legend()

    # RMP plot
    rmp_runs = []
    for run_i in runs:
        if run_i.log.bplasma:
            if run_i.log.bplasma_n > 1:
                rmp_runs.append(run_i)
    rmp_runs.sort(key=lambda x: x.log.rmp_phase)
    rmp_runs.sort(key=lambda x: x.log.rmp_current)
    rmp_runs.sort(key=lambda x: x.log.bplasma_n)
    rmp_runs.sort(key=lambda x: x.log.coil_set)
    print("len(rmp_runs) = ", len(rmp_runs))
    print("Expected num_rmp_runs = ", 9 * 2 * 3 * 2 * 2)
    bplasma_ns = np.array([2, 3, 4])
    rmp_responses = np.array([0, 1])
    coil_sets = np.array(["exterior_rmp", "interior_rmp"])

    total_fluxes_ratio = np.zeros((9, 2, 2, len(bplasma_ns), len(coil_sets)))
    for run_n in rmp_runs:
        j = np.where(rmp_responses == run_n.log.rmp_response)[0][0]
        l = np.where(bplasma_ns == run_n.log.bplasma_n)[0][0]
        m = np.where(coil_sets == run_n.log.coil_set)[0][0]
        if run_n.log.coil_set == "exterior_rmp":
            if run_n.log.bplasma_n == 2:
                rmp_current0 = 50
                phase0 = 61
            elif run_n.log.bplasma_n == 3:
                rmp_current0 = 90
                phase0 = 20
            elif run_n.log.bplasma_n == 4:
                rmp_current0 = 150
                phase0 = 321
        elif run_n.log.coil_set == "interior_rmp":
            if run_n.log.bplasma_n == 2:
                rmp_current0 = 30
                phase0 = 265
            elif run_n.log.bplasma_n == 3:
                rmp_current0 = 50
                phase0 = 173
            elif run_n.log.bplasma_n == 4:
                rmp_current0 = 80
                phase0 = 67
        rmp_currents = np.array([rmp_current0, 2 * rmp_current0])
        k = np.where(rmp_currents == run_n.log.rmp_current)[0][0]
        rmp_phases = np.arange(0, 360, 45)
        # append phase0
        rmp_phases = np.append(rmp_phases, phase0)
        i = np.where(rmp_phases == run_n.log.rmp_phase)[0][0]

        total_fluxes_ratio[i, j, k, l, m] = run_n.log.total_stopped_power \
                                          / run_n.log.pinj * 100
    # Find where total_flux is still 0
    # counter = -1
    # for i in range(len(rmp_phases)):
    #     for j in range(len(rmp_responses)):
    #         for k in range(len(rmp_currents)):
    #             for l in range(len(bplasma_ns)):
    #                 for m in range(len(coil_sets)):
    #                     counter += 1
    #                     if total_fluxes_ratio[i, j, k, l, m] == 0:
    #                         print(i, j, k, l, m, counter)
    counter = 27
    for i in range(len(coil_sets)):
        for j in range(len(bplasma_ns)):
            for k in range(len(rmp_currents)):
                for l in range(len(rmp_responses)):
                    for m in range(8):
                        counter += 1
                        if total_fluxes_ratio[m, l, k, j, i] == 0:
                            print(m, l, k, j, i, counter)

    # RWM plot
    rwm_runs = []
    for run_i in runs:
        if run_i.log.bplasma:
            if run_i.log.bplasma_n == 1:
                rwm_runs.append(run_i)
    rwm_runs.sort(key=lambda x: x.log.rwm_bscale)
    for run_i in rwm_runs:
        print(run_i.log.rwm_bscale)

def plot_simulation_time_over_runs(runs_dir):
    """
    Plot the simulation time over all runs in runs_dir.
    """

    runs = run.create_runs_list(runs_dir)
    for run_i in runs:
        run_i.update_log()
    simulation_times = np.array([])
    for run_i in runs:
        simulation_times = np.append(simulation_times, run_i.log.simulation_time)

    fig, ax = plt.subplots()
    ax.plot(simulation_times / 60**2)
    ax.set_xlabel("Run")
    ax.set_ylabel("Simulation time [h]")
    ax.set_title("Simulation time over runs")

def plot_axisymmetric_constant_zeff_vs_non_constant():
    """
    Plot the axisymmetric constant zeff vs the non-constant zeff.
    """
    repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    vary_zeff_dir = os.path.join(repository_path, "output_data", "FEC_2024_missing_31_45_processed",
                                  "gpu-q-27")
    vary_zeff_tag = "06-03-2024_23-28-10.025"
    vary_zeff_run = run.Run(vary_zeff_dir, vary_zeff_tag)
    cnst_zeff_dir = os.path.join(repository_path, "output_data", "FEC_2024_missing_31_45_processed",
                                 "axisymmetric_fixed_processed")
    cnst_zeff_tag = "26-03-2024_11-47-55.495"
    run_vary_zeff = run.Run(vary_zeff_dir, vary_zeff_tag)
    run_cnst_zeff = run.Run(cnst_zeff_dir, cnst_zeff_tag)

    run_vary_zeff.update_log()
    run_cnst_zeff.update_log()

    print(run_vary_zeff.log.total_stopped_power / run_vary_zeff.log.pinj * 100)
    print(run_cnst_zeff.log.total_stopped_power / run_vary_zeff.log.pinj * 100)
    print(run_vary_zeff.log.total_stopped_power_error / run_vary_zeff.log.pinj * 100)
    print(run_cnst_zeff.log.total_stopped_power_error / run_vary_zeff.log.pinj * 100)

if __name__ == "__main__":
    repository_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    runs_directory = os.path.join(repository_path, "output_data",
                                  "FEC_2024_missing_31_45_processed")
    # plot_total_flux_over_runs(runs_directory)
    # plot_simulation_time_over_runs(runs_directory)
    plot_axisymmetric_constant_zeff_vs_non_constant()
    plt.show()
