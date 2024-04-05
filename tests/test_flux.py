"""
Test the flux module.
"""
import os
import matplotlib.pyplot as plt
import numpy as np
from python_scripts import flux, run

def test_calc_total_energy():
    """
    Test the calc_total_energy function.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points=1000)

    assert np.isclose(test_run.flux.total_energy,
                      test_run.log.total_stopped_power,
                      atol=1e-2)

def test_calc_energy_flux_1d():
    """
    Test the calc_energy_flux_1d function.
    """

    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    output_dir = os.path.join(repo_path, 'tests', 'output_plots')
    os.makedirs(output_dir, exist_ok=True)
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    test_run.init_log()
    wall_path = os.path.join(repo_path, 'input_data', 'SPP-001_wall.dat')
    test_run.init_wall(wall_path)
    gfile_path = os.path.join(repo_path, 'input_data', 'SPR-045-16.eqdsk')
    test_run.init_gfile(gfile_path)
    test_run.init_markers()
    test_run.init_flux(num_grid_points=1000)
    flux.calc_energy_flux_1d(test_run, h_theta_1d=0.1)
    fig, ax = plt.subplots()
    ax.plot(test_run.flux.s_theta, test_run.flux.energy_1d)
    ax.set_xlabel('s_theta [m]')
    ax.set_ylabel('Energy flux [W/m^2]')
    fig.savefig(os.path.join(output_dir, 'axisymmetric_energy_flux_1d.png'))
    assert test_run.flux.h_theta_1d == 0.1
