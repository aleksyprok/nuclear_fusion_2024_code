"""
The module contains routines associated with the energy flux on the PFCs.
It also contains information about the distribution of the energy flux on the PFCs.
"""

from python_scripts import pkde

class Flux:
    """
    This class contains routines for calculating the energy flux on the PFCs
    """
    def __init__(self, num_grid_points = 1000):
        self.num_grid_points = num_grid_points
        self.s_phi = None
        self.s_theta = None
        self.h_theta_1d = None
        self.energy_fun_1d = None
        self.energy_arr_1d = None

def calc_energy_flux_1d(run,
                        h_theta_1d=None):
    """
    Calculate the energy flux on the PFCs using the distribution of the
    particle flux on the PFCs.

    Args:
        run: Run object
            Needs to have been initialized with the flux class attribute,
            also the wall and markers.stopped attributes need to
            have been populated. The flux object needs to have been
            initialized with the s_theta attribute and the
            h_theta_1d attribute.

    Returns:
        Flux object
            The Flux object with new attributes, namely the

    """

    if h_theta_1d is None:
        h_theta_1d = run.flux.h_theta_1d
    run.flux.energy_fun_1d = pkde.periodic_kde_1d(run.markers.stopped.s_theta,
                                                  run.flux.s_theta[0],
                                                  run.flux.s_theta[1],
                                                  run.markers.stopped.weights,
                                                  h_theta_1d,
                                                  kernel='gaussian',
                                                  num_grid_points=run.flux.num_grid_points)
    run.flux.energy_arr_1d = run.flux.energy_fun_1d(run.flux.s_theta)
