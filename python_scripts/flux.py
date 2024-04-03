"""
The module contains routines associated with the energy flux on the PFCs.
It also contains information about the distribution of the energy flux on the PFCs.
"""

from typing import Optional
import numpy as np
from python_scripts import pkde

class Flux:
    """
    This class contains routines for calculating the energy flux on the PFCs
    """
    def __init__(self, num_grid_points = 1000):
        self.num_grid_points = num_grid_points
        self.s_theta: Optional[np.ndarray] = None
        self.s_phi: Optional[np.ndarray] = None
        self.total_energy: Optional[float] = None # Store the total energy flux in MW on the PFCs
        self.h_theta_1d: Optional[float] = None
        self.energy_arr_1d: Optional[np.ndarray] = None

    def calc_s_theta_s_phi(self, run):
        """
        Calculate the s_theta and s_phi attributes assuming that
        the wall attribute of the run object is populated.
        """
        self.s_phi = np.linspace(run.wall.s_phi_min,
                                 run.wall.s_phi_max,
                                 self.num_grid_points)
        self.s_theta = np.linspace(run.wall.s_theta_min,
                                   run.wall.s_theta_max,
                                   self.num_grid_points)

    def calc_total_energy(self, run):
        """
        Calculate the total energy flux on the Plasma-Facing Components (PFCs) in megawatts (MW).

        This method relies on specific attributes and states of the provided `run` object. 
        Before calling this method, ensure that the `run` object meets the following conditions:
        
        1. It has been initialized with the Flux class attribute.
        2. The `markers.stopped` attribute has been populated with relevant data.
        3. The `markers.all` attribute has been populated with relevant data.
        4. The `log` attribute has been populated and includes the `pinj` attribute, 
        which specifies the total power injected in MW.

        Args:
            run: A Run object
                An instance of the Run class that meets the above requirements.

        Returns:
            float: The total energy flux on the PFCs in MW. If the requirements are not 
                met, the behavior of this method is undefined and may result in errors.
        """
        energy0 = run.markers.all.vr0**2 + run.markers.all.vphi0**2 + run.markers.all.vz0**2
        denominator = np.sum(energy0 * run.markers.all.weight)
        stopped_energy = run.markers.stopped.vr**2 + run.markers.stopped.vphi**2 \
                       + run.markers.stopped.vz**2
        self.total_energy = run.log.pinj * np.sum(stopped_energy * run.markers.stopped.weight) \
                          / denominator

    def calc_energy_flux_1d(self, run,
                            h_theta_1d=None):
        """
        Calculate the energy flux on the PFCs using the distribution of the
        particle flux on the PFCs.

        Args:
            run: Run object
                Needs to have been initialized with the flux class attribute,
                also the wall and markers.stopped attributes need to
                have been populated. The flux object needs to have been
                initialized with the s_theta attribute, num_grid_points and the
                h_theta_1d attribute (if h_theta_1d is None).

        Returns:
            Flux object
                The Flux object with new attributes, namely the energy_fun_1d
                and energy_arr_1d which give the energy flux along the wall
                in MW/m^2.

        """

        if h_theta_1d is None:
            h_theta_1d = run.flux.h_theta_1d
        kde_fun = pkde.periodic_kde_1d(run.markers.stopped.s_theta,
                                        run.wall.s_theta_min,
                                        run.wall.s_theta_max,
                                        run.markers.stopped.weight,
                                        h_theta_1d,
                                        kernel='gaussian',
                                        num_grid_points=run.flux.num_grid_points)
        kde_arr = kde_fun(run.flux.s_theta)

        # Calculate bias corrected KDE
        asymptotic_bias = \
            pkde.calc_asymptotic_bias_1d(kde_fun, run.flux.s_theta, run.flux.h_theta_1d)
        kde_array = kde_array - asymptotic_bias

        # Calculate unstretched KDE
        kde_array *= run.wall.scale_factor_1d(run.flux.s_theta)

        # Calculate energy flux
        energy_flux = kde_array * run.flux.total

        run.flux.energy_kde_1d = run.flux.energy_fun_1d(run.flux.s_theta)
