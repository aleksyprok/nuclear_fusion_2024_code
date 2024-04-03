"""
This file contains routines related to the Markers class which stored inforamtion
from reading a LOCUST FINAL_STATE*.dat file.
"""
import os
import time
import numpy as np
from python_scripts import wall

class ParticleGroup:
    """
    This class is used to store the information of a group of particles.
    """
    def __init__(self):
        # Initialize empty arrays for each attribute
        self.r = np.array([])
        self.phi = np.array([])
        self.z = np.array([])
        self.vr = np.array([])
        self.vphi = np.array([])
        self.vz = np.array([])
        self.t = np.array([])
        self.s = np.array([])
        self.particle_id = np.array([])
        self.r0 = np.array([])
        self.phi0 = np.array([])
        self.z0 = np.array([])
        self.vr0 = np.array([])
        self.vphi0 = np.array([])
        self.vz0 = np.array([])
        self.weight = np.array([])
        self.s_phi = np.array([])
        self.s_theta = np.array([])

    def add_particles(self, data):
        """
        Add particles to the group.

        Args:
            data (np.ndarray): The data to add to the group. Each row should have 16 columns
        """
        self.r = np.concatenate((self.r, data[:, 0]))
        self.phi = np.concatenate((self.phi, data[:, 1]))
        self.z = np.concatenate((self.z, data[:, 2]))
        self.vr = np.concatenate((self.vr, data[:, 3]))
        self.vphi = np.concatenate((self.vphi, data[:, 4]))
        self.vz = np.concatenate((self.vz, data[:, 5]))
        self.t = np.concatenate((self.t, data[:, 6]))
        self.s = np.concatenate((self.s, data[:, 7]))
        self.particle_id = np.concatenate((self.particle_id, data[:, 8]))
        self.r0 = np.concatenate((self.r0, data[:, 9]))
        self.phi0 = np.concatenate((self.phi0, data[:, 10]))
        self.z0 = np.concatenate((self.z0, data[:, 11]))
        self.vr0 = np.concatenate((self.vr0, data[:, 12]))
        self.vphi0 = np.concatenate((self.vphi0, data[:, 13]))
        self.vz0 = np.concatenate((self.vz0, data[:, 14]))
        self.weight = np.concatenate((self.weight, data[:, 15]))

class Markers:
    """
    This class contains routines for updating a Run object with information from a
    LOCUST FINAL_STATE*.dat file.
    """
    def __init__(self, fstate_path):
        self.fstate_path = fstate_path
        self.moving = ParticleGroup()
        self.thermal = ParticleGroup()
        self.stopped = ParticleGroup()
        self.unresolved = ParticleGroup()
        print(f"Processing file: {os.path.basename(self.fstate_path)}")
        start_time = time.time()
        self._process_file()
        end_time = time.time()
        print(f"Processed file: {os.path.basename(self.fstate_path)}")
        print(f"Time taken: {end_time - start_time:.2f} seconds")

    def _process_file(self):
        data = np.loadtxt(self.fstate_path)
        # Assuming the data format is consistent and each row has 16 columns
        for status_code in [-14, -9, -5, -3]:
            filtered_data = data[data[:, 7] == status_code]  # Assuming status is in the 8th column
            self._add_to_group(filtered_data, status_code)

    def _add_to_group(self, data, status_code):
        if status_code == -14:
            self.moving.add_particles(data)
        elif status_code == -9:
            self.thermal.add_particles(data)
        elif status_code == -5:
            self.stopped.add_particles(data)
        elif status_code == -3:
            self.unresolved.add_particles(data)

    # def calculate_s_phi_s_theta(self, r_wall, z_wall):
def get_s_phi_s_theta_from_r_z_phi(run):
    """
    Calculates and returns the s_phi and s_theta coordinates based on the r, z, and phi coordinates
    of the markers and the wall. This function relies on specific attributes of the `run` object.

    The `run` object must have the following attributes:
    - run.markers.stopped.r: np.ndarray
        The r coordinates of the markers.
    - run.markers.stopped.phi: np.ndarray
        The phi coordinates of the markers.
    - run.markers.stopped.z: np.ndarray
        The z coordinates of the markers.
    - run.wall.r: np.ndarray
        The r coordinates of the wall.
    - run.wall.z: np.ndarray
        The z coordinates of the wall.

    This function computes and assigns new values to the following attributes of `run`:
    - run.markers.stopped.s_phi: np.ndarray
        The computed s_phi coordinates.
    - run.markers.stopped.s_theta: np.ndarray
        The computed s_theta coordinates.

    Args:
        run: A custom object
            An object representing the run data, containing markers and wall information.

    Returns:
        None
    """
    run.markers.stopped.s_phi = wall.get_s_phi_from_phi(run.markers.stopped.phi,
                                                        run.wall.r, run.wall.z)
    run.markers.stopped.s_theta = wall.get_s_theta_from_rz(run.markers.stopped.r,
                                                           run.markers.stopped.z,
                                                           run.wall.r,
                                                           run.wall.z)
