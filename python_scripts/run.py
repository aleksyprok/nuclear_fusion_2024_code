"""
This module contains routines for updating a Run object which stores information
about a LOCUST run.

The run object has four classes as attributes:
- log: contains information obtained by reading the LOG*.out file
- wall: contains routines from converting from R, Phi, Z coordinates on the wall to
        wall coordinates, s_phi, s_theta.
- markers: contains routines related to the particle groups from reading the FINAL_STATE*.dat file.
- flux: contains routines related to the energy flux on the PFCs.
Note that the log wall and markers classes are standalone classes, but many of the methods in the
flux class require the log, wall and markers classes.
Note that to save memory we will often set the markers attribute to None.
"""
import glob
import os
from typing import Optional
from python_scripts import log, markers, my_gfile_reader, wall

class Run:
    """
    This class contains routines for updating a Run object which stores information
    about a LOCUST run.
    """
    def __init__(self, dir_path, tag):
        self.log: Optional[log.Log] = None
        self.wall: Optional[wall.Wall] = None
        self.markers: Optional[markers.Markers] = None
        self.gfile: Optional[my_gfile_reader.getGfile] = None
        self.dir_path = dir_path
        self.tag = tag
        self.log_path = self.dir_path + f'/LOG_{self.tag}.out'
        self.fstate_path = self.dir_path + f'/FINAL_STATE_{self.tag}.dat'

    def update_log(self):
        """
        Update the Run object with information from the LOCUST .log file.
        """
        self.log = log.Log(self.log_path)

    def update_wall(self, wall_path):
        """
        Update the Run object with information from the LOCUST .log file.
        """
        self.wall = wall.Wall(wall_path)

    def update_markers(self):
        """
        Update the Run object with information from the LOCUST FINAL_STATE*.dat file.
        """
        self.markers = markers.Markers(self.fstate_path)

    def update_gfile(self, gfile_path):
        """
        Update the Run object with information from the eqdsk file.
        """
        self.gfile = my_gfile_reader.getGfile(gfile_path)

def create_runs_list(dir_path):
    """
    Create a list of Run objects.

    Args:
        dir_path: str
            The path to the directory containing the runs.
    
    Returns:
        list
            A list of Run objects.
    """
    runs = []
    for file in glob.iglob(f'{dir_path}/**/FINAL_STATE*.dat', recursive=True):
        tag = os.path.splitext(os.path.basename(file))[0].split('FINAL_STATE_')[-1]
        runs.append(Run(os.path.dirname(file), tag))
    return runs
