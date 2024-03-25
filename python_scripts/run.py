"""
This module contains routines for updating a Run object which stores information
about a LOCUST run.

The run object contain three nested classes, log, stopped and grid.
The log class contains information obtained by reading the LOG*.out file.
The stopped class contains information obtained by reading the FINAL_STATE*.dat file
and only stores information about particles that have hit the wall.
The wall class contains information obtained by reading the wall 2D wall files, e.g.
input_data/SPP-001_wall.dat.
Note that the log class and grid classes are standalone classes, but many of the methods in the
stopped module require the log and grid classes.
"""
import glob
import os
from typing import Optional
from python_scripts import fstate, log

class Run:
    """
    This class contains routines for updating a Run object which stores information
    about a LOCUST run.
    """
    def __init__(self, dir_path, tag):
        self.log: Optional[log.Log] = None
        self.fstate: Optional[fstate.Fstate] = None
        self.dir_path = dir_path
        self.tag = tag
        self.log_path = self.dir_path + f'/LOG_{self.tag}.out'
        self.fstate_path = self.dir_path + f'/FINAL_STATE_{self.tag}.dat'
        self.einj = None

    def update_log(self):
        """
        Update the Run object with information from the LOCUST .log file.
        """
        self.log = log.Log(self.log_path)

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
