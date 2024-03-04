"""
Test module for the run.py module.
"""

import os
from python_scripts import run

def test_run_init():
    """
    Test the Run class initialisation.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    dir_path = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric")
    tag = '13-12-2023_16-51-52.811'
    test_run = run.Run(dir_path, tag)
    assert test_run.dir_path == dir_path
    assert test_run.tag == tag
    assert test_run.log is None
    assert test_run.fstate is None
    assert test_run.log_path == dir_path + f'/LOG_{tag}.out'
    assert test_run.fstate_path == dir_path + f'/FINAL_STATE_{tag}.dat'
