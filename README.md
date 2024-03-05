# nuclear_fusion_2024_code

**Warning:** This repository is >25GB so make sure you have enough space before cloning. Ensure that [Git LFS](https://git-lfs.com/) is installed to download the files in the `input_data` directory.

The codebase associated with our research, as detailed in the paper available for review and use at the following GitHub repository: <https://github.com/aleksyprok/Prokopyszyn_et_al_Nuclear_Fusion_IAEA_2023>

## Repository Structure

This repository is organized as follows:

`/input_data` - Input datasets required by the scripts are stored here.

`/python_scripts` - Core Python scripts for processing and analysis:

- `__init__.py` - Indicates that this directory is treated as a Python package.
  - `prepare_profiles.py` - The main script that processes profile data for simulations.
  - `my_gfile_reader.py` - Module for reading in eqdsk files (not written by me).

`/tests` - Contains automated test cases for the Python scripts:

- `test_prepare_profiles.py` - Unit tests for `prepare_profiles.py`.

`/.gitattributes` - Specifies attributes for pathnames and handling of line endings across different environments. This is also where we specify the large files for the Git Large File Storage (LFS).

`/.gitignore` - Lists files and directories that are to be ignored by version control.

`/README.md` - Provides an overview and documentation for the project.

 `/setup.py` - Used for packaging and distributing the project, defining the package name, version, and included packages.

## How we produced the data for the paper

### Step 1: Produce density/temperature profiles and alpha particle initial position.

The first thing we did was to produce the profile data. This was done using the `prepare_profiles.py` script. This script takes the raw data from the input directory and processes it into a format that can be used for simulations.
