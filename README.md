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

 ### Installation instructions:

  - Ensure `git` and `git-lfs` are installed on your machine.
  - Clone the repository. I reccommend checking that the the large files have downlaoded correctly by looking at the files in `input_data`.
  - To run the python scripts you need to have python3.9+ installed.
  - Setup a virtual environment with e.g. "python3 -m venv venv".
  - Activate the environment with e.g. "source venv/bin/activate"
  - Install python packages with "pip install -r requirements.txt"
  - Install the repo itself with "pip install -e .".

## How we produced the data for the paper

### Step 1: Produce density/temperature profiles and alpha particle initial position.

The first thing we did was to produce the profile data. This was done using the `prepare_profiles.py` script. This script takes the raw data from the input directory and processes it into a format that can be used for simulations.

### Step 2: ssh into GPU cluster

For this work I used the CSD3 cluster:
https://docs.hpc.cam.ac.uk/hpc/index.html

### Step 3: Clone locust code

For the shell scripts to work the `locust` code must be cloned into the home space and you need to run
``git checkout 1f28dcdc16a89e976832beaffd7e207b4d978c5b``
to get the same version of the code. The LOCUST repo is stored on GitLab here:
[https://git.ccfe.ac.uk/fast-particles/locust](https://git.ccfe.ac.uk/fast-particles/locust)

### Step 4: Clone this repository

Clone this repo to a directory of your choosing on CSD3. Make sure you have `git-lfs`, I load these modules in my `.bashrc` to ensure that I have `git-lfs`. Note that this also sets up `slurm` which you iwll need to run jobs:
```
module purge
module load gcc/5
module load rhel8/default-amp
module load slurm-current-gcc-5.4.0-6idu76o
module load git-lfs-2.3.0-gcc-5.4.0-oktvmkw
```

### Step 5: Create output directories

Create the file structure for LOCUST to output its files to. Your file strcuture needs to be:
- `~/locust.STEP/Outputfiles`
- `~/locust.STEP/InputFiles`
- `~/locust.STEP/CacheFiles`
I reccommend Using symbolic links to ensure you files are saved in e.g. `/rds/project/iris_vol2/rds-ukaea-ap001/` as you only have a limited amount of space in your home drive.

### Step 6: Compile the code.

Navigate to `nuclear_fusion_2024_code/shell_scripts/compile.sh`. Edit the the file and ensure `device="csd3"` (assuming you are running the code on CSD3). Then create a session onthe GPU nodes with e.g.:
```
sintr -t 36:0:0 -N1 -A UKAEA-AP001-GPU -p ampere --qos=INTR --gres=gpu:1
```
Then run `./compile.sh` to compile the code.

### Step 7: Run the code.

After the code has finished compiling. Edit `nuclear_fusion_2024_code/shell_scripts/run.sh` and ensure the setting are ok. You may need to edit the lines that read
```
#SBATCH --ntasks=248
#SBATCH --array=0-247
```
to 
```
#SBATCH --ntasks=124
#SBATCH --array=0-123
```
as sometime CSD3 won't let you execute so many jobs as a single array job. So you need to submit two seperate array jobs.
