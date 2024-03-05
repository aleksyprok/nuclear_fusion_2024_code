#!/bin/bash

# This shell script will compile and produce a range of binary files for LOCUST
# to execute. We compile every simulation used for the IAEA FEC 2024 paper.
# This includes simulations where we have an axisymmetric simulation,
# where we vary the ripple field, where we vary the RMP field and where
# we vary the amplitude of a RWM.

absolute_value() {
    if [ $1 -lt 0 ]; then
      echo $((-$1))
    else
      echo $1
    fi
}

device="csd3"
tokamak="STEP"
SPR_string="SPR-045-16"
input_dir="STEP_input_files/"$SPR_string
ptcle_file=$SPR_string"_markers_1000000.dat"
mesh_file="SPP-001-1.cdb.locust"
eqdsk_file=$SPR_string".eqdsk"

niter=2
threadsPerBlock=256
blocksPerGrid=256
timax="1.0_gpu"
UNBOR=1000
dt0="1.0e-06_gpu"

# Make a load of lists to store the values of the parameters we are varying
rcoils=()
ncoils=()
bripples=()
toroidal_modes=()
coil_sets=()
currents=()
responses=()
gain_values=()
bscales=()
bplasmas=()

# Axisymmetric simulations
rcoil=0
ncoil=0
bripple=0
toroidal_mode=0
coil_set="none"
current=0
response=0
phase=0
gain_value=0
bscale=0
bplasma=0
rwm_control=0
rmp=0
rcoils+=("$rcoil")
ncoils+=("$ncoil")
bripples+=("$bripple")
toroidal_modes+=("$toroidal_mode")
coil_sets+=("$coil_set")
currents+=("$current")
responses+=("$response")
phases+=("$phase")
gain_values+=("$gain_value")
bscales+=("$bscale")
bplasmas+=("$bplasma")
rwm_contrls+=("$rwm_control")
rmps+=("$rmp")

# Ripple simulations
bripple=1
toroidal_mode=0
coil_set="none"
current=0
response=0
phase=0
gain_value=0
bscale=0
bplasma=0
rwm_control=0
rmp=0
ncoils_unique=(12 16 18)
rcoils_unique=("7.0" "7.25" "7.5" "7.75" "8.0" "8.25" "8.5" "8.75" "9.0")
for ncoil in "${ncoils_unique[@]}"; do
    for rcoil in "${rcoils_unique[@]}"; do
        ncoils+=("$ncoil")
        rcoils+=("$rcoil")
        bripples+=("$bripple")
        toroidal_modes+=("$toroidal_mode")
        coil_sets+=("$coil_set")
        currents+=("$current")
        responses+=("$response")
        phases+=("$phase")
        gain_values+=("$gain_value")
        bscales+=("$bscale")
        bplasmas+=("$bplasma")
        rwm_contrls+=("$rwm_control")
        rmps+=("$rmp")
    done
done

# RMP simulations

# RMP special phases simulations

# RWM control simulations