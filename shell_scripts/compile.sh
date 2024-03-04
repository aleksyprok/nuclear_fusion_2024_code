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
BPLASMA=0
ncoil=16
