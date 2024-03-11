#!/bin/bash

device="aws_v100"
tokamak="STEP"

run_name="FEC_2024"

if [[ $device == "aws_v100" ]]; then
    root_dir="/home"
    user_id="qt4627"
    cc="70"
    cuda="11.0"
    ngpu=1
    nohdf5=1
else
    echo "Invalid device."
    exit 1
fi

export OMP_NUM_THREADS=$ngpu

echo "OMP_NUM_THREADS="\$OMP_NUM_THREADS

ulimit -s 2000000
export OMP_STACKSIZE=102400
export CUDA_CACHE_DISABLE=1

echo "OMP_NUM_THREADS="\$OMP_NUM_THREADS

# Clear CacheFiles
echo $HOSTNAME
rm -vf $root_dir"/"$user_id"/locust."$tokamak"/CacheFiles/"\$HOSTNAME"/"*
$root_dir"/"$user_id"/locust/locust_"$run_name"_"0