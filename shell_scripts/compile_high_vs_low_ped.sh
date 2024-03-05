#!/bin/bash
absolute_value() {
    if [ $1 -lt 0 ]; then
      echo $((-$1))
    else
      echo $1
    fi
}

device="csd3"
tokamak="STEP"

niter=2
threadsPerBlock=256
blocksPerGrid=256
timax="1.0_gpu"
UNBOR=1000
dt0="1.0e-06_gpu"
BPLASMA=0
ncoil=16

SPR_strings_unique=("SPR-045-14" "SPR-054-clean" "SPR-054-dirty")
rcoils_unique=("0.0" "7.0" "8.0" "9.0")
BRIPPLES_unique=(0 1 1 1)
rcoils=()
BRIPPLES=()
SPR_strings=()
for SPR_string in "${SPR_strings_unique[@]}"; do
    for rcoil in "${rcoils_unique[@]}"; do
        SPR_strings+=("$SPR_string")
        rcoils+=("$rcoil")
        if [[ $rcoil == "0.0" ]]; then
            BRIPPLES+=(0)
        else
            BRIPPLES+=(1)
        fi
    done
done
echo ${rcoils[@]}
echo ${BRIPPLES[@]}
echo ${SPR_strings[@]}


run_name="high_vs_low_ped"

module purge
if [[ $device == "csd3" ]]; then
    module load rhel8/default-amp
    module load nvhpc/22.3/gcc-9.4.0-ywtqynx
    module load hdf5/1.10.7/openmpi-4.1.1/nvhpc-22.3-strpuv5
    export HDF5_DIR="/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/nvhpc-22.3/hdf5-1.10.7-strpuv55e7ggr5ilkjrvs2zt3jdztwpv"
    user_id="ir-prok1"
    root_dir="/home"
    cc="80"
    cuda="11.6"
elif [[ $device == "leonardo" ]]; then
    module load nvhpc/23.1
    user_id="aprokopy"
    root_dir="/leonardo/home/userexternal"
    cc="80"
    cuda="11.8"
    nohdf5=1
else
    echo "Invalid device."
    exit 1
fi

export LIBRARY_PATH=$LIBRARY_PATH:"$HDF5_DIR/lib"
export CFLAGS="-I$HDF5_DIR/include"
export FFLAGS="-I$HDF5_DIR/include"
export LDFLAGS="-L$HDF5_DIR/lib"

if [[ $tokamak == "ITER" ]]; then
    DTOKAMAK=1
elif [[ $tokamak == "STEP" ]]; then
    DTOKAMAK=10
else
    echo "Invalid tokamak."
    exit 1
fi

home_dir=$root_dir"/"$user_id
prec_file=$home_dir"/locust/prec_mod.f90"
input_dir=$home_dir"/STEP_input_files/high_vs_low_pedestal"

cd $home_dir"/locust"

cp -vf \
"makefile_template" \
$home_dir"/locust/makefile"
SRC="ccxx,cudaxx\.x"
DST="cc"$cc",cuda"$cuda
sed -i "s/$SRC/$DST/g" "$home_dir/locust/makefile"
diff "makefile_template" "$home_dir/locust/makefile"

rsync -avh \
$input_dir"/"$mesh_file \
$home_dir"/locust."$tokamak"/InputFiles/."

num_runs=${#rcoils[@]}
for ((n=0; n<num_runs; n++)); do

    FLAGS_BASE="-DCONLY -DPFCMOD -DTOKHEAD -DFSTATE -DLEIID=6 -DSTDOUT \
                -DSMALLEQ -DOPENTRACK -DOPENTERM -DPSIT=0.7 -DTOKAMAK="$DTOKAMAK" \
                -DNOTUNE -DUNBOR="$UNBOR" \
                -DTETALL -DSOLCOL \
                -DRFORCE -DBP -DTIMAX="$timax" -DWREAL -DWLIST"

    rcoil=${rcoils[$n]}
    SPR_string=${SPR_strings[$n]}
    BRIPPLE=${BRIPPLES[$n]}
    echo "ncoil="$ncoil
    echo "rcoil="$rcoil
    echo "SPR_string="$SPR_string
    echo "BRIPPLE="$BRIPPLE

    ptcle_file="2200k_uniform_"$SPR_string"_with_noise.dat"
    mesh_file="SPP-001-1.cdb.locust"
    eqdsk_file=$SPR_string".eqdsk"

    cp -f \
    $input_dir"/Profiles_"$SPR_string"/profile_ne.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$SPR_string"_ne.dat"
    cp -f \
    $input_dir"/Profiles_"$SPR_string"/profile_Ti.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$SPR_string"_Ti.dat"
    cp -f \
    $input_dir"/Profiles_"$SPR_string"/profile_Te.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$SPR_string"_Te.dat"
    cp -f \
    $input_dir"/"$eqdsk_file \
    $home_dir"/locust."$tokamak"/InputFiles/."
    cp -f \
    $input_dir"/"$ptcle_file \
    $home_dir"/locust."$tokamak"/InputFiles/."

    if [[ $BPLASMA == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DB3D -DB3D_EX"
    fi
    if [[ $nohdf5 == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DNOHDF5"
    fi
    if [[ $BRIPPLE == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DBRIPPLE"
    fi
    # echo $FLAGS_BASE

	cp -vf \
	$input_dir"/locust_scripts/base.f90" \
	$prec_file 

    SRC="TF_Ncoil = 00"
    DST="TF_Ncoil = "$ncoil
    sed -i "s/$SRC/$DST/g" $prec_file
    SRC="TF_Rcoil = 0.00_gpu"
    DST="TF_Rcoil = "$rcoil"_gpu"
    sed -i "s/$SRC/$DST/g" $prec_file
	SRC="'eqdsk_file.eqdsk' ! apkp"
	DST="'"$eqdsk_file"'"
	sed -i "s/$SRC/$DST/g" $prec_file
	SRC="'mesh_file.cdb.locust' ! apkp"
	DST="'"$mesh_file"'"
	sed -i "s/$SRC/$DST/g" $prec_file
	SRC="niter  = 0 ! apkp - Needs changing"
	DST="niter  = "$niter
	sed -i "s/$SRC/$DST/g" $prec_file
	SRC="'ptcle_file.dat' ! apkp - Needs changing"
	DST="'"$ptcle_file"'"
	sed -i "s/$SRC/$DST/g" $prec_file
	SRC="threadsPerBlock = 256"
	DST="threadsPerBlock = "$threadsPerBlock
	sed -i "s/$SRC/$DST/g" $prec_file
	SRC="blocksPerGrid   = 512"
	DST="blocksPerGrid   = "$blocksPerGrid
	sed -i "s/$SRC/$DST/g" $prec_file
	SRC="dt0    = 0.0e-00_gpu ! apkp - Needs changing"
	DST="dt0    = "$dt0
	sed -i "s/$SRC/$DST/g" $prec_file
	SRC="root = '\/home' ! apkp - Needs changing on Marconi"
	DST="root = '"$root_dir"'"
	sed -i "s|$SRC|$DST|g" $prec_file
    SRC="file_profile_ne = 'profile_ne_file' ! apkp"
    DST="file_profile_ne = 'profile_"$SPR_string"_ne.dat'"
	sed -i "s|$SRC|$DST|g" $prec_file
    SRC="file_profile_Te = 'profile_Te_file' ! apkp"
    DST="file_profile_Te = 'profile_"$SPR_string"_Te.dat'"
	sed -i "s|$SRC|$DST|g" $prec_file
    SRC="file_profile_Ti = 'profile_Ti_file' ! apkp"
    DST="file_profile_Ti = 'profile_"$SPR_string"_Ti.dat'"
	sed -i "s|$SRC|$DST|g" $prec_file
	
	make clean
	make FLAGS="$FLAGS_BASE" -j
	mv -f "locust" "locust_"$run_name"_"$n

done
