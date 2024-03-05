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
SPR_string="SPR-045-14"
input_dir="STEP_input_files/"$SPR_string
ptcle_file="2200k_uniform_"$SPR_string".dat"
mesh_file="SPP-001-1.cdb.locust"
# eqdsk_file=$SPR_string".geqdsk"
eqdsk_file=$SPR_string".eqdsk"

niter=2
threadsPerBlock=256
blocksPerGrid=256
timax="1.0_gpu"
UNBOR=1000
dt0="1.0e-06_gpu"

BPLASMA=0
BRIPPLE=1

ncoils_unique=(12 16 18)
rcoils_unique=("7.0" "7.25" "7.5" "7.75" "8.0" "8.25" "8.5" "8.75" "9.0")
ncoils=()
rcoils=()
for ncoil in "${ncoils_unique[@]}"; do
    for rcoil in "${rcoils_unique[@]}"; do
        ncoils+=("$ncoil")
        rcoils+=("$rcoil")
    done
done

run_name="vary_ripple"

module purge
if [[ $device == "csd3" ]]; then
    module load rhel8/default-amp
    module load nvhpc/22.3/gcc-9.4.0-ywtqynx
    module load hdf5/1.10.7/openmpi-4.1.1/nvhpc-22.3-strpuv5
    export HDF5_DIR="/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/nvhpc-22.3/hdf5-1.10.7-strpuv55e7ggr5ilkjrvs2zt3jdztwpv"
    user_id="ir-prok1"
    root_dir="/home"
	sed -i '/CUDALIB = cc*/c\CUDALIB = cc80,cuda11.6' "makefile_template"
elif [[ $device == "leonardo" ]]; then
    module load nvhpc/23.1
    user_id="aprokopy"
    root_dir="/leonardo/home/userexternal"
	sed -i '/CUDALIB = cc*/c\CUDALIB = cc80,cuda11.8' "makefile_template"
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
input_dir=$home_dir"/STEP_input_files/"$SPR_string
prec_file=$home_dir"/locust/prec_mod.f90"

cp -vf \
"makefile_template" \
$home_dir"/locust/makefile"
cp -f \
$input_dir"/Profiles/"* \
$home_dir"/locust."$tokamak"/InputFiles/."
cp -f \
$input_dir"/"$eqdsk_file \
$home_dir"/locust."$tokamak"/InputFiles/."
rsync -avh \
$input_dir"/"$mesh_file \
$home_dir"/locust."$tokamak"/InputFiles/."
cp -f \
$input_dir"/"$ptcle_file \
$home_dir"/locust."$tokamak"/InputFiles/."

cd $home_dir"/locust"

FLAGS_BASE="-DCONLY -DPFCMOD -DTOKHEAD -DFSTATE -DLEIID=6 -DSTDOUT \
            -DSMALLEQ -DOPENTRACK -DOPENTERM -DPSIT=0.7 -DTOKAMAK="$DTOKAMAK" \
            -DNOTUNE -DUNBOR="$UNBOR" \
            -DTETALL -DSOLCOL \
            -DRFORCE -DBP -DTIMAX="$timax" -DWREAL -DWLIST"
if [[ $BPLASMA == 1 ]]; then
    FLAGS_BASE=$FLAGS_BASE" -DB3D -DB3D_EX"
fi
if [[ $nohdf5 == 1 ]]; then
    FLAGS_BASE=$FLAGS_BASE" -DNOHDF5"
fi
if [[ $BRIPPLE == 1 ]]; then
    FLAGS_BASE=$FLAGS_BASE" -DBRIPPLE"
fi
echo $FLAGS_BASE

length_of_tf_coil_sets=${#ncoils[@]}
for ((n=0; n<length_of_tf_coil_sets; n++)); do
	phase=${phases[$n]}
	echo "n="$n
    rcoil=${rcoils[$n]}
    ncoil=${ncoils[$n]}
    echo "ncoil="$ncoil
    echo "rcoil="$rcoil

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
	
	make clean
	make FLAGS="$FLAGS_BASE" -j
	mv -f "locust" "locust_"$run_name"_"$n

done
