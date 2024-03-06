#!/bin/bash
absolute_value() {
    if [ $1 -lt 0 ]; then
      echo $((-$1))
    else
      echo $1
    fi
}

device="leonardo"
tokamak="STEP"
SPR_string="SPR-045-14"
input_dir="STEP_input_files/"$SPR_string
ptcle_file="2200k_uniform_"$SPR_string"_with_noise.dat"
mesh_file="SPP-001-1.cdb.locust"
# eqdsk_file=$SPR_string".geqdsk"
eqdsk_file=$SPR_string".eqdsk"

niter=1
threadsPerBlock=256
blocksPerGrid=256
timax="1.0_gpu"
UNBOR=1000
dt0="1.0e-06_gpu"

# BPLASMA variables:
BPLASMA=1
toroidal_mode=-1
toroidal_mode_abs=$(absolute_value $toroidal_mode)
bscales_unique=(1 10)
gain_values_unique=("0.0" "1.0" "2.0" "2.6")
bscales=()
gain_values=()
for gain_value in "${gain_values_unique[@]}"; do
    for bscale in "${bscales_unique[@]}"; do
        gain_values+=($gain_value)
        bscales+=($bscale)
    done
done
echo ${bscales[@]}
echo ${gain_values[@]}
run_name="vary_rwm_tesla"
num_runs=${#bscales[@]}

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
input_dir=$home_dir"/STEP_input_files/"$SPR_string
prec_file=$home_dir"/locust/prec_mod.f90"

cp -vf \
"makefile_template" \
$home_dir"/locust/makefile"
SRC="ccxx,cudaxx\.x"
DST="cc"$cc",cuda"$cuda
sed -i "s/$SRC/$DST/g" "$home_dir/locust/makefile"
diff "makefile_template" "$home_dir/locust/makefile"

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

for ((n=0; n<num_runs; n++)); do
# for ((n=0; n<2; n++)); do
	coil_set=${coil_sets[$n]}
    bscale=${bscales[$n]}
    gain_value=${gain_values[$n]}
	echo "n="$n
    echo "gain_value="$gain_value
    echo "bscale="$bscale

	cp -vf \
	$input_dir"/locust_scripts/base.f90" \
	$prec_file 
	
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
	
	# BPLASMA code
	
	SRC="nnum   = \[16, -3\] ! apkp - Needs changing"
	DST="nnum   = ["$toroidal_mode"]"
	sed -i "s/$SRC/$DST/g" $prec_file
	
	SRC="nmde   = 2 ! apkp - Needs changing"
	DST="nmde   = 1"
	sed -i "s/$SRC/$DST/g" $prec_file
	
	SRC="phase  = \[0.0e0_gpu, 0.0e0_gpu\] ! apkp - Needs changing"
	DST="phase  = [0.0e0_gpu]"
	sed -i "s/$SRC/$DST/g" $prec_file

	BPLASMA_parts=(
		"BPLASMA_cylindrical_tesla_G="$gain_value
		"_bscale="$bscale
		"_n")
	bplasma_file=$(printf "%s" "${BPLASMA_parts[@]}")
	SRC="'bplasma_file' ! apkp"
	DST="'"$bplasma_file"'"
	sed -i "s|$SRC|$DST|g" $prec_file

	BPLASMA_directory="RWM_control/tesla"
	cp -vf \
	$input_dir"/BPLASMA/"$BPLASMA_directory"/"$bplasma_file$toroidal_mode_abs \
	$home_dir"/locust."$tokamak"/InputFiles/."
    
    if [[ $n == 0 ]]; then
        diff $prec_file $input_dir"/locust_scripts/base.f90"
    fi
	
	make clean
	make FLAGS="$FLAGS_BASE" -j
	mv -f "locust" "locust_"$run_name"_"$n

done
