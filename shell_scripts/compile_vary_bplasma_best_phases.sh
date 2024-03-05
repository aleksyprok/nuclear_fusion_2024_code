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

niter=3
threadsPerBlock=128
blocksPerGrid=128
timax="1.0_gpu"
UNBOR=1000
dt0="1.0e-06_gpu"

# BPLASMA variables:
BPLASMA=1
ripple=0
sideband=0
toroidal_modes_unique=(-2 -3 -4)
coil_sets_unique=("efcc" "rwm")
responses_unique=(0 1)
current_multipliers_unique=(1 2)
coil_sets=()
toroidal_modes=()
currents=()
responses=()
phases=()
for coil_set in "${coil_sets_unique[@]}"; do
    for toroidal_mode in "${toroidal_modes_unique[@]}"; do
        for current_multiplier in "${current_multipliers_unique[@]}"; do
            if [[ $coil_set == "efcc" ]]; then
                if [[ $(($toroidal_mode * $toroidal_mode)) -eq 4 ]]; then
                    current=$(($current_multiplier * 50))
                    phase="061.0"
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 9 ]]; then
                    current=$(($current_multiplier * 90))
                    phase="020.0"
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 16 ]]; then
                    current=$(($current_multiplier * 150))
                    phase="321.0"
                else
                    echo "Invalid toroidal number"
                    exit 1
                fi
            elif [[ $coil_set == "rwm" ]]; then
                if [[ $(($toroidal_mode * $toroidal_mode)) -eq 4 ]]; then
                    current=$(($current_multiplier * 30))
                    phase="265.0"
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 9 ]]; then
                    current=$(($current_multiplier * 50))
                    phase="173.0"
                elif [[ $(($toroidal_mode * $toroidal_mode)) -eq 16 ]]; then
                    current=$(($current_multiplier * 80))
                    phase="067.0"
                else
                    echo "Invalid toroidal number"
                    exit 1
                fi
            else
                echo "Invalid coil set"
                exit 1
            fi
			for response in "${responses_unique[@]}"; do
                coil_sets+=("$coil_set")
                toroidal_modes+=("$toroidal_mode")
                responses+=("$response")
                currents+=("$current")
                phases+=("$phase")
			done
        done
    done
done
length_of_coil_sets=${#coil_sets[@]}
echo "Length of coil_sets: $length_of_coil_sets"
run_name="vary_dave_bplasma_best_phases"

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
echo $FLAGS_BASE

for ((n=0; n<length_of_coil_sets; n++)); do
	coil_set=${coil_sets[$n]}
	toroidal_mode=${toroidal_modes[$n]}
	toroidal_mode_abs=$(absolute_value $toroidal_mode)
	echo "toroidal_mode_abs="$toroidal_mode_abs
	current=${currents[$n]}
	response=${responses[$n]}
	phase=${phases[$n]}
	echo "n="$n

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
	
	nnums=($toroidal_mode)
	if [[ $ripple == 1 ]]; then
	    nnums+=($Ncoil)
	fi
	if [[ $sideband == 1 ]]; then
	    nnums+=($sideband_mode)
	fi
	SRC="nnum   = \[16, -3\] ! apkp - Needs changing"
	DST="nnum   = ["
	if [[ $BPLASMA == 1 ]]; then
		DST+=$toroidal_mode
	fi
	if [[ $ripple == 1 ]]; then
		DST+=", "$Ncoil
	fi
	if [[ $sideband == 1 ]]; then
		DST+=", "$sideband_mode
	fi
	DST+="]"
	sed -i "s/$SRC/$DST/g" $prec_file
	echo "DST="$DST
	
	nmde=$((BPLASMA+ripple+sideband))
	SRC="nmde   = 2 ! apkp - Needs changing"
	DST="nmde   = "$nmde
	sed -i "s/$SRC/$DST/g" $prec_file
	
	SRC="phase  = \[0.0e0_gpu, 0.0e0_gpu\] ! apkp - Needs changing"
	DST="phase  = ["
	if [[ $nmde > 0 ]]; then
		DST+="0.0e0_gpu"
	fi
	for (( i=1; i<nmde; i++ )); do
		DST+=", 0.0e0_gpu"
	done 		
	DST+="]"
	sed -i "s/$SRC/$DST/g" $prec_file
	echo "DST="$DST			
	BPLASMA_parts=(
		"BPLASMA_"$coil_set
		"_response="$response
	    "_current="$(printf "%03d" $current)
		"_100x200"
		"_phase="$(printf "%05.1f" $phase)
		"_n")
	bplasma_file=$(printf "%s" "${BPLASMA_parts[@]}")
	SRC="'bplasma_file' ! apkp"
	DST="'"$bplasma_file"'"
	sed -i "s|$SRC|$DST|g" $prec_file

	BPLASMA_directory="DAVE_v2/"$coil_set"_n"$toroidal_mode
	cp -vf \
	$input_dir"/BPLASMA/"$BPLASMA_directory"/"$bplasma_file$toroidal_mode_abs \
	$home_dir"/locust."$tokamak"/InputFiles/."
	
	make clean
	make FLAGS="$FLAGS_BASE" -j
	mv -f "locust" "locust_"$run_name"_"$n

done
