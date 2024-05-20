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
input_dir="$(pwd)/../input_data"
home_dir=$HOME
prec_file=$home_dir"/locust/prec_mod.f90"
mesh_file="SPP-001-1.cdb.locust"
run_name="francis_q_profile"
spr_strings_unique=("SPR-045-14" "SPR-045-16" "SPR046_16b")

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
phases=()
gain_values=()
bscales=()
bplasmas=()
rwm_controls=()
rmps=()
spr_strings=()

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
spr_string="none"
for spr_string in "${spr_strings_unique[@]}"; do
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
    rwm_controls+=("$rwm_control")
    rmps+=("$rmp")
    spr_strings+=("$spr_string")
done

# Ripple simulations
rcoil="7.25"
ncoil=16
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
spr_string="none"
for spr_string in "${spr_strings_unique[@]}"; do
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
    rwm_controls+=("$rwm_control")
    rmps+=("$rmp")
    spr_strings+=("$spr_string")
done

# RMP simulations
rcoil=0
ncoil=0
bripple=0
toroidal_mode=-3
coil_set="efcc"
current=90
response=1
phase="020.00"
gain_value=0
bscale=0
bplasma=1
rwm_control=0
rmp=1
spr_string="none"
for spr_string in "${spr_strings_unique[@]}"; do
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
    rwm_controls+=("$rwm_control")
    rmps+=("$rmp")
    spr_strings+=("$spr_string")
done
spr_string="none"

if [[ $device == "csd3" ]]; then
    module purge
    module load rhel8/default-amp
    module load nvhpc/22.3/gcc-9.4.0-ywtqynx
    module load hdf5/1.10.7/openmpi-4.1.1/nvhpc-22.3-strpuv5
    export HDF5_DIR="/usr/local/software/spack/spack-rhel8-20210927/opt/spack/linux-centos8-zen2/nvhpc-22.3/hdf5-1.10.7-strpuv55e7ggr5ilkjrvs2zt3jdztwpv"
    export LIBRARY_PATH=$LIBRARY_PATH:"$HDF5_DIR/lib"
    export CFLAGS="-I$HDF5_DIR/include"
    export FFLAGS="-I$HDF5_DIR/include"
    export LDFLAGS="-L$HDF5_DIR/lib"
    root_dir="/home"
    cc="80"
    cuda="11.6"
    nohdf5=0
elif [[ $device == "leonardo" ]]; then
    module purge
    module load nvhpc/23.1
    root_dir="/leonardo/home/userexternal"
    cc="80"
    cuda="11.8"
    nohdf5=1
elif [[ $device == "aws_v100" ]]; then
    export PATH=$PATH:/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/compilers/bin
    root_dir="/home"
    cc="70"
    cuda="11.0"
    nohdf5=1
else
    echo "Invalid device."
    exit 1
fi

cp -vf \
$input_dir/"makefile_template" \
$home_dir"/locust/makefile"
SRC="ccxx,cudaxx\.x"
DST="cc"$cc",cuda"$cuda
echo "nohdf5 = "$nohdf5
sed -i "s/$SRC/$DST/g" "$home_dir/locust/makefile"
if [[ $nohdf5 == 1 ]]; then
    SRC="-lstdc++ -lhdf5_fortran -lhdf5"
    DST="-lstdc++"
    sed -i "s/$SRC/$DST/g" "$home_dir/locust/makefile"
fi
diff $input_dir/"makefile_template" "$home_dir/locust/makefile"

rsync -avh \
$input_dir"/"$mesh_file \
$home_dir"/locust."$tokamak"/InputFiles/."

cd $home_dir"/locust"

num_runs=${#rcoils[@]}
for ((n=0; n<num_runs; n++)); do

    rcoil=${rcoils[$n]}
    ncoil=${ncoils[$n]}
    bripple=${bripples[$n]}
    toroidal_mode=${toroidal_modes[$n]}
    coil_set=${coil_sets[$n]}
    current=${currents[$n]}
    response=${responses[$n]}
    phase=${phases[$n]}
    gain_value=${gain_values[$n]}
    bscale=${bscales[$n]}
    bplasma=${bplasmas[$n]}
    rwm_control=${rwm_controls[$n]}
    rmp=${rmps[$n]}
    spr_string=${spr_strings[$n]}
    echo "ncoil="$ncoil
    echo "rcoil="$rcoil
    echo "bripple="$bripple
    echo "toroidal_mode="$toroidal_mode
    echo "coil_set="$coil_set
    echo "current="$current
    echo "response="$response
    echo "phase="$phase
    echo "gain_value="$gain_value
    echo "bscale="$bscale
    echo "bplasma="$bplasma
    echo "rwm_control="$rwm_control
    echo "rmp="$rmp
    echo "spr_string="$spr_string

    toroidal_mode_abs=$(absolute_value $toroidal_mode)
    echo "toroidal_mode_abs="$toroidal_mode_abs


    FLAGS_BASE="-DCONLY -DPFCMOD -DTOKHEAD -DFSTATE -DLEIID=6 -DSTDOUT \
                -DSMALLEQ -DOPENTRACK -DOPENTERM -DPSIT=0.7 \
                -DNOTUNE -DUNBOR="$UNBOR" \
                -DTETALL -DSOLCOL \
                -DRFORCE -DBP -DTIMAX="$timax" -DWREAL -DWLIST"
    if [[ $bplasma == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DB3D -DB3D_EX"
    fi
    if [[ $nohdf5 == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DNOHDF5"
    fi
    if [[ $tokamak == "ITER" ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DTOKAMAK=1"
    elif [[ $tokamak == "STEP" ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DTOKAMAK=10"
    else
        echo "Invalid tokamak."
        exit 1
    fi
    if [[ $bripple == 1 ]]; then
        FLAGS_BASE=$FLAGS_BASE" -DBRIPPLE"
    fi
    echo $FLAGS_BASE

    cp -vf \
    $input_dir"/base.f90" \
    $prec_file 
    cp -f \
    $input_dir"/profile_"$spr_string"_ne.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$spr_string"_ne.dat"
    cp -f \
    $input_dir"/profile_"$spr_string"_Ti.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$spr_string"_Ti.dat"
    cp -f \
    $input_dir"/profile_"$spr_string"_Te.dat" \
    $home_dir"/locust."$tokamak"/InputFiles/profile_"$spr_string"_Te.dat"
    eqdsk_file=$spr_string".eqdsk"
    cp -f \
    $input_dir"/"$eqdsk_file \
    $home_dir"/locust."$tokamak"/InputFiles/."
    ptcle_file=$spr_string"_markers_1000000.dat"
    cp -f \
    $input_dir"/"$ptcle_file \
    $home_dir"/locust."$tokamak"/InputFiles/."

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
    DST="file_profile_ne = 'profile_"$spr_string"_ne.dat'"
    sed -i "s|$SRC|$DST|g" $prec_file
    SRC="file_profile_Te = 'profile_Te_file' ! apkp"
    DST="file_profile_Te = 'profile_"$spr_string"_Te.dat'"
    sed -i "s|$SRC|$DST|g" $prec_file
    SRC="file_profile_Ti = 'profile_Ti_file' ! apkp"
    DST="file_profile_Ti = 'profile_"$spr_string"_Ti.dat'"
    sed -i "s|$SRC|$DST|g" $prec_file

    # BPLASMA code

    if [[ $bplasma == 1 ]]; then

        SRC="nnum   = \[16, -3\] ! apkp - Needs changing"
        DST="nnum   = ["$toroidal_mode"]"
        sed -i "s/$SRC/$DST/g" $prec_file

        SRC="nmde   = 2 ! apkp - Needs changing"
        DST="nmde   = 1"
        sed -i "s/$SRC/$DST/g" $prec_file

        SRC="phase  = \[0.0e0_gpu, 0.0e0_gpu\] ! apkp - Needs changing"
        DST="phase  = [0.0e0_gpu]"
        sed -i "s/$SRC/$DST/g" $prec_file

        # rwm_control and rmp cannot both equal 1
        # Check this:
        if [[ $rwm_control == 1 ]] && [[ $rmp == 1 ]]; then
            echo "rwm_control and rmp cannot both equal 1."
            exit 1
        fi
        if [[ $rwm_control == 1 ]]; then
            BPLASMA_parts=(
                "BPLASMA_cylindrical_tesla_G="$gain_value
                "_bscale="$bscale
                "_n")
        elif [[ $rmp == 1 ]]; then
            BPLASMA_parts=(
                "BPLASMA_"$coil_set
                "_response="$response
                "_current="$(printf "%03d" $current)
                "_100x200"
                "_phase="$(printf "%05.1f" $phase)
                "_n")
        fi
        bplasma_file=$(printf "%s" "${BPLASMA_parts[@]}")
        SRC="'bplasma_file' ! apkp"
        DST="'"$bplasma_file"'"
        sed -i "s|$SRC|$DST|g" $prec_file

        if [[ $rwm_control == 1 ]]; then
            BPLASMA_directory="RWM_control"
        elif [[ $rmp == 1 ]]; then
            BPLASMA_directory="DAVE_v2_100x200/"$coil_set"_n"$toroidal_mode
        fi
        cp -vf \
        $input_dir"/BPLASMA/"$BPLASMA_directory"/"$bplasma_file$toroidal_mode_abs \
        $home_dir"/locust."$tokamak"/InputFiles/."

    fi

    # Read and process pdep_fi.dat file
    pdep_fi_file="$input_dir/ion_info_"$spr_string".dat"
    nion=$(sed -n '1p' "$pdep_fi_file")
    pdep=$(sed -n '2p' "$pdep_fi_file")
    echo $pdep

    # Prepare arrays for fi, Ai, and Zi
    declare -a fi_array=() Ai_array=() Zi_array=()
    while IFS=' ' read -r zi fi ion_name; do
        fi_array+=("$fi")
        Zi_array+=("$zi")
        Ai_array+=("$ion_name") # Directly use the ion name string
    done < <(tail -n +3 "$pdep_fi_file")

    # Convert arrays to comma-separated strings for Fortran
    fi_str=$(IFS=,; echo "${fi_array[*]}")
    Ai_str=$(IFS=,; echo "${Ai_array[*]}")
    Zi_str=$(IFS=,; echo "${Zi_array[*]}")

    # Update prec_file with the new values
    sed -i "s/real( gpu )                                   :: Pdep   = x/real( gpu )                                   :: Pdep   = $pdep/" $prec_file
    sed -i "s/integer,                            parameter :: nion   = x/integer,                            parameter :: nion   = $nion/" $prec_file
    sed -i "s/real( gpu ),  dimension(nion),      parameter :: fi     = \[x\]/real( gpu ),  dimension(nion),      parameter :: fi     = \[$fi_str\]/" $prec_file
    sed -i "s/real( gpu ),  dimension(nion),      parameter :: Ai     = \[x\]/real( gpu ),  dimension(nion),      parameter :: Ai     = \[$Ai_str\]/" $prec_file
    sed -i "s/real( gpu ),  dimension(nion),      parameter :: Zi     = \[x\]/real( gpu ),  dimension(nion),      parameter :: Zi     = \[$Zi_str\]/" $prec_file

    make clean
    make FLAGS="$FLAGS_BASE" -j
    mv -f "locust" "locust_"$run_name"_"$n

done

# Print length of arrays
echo "Length of rcoils: ${#rcoils[@]}"
echo "Length of ncoils: ${#ncoils[@]}"
echo "Length of bripples: ${#bripples[@]}"
echo "Length of toroidal_modes: ${#toroidal_modes[@]}"
echo "Length of coil_sets: ${#coil_sets[@]}"
echo "Length of currents: ${#currents[@]}"
echo "Length of responses: ${#responses[@]}"
echo "Length of phases: ${#phases[@]}"
echo "Length of gain_values: ${#gain_values[@]}"
echo "Length of bscales: ${#bscales[@]}"
echo "Length of bplasmas: ${#bplasmas[@]}"
echo "Length of rwm_controls: ${#rwm_controls[@]}"
echo "Length of rmps: ${#rmps[@]}"
echo "Length of spr_strings: ${#spr_strings[@]}"
