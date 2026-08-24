#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=70:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=12000m
#SBATCH --job-name=HP_csbr_35m156_tip4p
#SBATCH --output=slurm_codes/HP_csbr_35m156_tip4p.%j.log

module purge
ml anaconda
conda activate polymerist-env

# CASE DETAILS (need to edit each time)
export RESTRAINT=HP
export SALT="Cs Br"
export CENTERIONS="Cs Br"
export WATER="TIP4P_EW"
export CONCENTRATIONS="3.5"
export REPNUM=6

# k (single); if you want multiple ks, replace with --k_values in the call below
export FORCEK=1.56

# New: box geometry
export LZ_NM=10.0     # total Z box length (nm)
export ZCENTER_NM=5.0  # restraint center (nm); default is LZ_NM/2 if omitted

# Force fields
export FF_WATER='tip4p_ew.offxml'
export FF_IONS='ionsjc_tip4pew.offxml'

# shouldn't need to change
export CASE=${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}_156
export FF_FILES="$FF_WATER $FF_IONS"

python osmotic_sim_dispatch_membar_cont.py \
  -n "${CASE}" \
  -s ${SALT} \
  -ci ${CENTERIONS} \
  -w "${WATER}" \
  -rn "${REPNUM}" \
  -m ${CONCENTRATIONS} \
  -ff ${FF_FILES} \
  -r ${RESTRAINT} \
  -k ${FORCEK} \
  -du nanometer \
  --lz_nm "${LZ_NM}" \
  --z_center_nm "${ZCENTER_NM}" \
  --temperature_K 298.15 \
  --equil_nvt_ps 50 \
  --equil_npt_ns 25 \
  --prod_ns 200 \
  --prod_samples 2000\
  --pressure_bar 1.01325 \
  --surface_tension_bar_nm 0.0

echo "CASE=${CASE}  k=${FORCEK}  Lz=${LZ_NM} nm  z_center=${ZCENTER_NM} nm  DATE=$(date -Iseconds)"

sacct --format=jobid,jobname,cputime,elapsed