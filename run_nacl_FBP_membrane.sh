#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=160:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=10000m
#SBATCH --job-name=FBP_nacl_4184
#SBATCH --output=slurm_codes/FBP_nacl_4184.%j.log

module purge
module avail
ml anaconda
conda activate polymerist-env

# CASE DETAILS (need to edit each time)
export RESTRAINT=FBP
export SALT="Na Cl"
export CENTERIONS="Na Cl"
export WATER="TIP3P"
export CONCENTRATIONS="0.5 1.0 1.4 2.0 2.5 3.0 3.5 4.0"
export REPNUM=20

# k (single); if you want multiple ks, replace with --k_values in the call below
export FORCEK=4184

# New: box geometry
export LZ_NM=14.4     # total Z box length (nm)
export ZCENTER_NM=7.2  # restraint center (nm); default is LZ_NM/2 if omitted

# Force fields
export FF_WATER='tip3p.offxml'
export FF_IONS='openff-2.3.0.offxml'

# shouldn't need to change
export CASE=${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}_4184
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
  --delta_z 2.4

echo "CASE=${CASE}  k=${FORCEK}  Lz=${LZ_NM} nm  z_center=${ZCENTER_NM} nm  DATE=$(date -Iseconds)"

sacct --format=jobid,jobname,cputime,elapsed