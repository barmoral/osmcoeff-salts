#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=12000m
#SBATCH --job-name=HP_nacl_db35m050smoke
#SBATCH --output=slurm_codes/HP_nacl_db35m050smoke.%j.log

module purge
module avail
ml anaconda
conda activate polymerist-env

# CASE DETAILS (need to edit each time)
export RESTRAINT=HP
export SALT="Na Cl"
export CENTERIONS="Na Cl"
export WATER="TIP3P"
export CONCENTRATIONS="3.5"
export REPNUM=1
export PDB_SUFFIX="d"

# k (single); if you want multiple ks, replace with --k_values in the call below
export FORCEK=0.50

# New: box geometry
export LZ_NM=28.8     # total Z box length (nm)
export ZCENTER_NM=14.4  # restraint center (nm); default is LZ_NM/2 if omitted

# Force fields
export FF_WATER='tip3p.offxml'
export FF_IONS='openff-2.3.0.offxml'

# shouldn't need to change
export CASE=${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}_050_${PDB_SUFFIX}box_smoke
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
  --pdb_suffix "${PDB_SUFFIX}" \
  --temperature_K 298.15 \
  --equil_npt_ns 0.5 \
  --prod_ns 1 \
  --prod_samples 20 \
  --pressure_bar 1.01325 \
  --surface_tension_bar_nm 0.0

echo "CASE=${CASE}  k=${FORCEK}  Lz=${LZ_NM} nm  z_center=${ZCENTER_NM} nm  DATE=$(date -Iseconds)"
sacct --format=jobid,jobname,cputime,elapsed