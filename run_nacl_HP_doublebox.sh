#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=70:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=6000m
#SBATCH --job-name=HP_nacl_db35m
#SBATCH --output=slurm_codes/HP_nacl_db35m.%j.log

module purge
ml anaconda
conda activate polymerist-env

# CASE DETAILS (need to edit each time)
export RESTRAINT=HP
export SALT="Na Cl"
export CENTERIONS="Na Cl"
export WATER="TIP3P"
export CONCENTRATIONS="3.5"
export REPNUM=6
export SALT_DATA=structures/salt_data_48x48x288_clone.py
export PDB_SUFFIX="d"


# Force fields
export FF_WATER='tip3p.offxml'
export FF_IONS='openff-2.3.0.offxml'

# shouldn't need to change
export CASE=${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}_${PDB_SUFFIX}box
export FF_FILES="$FF_WATER $FF_IONS"

python -u osmotic_sim_dispatch_membar_cont.py \
  -n "${CASE}" \
  -s ${SALT} \
  -ci ${CENTERIONS} \
  -w "${WATER}" \
  -rn "${REPNUM}" \
  -m ${CONCENTRATIONS} \
  -ff ${FF_FILES} \
  -r ${RESTRAINT} \
  -du nanometer \
  --salt_data "${SALT_DATA}" \
  --pdb_suffix "${PDB_SUFFIX}" \
  --temperature_K 298.15 \
  --equil_nvt_ps 50 \
  --equil_npt_ns 5 \
  --prod_ns 20 \
  --prod_samples 1000 \
  --pressure_bar 1.01325 \
  --surface_tension_bar_nm 0.0

echo "CASE=${CASE} Lz=${LZ_NM} nm  z_center=${ZCENTER_NM} nm  DATE=$(date -Iseconds)"
sacct --format=jobid,jobname,cputime,elapsed