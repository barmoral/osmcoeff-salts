#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=8000m
#SBATCH --job-name=HP_nacl_35m_smoke
#SBATCH --output=slurm_codes/HP_nacl_35m_smoke.%j.log

module purge
ml anaconda
conda activate polymerist-env

# CASE DETAILS (need to edit each time)
export RESTRAINT=HP
export SALT="Na Cl"
export CENTERIONS="Na Cl"
export WATER="TIP3P"
export CONCENTRATIONS="3.5"
export REPNUM=1
export SALT_DATA=structures/salt_data_48x48x144.py
export PDB_SUFFIX=""

# New: box geometry
export LZ_NM=14.4     # total Z box length (nm)
export ZCENTER_NM=7.2  # restraint center (nm); default is LZ_NM/2 if omitted

# Force fields
export FF_WATER='tip3p.offxml'
export FF_IONS='openff-2.3.0.offxml'
export FF_FILES="$FF_WATER $FF_IONS"

# shouldn't need to change
export CASE=${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}_smoke

mkdir -p slurm_codes
 
# ------------------------------------------------------- pre-flight checks
# Fail fast rather than 20 minutes in.
for f in "${SALT_DATA}" osmotic_sim_dispatch_membar_cont.py \
         structures/nacl_35m_r0.pdb structures/na.sdf structures/cl.sdf \
         ${FF_WATER} ${FF_IONS}; do
    if [ ! -f "$f" ]; then
        echo "MISSING: $f"
        echo
        echo "If it is the PDB, generate it first:"
        echo "  python build_salt_data.py packmol --salt-data ${SALT_DATA} \\"
        echo "      --salt NaCl --molality 3.5 --replicates 1 --outdir packmol_inputs"
        echo "  cd packmol_inputs && packmol < build_nacl_35m_r0.inp && bash add_cryst1.sh"
        echo "  mv nacl_35m_r0.pdb ../structures/"
        exit 1
    fi
done
 
echo "=== structure check ==="
awk '/^ATOM|^HETATM/ {r=substr($0,18,4); gsub(/ /,"",r); c[r]++}
     END {for (k in c) printf "  %-5s %6d atoms\n", k, c[k]}' \
    structures/nacl_35m_r0.pdb
echo "  (expect HOH 33267 atoms = 11089 waters, NA 236, CL 236)"
echo "  CRYST1: $(grep -m1 '^CRYST1' structures/nacl_35m_r0.pdb)"
echo
 
# ------------------------------------------------------------------- run
python -u osmotic_sim_dispatch_membar_cont.py \
  -n "${CASE}" \
  -s ${SALT} \
  -ci ${CENTERIONS} \
  -w "${WATER}" \
  -rn "${REPNUM}" \
  -m ${CONCENTRATIONS} \
  -ff ${FF_FILES} \
  -r ${RESTRAINT} \
  --salt_data "${SALT_DATA}" \
  --pdb_suffix "${PDB_SUFFIX}" \
  -du nanometer \
  --equil_nvt_ps 50 \
  --equil_npt_ns 0.5 \
  --prod_ns 1 \
  --prod_samples 20 \
  --pressure_bar 1.01325 \
  --surface_tension_bar_nm 0.0
 
echo
echo "CASE=${CASE} salt_data=${SALT_DATA}  DATE=$(date -Iseconds)"
sacct --format=jobid,jobname,cputime,elapsed