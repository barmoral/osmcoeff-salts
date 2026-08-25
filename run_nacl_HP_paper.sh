#!/bin/bash
#
# PRODUCTION - Hosseini & Ashbaugh replication
# HP, NaCl, 3.0 x 3.0 x 10.0 nm, TIP4P-Ew + Joung-Cheatham(TIP4P-Ew) ions
#
# Paper procedure                        this run
# -------------------------------------  --------------------------------------
# box 30 x 30 x 100 A                    same (salt_data_paper_30x30x100.py)
# 60 cations + 60 anions, ~3000 waters   60 pairs, 3008 waters
# k tuned so C(0) ~ 3.5 M                k = 1.5567 (--edge-tol 3.9e-4)
# TIP4P/2005 water                       TIP4P-Ew  (JC ions were fit to it)
# Joung-Cheatham ions, LB cross terms    ionsjc_tip4pew.offxml, chi = 0
# 25 C                                   298.15 K
# 2 fs, LINCS water constraints          2 fs, constraints from tip4p_ew.offxml
# LJ 9 A + long-range correction, PME    cutoff 9 A / PME (set in the offxml)
# 25 ns equilibration                    50 ps NVT + 25 ns
# 200 ns production                      200 ns, 2000 frames (100 ps)
# NVT, water count hand-tuned to         MonteCarloMembraneBarostat, XY iso,
#   rho_reservoir = 0.997                  Z fixed, 1 atm, zero surface tension
#
# Deliberate deviations:
#   * TIP4P-Ew instead of TIP4P/2005. JC ions were parameterised against
#     TIP4P-Ew, so this pairing is internally consistent; the paper transferred
#     them into TIP4P/2005. Expect small differences in absolute Pi.
#   * Membrane barostat instead of hand-tuning the water count.
#   * 6 independent replicates instead of one long trajectory.
#
# THIS SCRIPT IS DESIGNED TO BE SUBMITTED 6 TIMES AS A CHAIN - see the bottom.
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=4000m
#SBATCH --job-name=HP_nacl_paperbox_t4p
#SBATCH --output=slurm_codes/HP_nacl_paperbox_t4p.%j.log

module purge
ml anaconda
conda activate polymerist-env

export SALT_DATA=structures/salt_data_30x30x100.py
export RESTRAINT=HP
export SALT="Na Cl"
export CENTERIONS="Na Cl"
export WATER="TIP4PEW"                # must NOT be "TIP3P"
export CONCENTRATIONS="3.5"
export REPNUM=6                       # ALL replicates; each job does what fits
export PDB_SUFFIX="p"
export PROD_FRAMES=2000               # must equal --prod_samples below

export FF_IONS='ionsjc_tip4pew.offxml'
export FF_WATER='tip4p_ew.offxml'
export FF_FILES="$FF_IONS $FF_WATER"

export CASE=${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}_paperbox

mkdir -p slurm_codes

# --------------------------------------------------------- pre-flight
for f in "${SALT_DATA}" "${FF_IONS}" "${FF_WATER}" \
         structures/na.sdf structures/cl.sdf; do
    [ -f "$f" ] || { echo "MISSING: $f"; exit 1; }
done
NPDB=$(ls structures/nacl_35m_${PDB_SUFFIX}r?.pdb 2>/dev/null | wc -l)
if [ "${NPDB}" -lt 6 ]; then
    echo "Only ${NPDB} of 6 replicate pdbs found. Build them all first:"
    echo "  python build_salt_data.py packmol --salt-data ${SALT_DATA} \\"
    echo "      --salt NaCl --molality 3.5 --replicates 6 --suffix ${PDB_SUFFIX}"
    exit 1
fi

# ---------------------------------------------------- resume cleanup
# A replicate killed by the wall clock leaves all three of its folders behind,
# and simulation_exists() would then treat it as finished. Delete any replicate
# whose production trajectory is short, so the next job in the chain redoes it.
echo "=== resume check ==="
python resume_cleanup.py "${CASE}" "${PROD_FRAMES}"
echo

echo "=== inputs ==="
for f in structures/nacl_35m_${PDB_SUFFIX}r?.pdb; do
    printf "  %-30s %s\n" "$(basename $f)" \
        "$(awk '/^ATOM|^HETATM/ {r=substr($0,18,4); gsub(/ /,"",r); c[r]++}
                END {printf "HOH %d  NA %d  CL %d", c["HOH"], c["NA"], c["CL"]}' $f)"
done
echo "  (expect HOH 9024, NA 60, CL 60 in every file)"
echo

# ---------------------------------------------------------------- run
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
  --time_step_fs 2.0 \
  --equil_nvt_ps 50 \
  --equil_npt_ns 25 \
  --prod_ns 200 \
  --prod_samples ${PROD_FRAMES} \
  --pressure_bar 1.01325 \
  --surface_tension_bar_nm 0.0

echo
echo "=== state at end of this job ==="
python resume_cleanup.py "${CASE}" "${PROD_FRAMES}" --report-only
echo "CASE=${CASE}  DATE=$(date -Iseconds)"
sacct -j ${SLURM_JOB_ID} --format=jobid,jobname,cputime,elapsed,maxrss

# =====================================================================
# HOW TO CHAIN THIS
# =====================================================================
#
# Submit six copies of this same file. SLURM holds five of them and releases
# them one at a time, each starting only after the previous finishes.
#
#   jid=$(sbatch --parsable run_nacl_HP_paperbox_tip4pew_PROD.sh)
#   echo "job 1: $jid"
#   for i in 2 3 4 5 6; do
#       jid=$(sbatch --parsable --dependency=afterany:$jid \
#             run_nacl_HP_paperbox_tip4pew_PROD.sh)
#       echo "job $i: $jid"
#   done
#
# --parsable          sbatch prints just the number, so $jid captures it
# --dependency=afterany:$jid   wait for that job to end, whatever its exit code
#                              (afterok would stop the chain on any crash)
#
# Nothing is edited between submissions. REPNUM=6 means each job walks r0..r5,
# skips whatever is already finished, works on the next one, and gets killed by
# the wall clock partway through a replicate. The resume-cleanup block deletes
# that partial replicate so the following job restarts it cleanly.
#
# MONITOR
#   squeue -u $USER                  # PD + reason (Dependency) = queued behind
#   python resume_cleanup.py ${CASE} 2000 --report-only
#   grep -c "Completed" slurm_codes/HP_nacl_paperbox_t4p.*.log
#
# CANCEL THE WHOLE CHAIN
#   scancel -u $USER --name=HP_nacl_paperbox_t4p
#
# DONE when the report shows 6 replicates at 2000 frames. Then:
#   conda activate osmotic-analysis
#   python check_traj.py ${CASE}     # expect meas/ideal ~0.85, not 0.71
#   python run_HP_analysis.py