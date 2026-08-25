#!/bin/bash
#
# SMOKE + VERIFICATION - HP, NaCl, paper box 3.0 x 3.0 x 10.0 nm, TIP4P-Ew
#
# First run of the 4-point-water path. 1 replicate, 0.5 ns equil + 1 ns prod.
# The point is the verification block at the end, not the physics.
#
# The failure this catches: if `-w` is left as "TIP3P", build_system_HP takes
#   if water == 'TIP3P': create_interchange(..., charge_from_molecules=[water_TIP3P])
# which assigns 3-site TIP3P charges while tip4p_ew.offxml still adds the M
# virtual site -> water with TIP3P charges AND a charged M-site. Nothing errors.
# The charge check below is the only thing that would notice.
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --partition=blanca-shirts
#SBATCH --qos=blanca-shirts
#SBATCH --account=blanca-shirts
#SBATCH --gres=gpu
#SBATCH --mem=4000m
#SBATCH --job-name=smoke_nacl_paperbox_t4p
#SBATCH --output=slurm_codes/smoke_nacl_paperbox_t4p.%j.log

module purge
ml anaconda
conda activate polymerist-env

export SALT_DATA=structures/salt_data_30x30x100.py
export RESTRAINT=HP
export SALT="Na Cl"
export CENTERIONS="Na Cl"
export WATER="TIP4PEW"                 # must NOT be "TIP3P"
export CONCENTRATIONS="3.5"
export REPNUM=1
export PDB_SUFFIX="p"

export FF_IONS='ionsjc_tip4pew.offxml'
export FF_WATER='tip4p_ew.offxml'
export FF_FILES="$FF_IONS $FF_WATER"

export CASE=smoke_${RESTRAINT}_$(echo ${SALT} | tr -d ' ')_${WATER}_paperbox
export PDB=structures/nacl_35m_${PDB_SUFFIX}r0.pdb

mkdir -p slurm_codes

# ------------------------------------------------------------ pre-flight
for f in "${SALT_DATA}" "${FF_IONS}" "${FF_WATER}" "${PDB}" \
         structures/na.sdf structures/cl.sdf; do
    [ -f "$f" ] || { echo "MISSING: $f"; exit 1; }
done

echo "=== 1. structure on disk (3-site water; M-sites added at runtime) ==="
awk '/^ATOM|^HETATM/ {r=substr($0,18,4); gsub(/ /,"",r); c[r]++; n++}
     END {for (k in c) printf "  %-5s %6d atoms\n", k, c[k]; printf "  TOTAL %6d\n", n}' "${PDB}"
echo "  expect: HOH 9024 (=3008 waters), NA 60, CL 60, TOTAL 9144"
echo "  CRYST1: $(grep -m1 '^CRYST1' ${PDB})"
echo "  expect:  30.000   30.000  100.000"
echo

# ------------------------------------------------------------------ run
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

RC=$?
echo; echo "dispatch exit code: ${RC}"; [ ${RC} -ne 0 ] && exit ${RC}

# ============================================================== VERIFY
echo
echo "################  VERIFICATION  ################"

python - "${CASE}" << 'PYEOF'
import sys, pickle, glob, os
import numpy as np
import openmm

case = sys.argv[1]
fails, warns = [], []

def ok(msg):   print(f"  PASS  {msg}")
def bad(msg):  fails.append(msg); print(f"  FAIL  {msg}")
def warn(msg): warns.append(msg); print(f"  WARN  {msg}")

# ---- 2. built system: particle count, virtual sites, charges --------------
pkl = f"{case}/omm_modified_35.pkl"
if not os.path.exists(pkl):
    bad(f"{pkl} not found - system was never built"); sys.exit(1)

d = pickle.load(open(pkl, "rb"))
system = d['r0']['system']
top    = d['r0']['topology']

print("\n2. built system")
n_part = system.getNumParticles()
print(f"  particles = {n_part}   (expect 12152 = 3008*4 + 120)")
if n_part == 12152: ok("particle count matches a 4-site model")
elif n_part == 9144: bad("9144 particles -> NO virtual sites; tip4p_ew.offxml did not apply")
else: bad(f"unexpected particle count {n_part}")

n_vsite = sum(1 for i in range(n_part) if system.isVirtualSite(i))
print(f"  virtual sites = {n_vsite}   (expect 3008)")
ok("one M-site per water") if n_vsite == 3008 else bad(f"{n_vsite} virtual sites, expected 3008")

nb = [f for f in system.getForces() if isinstance(f, openmm.NonbondedForce)][0]
EC = openmm.unit.elementary_charge
qs = np.array([nb.getParticleParameters(i)[0].value_in_unit(EC) for i in range(n_part)])

# Order-independent: count charge populations rather than assuming O,H,H,M.
# This topology is H,H,O per water with virtual sites appended at the end.
n_H   = int(np.isclose(qs,  0.52422, atol=1e-4).sum())
n_M   = int(np.isclose(qs, -1.04844, atol=1e-4).sum())
n_O   = int(np.isclose(qs,  0.0,     atol=1e-6).sum())
n_cat = int(np.isclose(qs, +1.0,     atol=1e-6).sum())
n_an  = int(np.isclose(qs, -1.0,     atol=1e-6).sum())
n_t3  = int(np.isclose(qs, -0.834,   atol=1e-3).sum())

print(f"  charge populations:")
print(f"    +0.52422 (H)  {n_H:6d}   expect 6016")
print(f"    -1.04844 (M)  {n_M:6d}   expect 3008")
print(f"     0.00000 (O)  {n_O:6d}   expect 3008")
print(f"    +1.0     (Na) {n_cat:6d}   expect 60")
print(f"    -1.0     (Cl) {n_an:6d}   expect 60")
print(f"  first 4 particles = {[round(x,5) for x in qs[:4]]}"
      f"   (H,H,O,H - virtual sites are at the end)")

if n_t3:
    bad(f"{n_t3} particles carry -0.834 -> TIP3P charges; the `-w TIP3P` branch fired")
elif (n_H, n_M, n_O, n_cat, n_an) == (6016, 3008, 3008, 60, 60):
    ok("TIP4P-Ew charges, all populations correct")
else:
    bad(f"charge populations {(n_H, n_M, n_O, n_cat, n_an)} != (6016, 3008, 3008, 60, 60)")

qtot = float(qs.sum())
print(f"  total system charge = {qtot:+.4e}")
ok("electroneutral") if abs(qtot) < 1e-4 else bad(f"net charge {qtot:+.4e}")

# ---- 3. restraint --------------------------------------------------------
print("\n3. restraint")
cef = [f for f in system.getForces() if isinstance(f, openmm.CustomExternalForce)]
if not cef:
    bad("no CustomExternalForce in the system")
else:
    f = cef[0]
    pars = {f.getGlobalParameterName(i): f.getGlobalParameterDefaultValue(i)
            for i in range(f.getNumGlobalParameters())}
    print(f"  restrained particles = {f.getNumParticles()}   (expect 120)")
    print(f"  global parameters    = { {k: round(v,4) for k,v in pars.items()} }")
    ok("120 ions restrained") if f.getNumParticles() == 120 else bad(
        f"{f.getNumParticles()} restrained, expected 120")
    if abs(pars.get('k', 0) - 1.5567) < 1e-3: ok("k = 1.5567 from salt_data")
    else: bad(f"k = {pars.get('k')}, expected 1.5567")
    if abs(pars.get('z0', 0) - 5.0) < 1e-6: ok("z0 = 5.0 nm")
    else: bad(f"z0 = {pars.get('z0')}, expected 5.0")

    # ions must be ions, not M-sites: check masses are Na/Cl-like
    idx = [f.getParticleParameters(i)[0] for i in range(f.getNumParticles())]
    masses = [system.getParticleMass(i).value_in_unit(openmm.unit.dalton) for i in idx]
    print(f"  restrained masses: min {min(masses):.2f}  max {max(masses):.2f}"
          f"   (expect ~22.99 and ~35.45)")
    if min(masses) > 20: ok("restrained particles are ions, not virtual sites")
    else: bad("a restrained particle has ~0 mass - an M-site got restrained")

# ---- 4. trajectory: needs MDAnalysis, which lives in the OTHER env ------
print("\n4. trajectory")
print("  skipped here (MDAnalysis is in osmotic-analysis, not polymerist-env)")
print(f"  run:  conda activate osmotic-analysis && python check_traj.py {case}")

# ---- summary -------------------------------------------------------------
print("\n################  SUMMARY  ################")
if fails:
    print(f"{len(fails)} FAILURE(S) - do not launch production:")
    for m in fails: print(f"  - {m}")
    sys.exit(1)
print("all checks passed" + (f" ({len(warns)} warning(s))" if warns else ""))
for m in warns: print(f"  - {m}")
PYEOF

echo
echo "throughput: 1.55 ns of MD in this job; pull the exact rate with"
echo "  grep -E 'Starting simulation|Completed' slurm_codes/smoke_nacl_paperbox_t4p.*.log"
echo "225 ns/replicate at that rate sets the production wall time."
sacct --format=jobid,jobname,cputime,elapsed
