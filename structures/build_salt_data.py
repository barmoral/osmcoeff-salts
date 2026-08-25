
#!/usr/bin/env python3
"""
build_salt_data.py
==================
Turns salt_reference.py (experimental data, box-independent) into one
salt_data_<tag>.py per simulation box, holding everything the HP and FBP
pipelines need. Can also emit the matching packmol inputs.
 
    salt_reference.py  ->  build  ->  salt_data_<tag>.py
                                          |
                                       packmol  ->  *.inp  ->  structures/*.pdb
 
Usage
-----
Build a box file:
 
    python build_salt_data.py build --box 4.8 4.8 14.4 --tag 48x48x144
 
Build + packmol inputs for one system (6 replicates):
 
    python build_salt_data.py build   --box 4.8 4.8 14.4 --tag 48x48x144
    python build_salt_data.py packmol --salt-data salt_data_48x48x144.py \
        --salt NaCl --molality 3.5 --replicates 6
    cd packmol_inputs && for f in build_*.inp; do packmol < $f; done
    bash add_cryst1.sh && mv nacl_35m_r*.pdb ../structures/
 
Longer box, its own restraint (normal build):
 
    python build_salt_data.py build --box 4.8 4.8 28.8 --tag 48x48x288
 
Longer box, IDENTICAL restraint and N (controlled box-length experiment only):
 
    python build_salt_data.py build --box 4.8 4.8 28.8 \
        --clone-from salt_data_48x48x144.py --tag 48x48x288_clone
 
Thicker reservoirs / fewer ions (tighter edge tolerance):
 
    python build_salt_data.py build --box 4.8 4.8 28.8 --edge-tol 1e-6 --tag dbox_tight
 
Second replicate set with different seeds, own pdb_suffix:
 
    python build_salt_data.py packmol --salt-data salt_data_48x48x288.py \
        --salt NaCl --molality 3.5 --replicates 6 --suffix d --seed0 5000
 
Design
------
Targets are MOLALITY (mol/kg water), whole-box basis - the same convention as
convert_profile_to_molal() in the analysis, so design and results are directly
comparable. Molality also makes the design immune to the barostat squeezing XY:
m(0) = C(0)/f, and both scale as 1/A, so the area cancels.
 
Per box, ONE k_HP and ONE delta_z_FBP serve every concentration; only the
particle count changes:
 
  k_HP     from the box alone. Require C(edge)/C(0) <= edge_tol:
               K = -ln(edge_tol) / (Lz/2)^2 ,  k = 2RT*K
           Default 1e-3 reproduces values already in use (Lz=14.4 -> 0.6606,
           close to the historical 0.68; Lz=28.8 -> 0.165).
 
  N(m)     linear in molality, from the Gaussian normalisation (erf included):
               N = m * 2w * MOLAR_TO_NM3 * mass_water / (Lz * 1e-24) / erf(...)
           with w = int_0^inf exp(-K z^2) dz = 0.5*sqrt(pi/K).
 
  delta_z  equals w exactly, for the same N and target. HP spreads N ions over
           half-width w as a Gaussian; FBP spreads the same N over half-width w
           as a top hat. That is why ONE .pdb serves both methods - same box,
           same N, only the restraint differs.
 
Notes
-----
* With k_FBP_wall = 4184 the decay length is 0.024 nm, <1% of ions outside.
* HP uses the ideal-solution Gaussian. Measured peaks run ~10-15% below target
  (CsBr: measured/ideal = 0.86). Target m/0.86, or calibrate from a short run.
* T = 298.15 K only. salt_reference.py declares REFERENCE_TEMPERATURE and the
  build refuses a mismatch: phi and density are 25 C measurements.
"""

import argparse
import re
from pathlib import Path
from math import erf, exp, log, pi, sqrt

NA = 6.02214076e23
R_KJ = 8.31446261815324e-3      # kJ/(mol K)
M_WATER = 18.01528              # g/mol
MOLAR_TO_NM3 = NA * 1e-24       # mol/L -> nm^-3

MOLAR_MASS = {
    "CsBr": 212.809, "CsCl": 168.358, "CsI": 259.810, "CsNO3": 194.912,
    "KBr": 119.002, "KCl": 74.551, "KI": 166.003, "KNO3": 101.103,
    "LiBr": 86.845, "LiCl": 42.394, "LiClO4": 106.392, "LiI": 133.845,
    "LiNO3": 68.946, "NaBr": 102.894, "NaCl": 58.443, "NaI": 149.894,
    "NaNO3": 84.995, "RbBr": 165.372, "RbCl": 120.921, "RbNO3": 147.473,
    "NH4Cl": 53.491, "NH4NO3": 80.043,
    "Na2SO4": 142.042, "NH42SO4": 132.140, "MgCl": 95.211,
    "Na2HP": 141.958, "NH42HP": 132.056,
}
VANT_HOFF = {s: 2 for s in MOLAR_MASS}
for _s in ("Na2SO4", "NH42SO4", "MgCl", "Na2HP", "NH42HP"):
    VANT_HOFF[_s] = 3

K_FBP_WALL = 4184.0             # kJ/mol/nm^2, flat-bottom wall stiffness


# ---------------------------------------------------------------------------
def box_properties(lx, ly, lz, temperature, edge_tol, water_molarity):
    """Everything that depends on the box but not on the salt."""
    volume_nm3 = lx * ly * lz
    n_water = round(volume_nm3 * water_molarity * MOLAR_TO_NM3)
    mass_water_kg = n_water * M_WATER / NA * 1e-3
    volume_L = volume_nm3 * 1e-24

    K = -log(edge_tol) / (lz / 2) ** 2               # nm^-2
    k_hp = 2 * R_KJ * temperature * K                # kJ/mol/nm^2
    sigma = sqrt(1 / (2 * K))
    w = 0.5 * sqrt(pi / K)                           # int_0^inf exp(-K z^2) dz
    erf_corr = erf(sqrt(K) * lz / 2)

    return dict(
        lx=lx, ly=ly, lz=lz, area=lx * ly, volume_nm3=volume_nm3,
        volume_L=volume_L, n_water=n_water, mass_water_kg=mass_water_kg,
        temperature=temperature, edge_tol=edge_tol,
        reference_temperature=298.15,
        water_molarity=water_molarity,
        k_HP=k_hp, delta_z_FBP=w, k_FBP_wall=K_FBP_WALL,
        z_center=lz / 2, sigma=sigma, four_sigma=4 * sigma,
        erf_correction=erf_corr,
        molal_factor=mass_water_kg / volume_L,       # kg water / L solution
    )


def n_pairs(molality, box):
    """Ion pairs giving the target peak (HP) / plateau (FBP) molality.

    Independent of cross-sectional area - see module docstring.
    """
    w = box["delta_z_FBP"]
    n = (molality * 2 * w * MOLAR_TO_NM3 * box["mass_water_kg"]
         / (box["lz"] * 1e-24) / box["erf_correction"])
    return round(n)


def realised_molality(n, box):
    """Inverse of n_pairs - what a given N actually delivers."""
    w = box["delta_z_FBP"]
    return (n * box["erf_correction"] * box["lz"] * 1e-24
            / (2 * w * MOLAR_TO_NM3 * box["mass_water_kg"]))


# ---------------------------------------------------------------------------
def load_clone_source(path):
    """Import a previously generated salt_data_*.py and return (BOX, rows)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_clone_src", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.BOX, {(d["salt"], d["molality"]): d for d in mod.salt_infos}


def load_reference(path):
    """Import salt_reference.py (or a legacy salt_data.py) and return its rows.

    A proper import for salt_reference.py; falls back to regex-scraping for the
    old hand-written format so the original file can still be re-read once.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("_ref", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if hasattr(mod, "SALT_REFERENCE"):
        return mod.SALT_REFERENCE, getattr(mod, "MOLAR_MASS", MOLAR_MASS), \
               getattr(mod, "VANT_HOFF", VANT_HOFF)
    rows = [d for d in getattr(mod, "salt_infos", []) if "salt" in d]
    return rows, MOLAR_MASS, VANT_HOFF


def parse_legacy(path):
    rows, _, _ = load_reference(path)
    return [{k: str(v) for k, v in r.items()} for r in rows]


HEADER = '''"""
{fname}  (generated by build_salt_data.py - do not edit by hand)
{under}

Box:            {lx} x {ly} x {lz} nm    A = {area:.3f} nm^2
Waters to pack: {n_water}   (at {water_molarity} M)
Temperature:    {temperature} K
Edge tolerance: C(edge)/C(0) <= {edge_tol}

Restraint parameters - fixed for this box, valid at EVERY target molality:

    HP   k_HP        = {k_HP:.4f} kJ/mol/nm^2     z0 = {z_center} nm
         sigma       = {sigma:.3f} nm   (4 sigma = {four_sigma:.3f} nm)
    FBP  delta_z     = {delta_z_FBP:.4f} nm            z0 = {z_center} nm
         k_FBP_wall  = {k_FBP_wall} kJ/mol/nm^2

Only `num_particles` changes with the target concentration, so ONE packmol
build per (salt, molality) serves BOTH methods: same box, same ion count, the
restraint is applied afterwards by the dispatch script.

`molality` is the TARGET peak concentration (HP: at z0) or plateau
concentration (FBP: between the walls), NOT the nominal box composition.
Whole-box molal basis, matching convert_profile_to_molal() in the analysis.

{clone_note}`osmotic_coefficient` and `density` are experimental values for that molality
(Hamer & Wu / Robinson & Stokes; densities IAPWS-consistent, 25 C).
`molarity` is derived, not stored: C = 1000*m*rho/(1000 + m*MM).
"""

from dataclasses import dataclass

NA = 6.02214076e23
M_WATER = 18.01528
MOLAR_TO_NM3 = NA * 1e-24

BOX = {box_repr}


@dataclass(frozen=True)
class SaltData:
    salt: str
    molality: float               # target peak / plateau, mol/kg
    osmotic_coefficient: float    # experimental at this molality
    density: float                # g/cm3, experimental solution density
    num_particles: int            # ion pairs for THIS box
    k_HP: float                   # kJ/mol/nm^2, harmonic restraint
    delta_z_FBP: float            # nm, flat-bottom half width
    k_FBP_wall: float             # kJ/mol/nm^2, flat-bottom wall stiffness

    @property
    def molarity(self) -> float:
        mm = MOLAR_MASS[self.salt]
        return 1000.0 * self.molality * self.density / (1000.0 + self.molality * mm)


def lookup(salt, molality, tol=1e-6):
    """Return the SaltData row for a salt at a target molality."""
    for d in salt_infos:
        if d["salt"] == salt and abs(d["molality"] - molality) < tol:
            return SaltData(**d)
    raise KeyError(f"{{salt}} at {{molality}} mol/kg not in {fname}")


'''


def _verify_output(fname, box):
    """Re-import the file we just wrote and check it round-trips.

    Catches the whole class of "serialiser mangled a number" bugs: rounding a
    1e-22 quantity to 0.0, truncating precision, a mis-typed key, a row whose
    num_particles does not match what n_pairs() would compute now.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location("_verify_probe", fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    problems = []
    for key, val in box.items():
        if not isinstance(val, float):
            continue
        got = mod.BOX.get(key)
        if got is None:
            problems.append(f"BOX['{key}'] missing from file")
        elif val == 0.0:
            if got != 0.0:
                problems.append(f"BOX['{key}']: file {got!r} != computed 0.0")
        elif abs(got - val) / abs(val) > 1e-12:
            extra = "  <-- collapsed to zero" if got == 0.0 else ""
            problems.append(f"BOX['{key}']: file {got!r} != computed {val!r}{extra}")

    if mod.salt_infos:
        row = mod.SaltData(**mod.salt_infos[0])
        want = n_pairs(row.molality, box)
        if row.num_particles != want:
            problems.append(
                f"num_particles {row.salt}@{row.molality}: file {row.num_particles} != {want}")
        if abs(row.k_HP - box["k_HP"]) > 1e-12:
            problems.append(f"row k_HP {row.k_HP} != BOX k_HP {box['k_HP']}")
        if abs(row.delta_z_FBP - box["delta_z_FBP"]) > 1e-12:
            problems.append("row delta_z_FBP != BOX delta_z_FBP")
        _ = row.molarity          # exercises MOLAR_MASS lookup + the property

    if problems:
        raise SystemExit("\n*** {} FAILED ROUND-TRIP CHECK ***\n  {}\n"
                         .format(fname, "\n  ".join(problems)))
    print(f"  round-trip  OK ({len(mod.salt_infos)} rows re-imported, BOX values exact)")


def build(args):
    global _ref_mod
    import importlib.util as _ilu
    _spec = _ilu.spec_from_file_location("_ref_probe", args.legacy)
    _ref_mod = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_ref_mod)

    box = box_properties(*args.box, args.temperature, args.edge_tol,
                         args.water_molarity)

    T_ref = getattr(_ref_mod, "REFERENCE_TEMPERATURE", None)
    if T_ref is not None and abs(args.temperature - T_ref) > 0.5:
        msg = (f"\n*** TEMPERATURE MISMATCH ***\n"
               f"    simulation T = {args.temperature} K  (sets k_HP)\n"
               f"    reference  T = {T_ref} K  (osmotic coefficients and densities)\n"
               f"    k_HP will be correct, but every experimental comparison in the\n"
               f"    analysis will be against {T_ref} K data. To run at another\n"
               f"    temperature you need a new salt_reference.py: densities are\n"
               f"    tabulated at other T, but phi(T) requires a Pitzer model -\n"
               f"    Hamer & Wu is 25 C only.\n")
        if not args.allow_temperature_mismatch:
            raise SystemExit(msg + "    Re-run with --allow-temperature-mismatch to override.\n")
        print(msg)

    clone_rows = None
    if args.clone_from:
        src_box, clone_rows = load_clone_source(args.clone_from)
        # Keep the restraint EXACTLY as in the source box; only geometry changes.
        K_src = src_box["k_HP"] / (2 * R_KJ * box["temperature"])
        box["k_HP"] = src_box["k_HP"]
        box["delta_z_FBP"] = src_box["delta_z_FBP"]
        box["sigma"] = sqrt(1 / (2 * K_src))
        box["four_sigma"] = 4 * box["sigma"]
        box["erf_correction"] = erf(sqrt(K_src) * box["lz"] / 2)
        box["edge_tol"] = exp(-K_src * (box["lz"] / 2) ** 2)
        box["cloned_from"] = args.clone_from
        box["cloned_from_lz"] = src_box["lz"]
    tag = args.tag or f"{args.box[0]}x{args.box[1]}x{args.box[2]}".replace(".", "")
    fname = args.out or f"salt_data_{tag}.py"

    # NOTE: do NOT round() here - mass_water_kg ~ 1e-22 and volume_L ~ 1e-22
    # would both collapse to 0.0. repr() keeps full float precision.
    box_repr = "{\n" + "".join(
        f"    {k!r}: {v!r},\n" for k, v in box.items()) + "}"

    clone_note = ""
    if clone_rows is not None:
        clone_note = (
            "!! CONTROLLED-EXPERIMENT FILE - NOT A NORMAL BUILD !!\n"
            f"k_HP, delta_z_FBP and num_particles were copied verbatim from\n"
            f"{args.clone_from} (Lz = {box['cloned_from_lz']} nm) so that ONLY the\n"
            "box length differs. Use it to test whether a longer box changes the HP\n"
            "result or the FBP interface artefacts. Do not use it as a template for\n"
            "new systems - regenerate without --clone-from for those.\n\n")
    head = HEADER.format(fname=fname, under="=" * (len(fname) + 44),
                         box_repr=box_repr, clone_note=clone_note, **box)

    lines = [head, "MOLAR_MASS = {"]
    lines += [f'    "{s}": {mm},' for s, mm in sorted(MOLAR_MASS.items())]
    lines += ["}\n", "VANT_HOFF = {"]
    lines += [f'    "{s}": {VANT_HOFF[s]},' for s in sorted(MOLAR_MASS)]
    lines += ["}\n", "salt_infos = ["]

    recs = parse_legacy(args.legacy)
    skipped = 0
    for r in recs:
        salt, m = r["salt"], float(r["molality"])
        if clone_rows is not None:
            src = clone_rows.get((salt, m))
            if src is None:
                skipped += 1
                continue
            n = src["num_particles"]
        else:
            n = n_pairs(m, box)
            if n < args.min_pairs:
                skipped += 1
                continue
        lines += [
            "    {",
            f'        "salt": "{salt}",',
            f'        "molality": {m},',
            f'        "osmotic_coefficient": {float(r["osmotic_coefficient"])},',
            f'        "density": {float(r["density"])},',
            f'        "num_particles": {n},',
            f'        "k_HP": {box["k_HP"]!r},',
            f'        "delta_z_FBP": {box["delta_z_FBP"]!r},',
            f'        "k_FBP_wall": {box["k_FBP_wall"]},',
            "    },",
        ]
    lines.append("]")

    with open(fname, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    _verify_output(fname, box)
    print(f"wrote {fname}")
    print(f"  box        {box['lx']} x {box['ly']} x {box['lz']} nm, "
          f"A = {box['area']:.3f} nm^2, {box['n_water']} waters")
    print(f"  HP         k = {box['k_HP']:.4f} kJ/mol/nm^2, sigma = {box['sigma']:.3f} nm, "
          f"4sigma = {box['four_sigma']:.3f} nm (half-box {box['lz']/2:.2f})")
    print(f"  FBP        delta_z = {box['delta_z_FBP']:.4f} nm, wall k = {box['k_FBP_wall']}")
    print(f"  rows       {len(recs) - skipped} kept, {skipped} dropped (< {args.min_pairs} pairs)")
    print(f"  z0         {box['z_center']} nm")
    if clone_rows is not None:
        print(f"  CLONED     restraint + N copied verbatim from {args.clone_from}")
        print(f"             (Lz {box['cloned_from_lz']} -> {box['lz']} nm; "
              f"effective edge tol now {box['edge_tol']:.2e})")
        print(f"             HP reservoir  = {box['lz']/2 - box['four_sigma']:.2f} nm/side")
        print(f"             FBP reservoir = {box['lz']/2 - box['delta_z_FBP']:.2f} nm/side")
    print()
    print(f"  {'target m':>9} {'N pairs':>8}   (same N for HP and FBP)")
    for m in (0.5, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0):
        n = n_pairs(m, box)
        print(f"  {m:>9} {n:>8}   -> realised {realised_molality(n, box):.3f} mol/kg")


PACKMOL_TEMPLATE = """#
# {salt} at {molality} mol/kg target  ({method_note})
# box {lx} x {ly} x {lz} nm | {n_water} waters | {n_pairs} ion pairs
# generated by build_salt_data.py from {sdfile} -- do not hand-edit
#
tolerance {tolerance}
filetype pdb
seed {seed}
output {outfile}

{water_blocks}
structure {cat_pdb}
  number {n_pairs}
  inside box {m} {m} {ion_lo:.1f} {xmax:.1f} {ymax:.1f} {ion_hi:.1f}
end structure

structure {an_pdb}
  number {n_pairs}
  inside box {m} {m} {ion_lo:.1f} {xmax:.1f} {ymax:.1f} {ion_hi:.1f}
end structure
"""

WATER_BLOCK = """structure {water_pdb}
  number {n}
  inside box {m} {m} {zlo:.1f} {xmax:.1f} {ymax:.1f} {zhi:.1f}
end structure
"""

CRYST1 = "CRYST1{a:>9.3f}{b:>9.3f}{c:>9.3f}{al:>7.2f}{be:>7.2f}{ga:>7.2f} P 1           1"


def packmol(args):
    """Emit one packmol .inp per replicate, straight from a salt_data file."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("_sd", args.salt_data)
    sd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sd)
    box = sd.BOX

    row = sd.lookup(args.salt, args.molality)
    n = row.num_particles
    lx_a, ly_a, lz_a = box["lx"] * 10, box["ly"] * 10, box["lz"] * 10
    margin = args.margin

    # Ions start uniformly over +/- 1.5 * delta_z: wider than the FBP slab and
    # than the HP core, so the initial density is below BOTH targets and packmol
    # has room. Both restraints pull the cloud into shape within a few hundred ps.
    half = min(1.5 * row.delta_z_FBP * 10, lz_a / 2 - margin)
    ion_lo, ion_hi = lz_a / 2 - half, lz_a / 2 + half

    # water in ~50 A slabs so packmol converges quickly
    n_slabs = max(1, round(lz_a / 50.0))
    per = [box["n_water"] // n_slabs] * n_slabs
    per[-1] += box["n_water"] - sum(per)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mi = _mi(args.molality)
    written = []
    for r in range(args.replicates):
        blocks = ""
        for i, npw in enumerate(per):
            zlo = margin if i == 0 else i * lz_a / n_slabs
            zhi = lz_a - margin if i == n_slabs - 1 else (i + 1) * lz_a / n_slabs
            blocks += WATER_BLOCK.format(water_pdb=args.water_pdb, n=npw, m=margin,
                                         xmax=lx_a - margin, ymax=ly_a - margin,
                                         zlo=zlo, zhi=zhi)
        name = f"{args.salt.lower()}_{mi}m_{args.suffix}r{r}.pdb"
        inp = PACKMOL_TEMPLATE.format(
            salt=args.salt, molality=args.molality,
            method_note="same file serves HP and FBP",
            lx=box["lx"], ly=box["ly"], lz=box["lz"],
            n_water=box["n_water"], n_pairs=n, sdfile=Path(args.salt_data).name,
            tolerance=args.tolerance, seed=args.seed0 + r, outfile=name,
            water_blocks=blocks, cat_pdb=args.cation_pdb, an_pdb=args.anion_pdb,
            m=margin, xmax=lx_a - margin, ymax=ly_a - margin,
            ion_lo=ion_lo, ion_hi=ion_hi)
        f = outdir / f"build_{args.salt.lower()}_{mi}m_{args.suffix}r{r}.inp"
        f.write_text(inp)
        written.append((f, name))

    cryst = CRYST1.format(a=lx_a, b=ly_a, c=lz_a, al=90, be=90, ga=90)
    (outdir / "add_cryst1.sh").write_text(
        "#!/bin/bash\n# insert the correct CRYST1 record into every packmol output\n"
        f'CRYST="{cryst}"\n'
        f'for f in {args.salt.lower()}_{mi}m_{args.suffix}r*.pdb; do\n'
        '  grep -q "^CRYST1" "$f" && sed -i "/^CRYST1/d" "$f"\n'
        '  sed -i "1i $CRYST" "$f"\n'
        '  echo "$f: $(grep -c HOH "$f") HOH  $(head -1 "$f")"\n'
        "done\n")

    print(f"{len(written)} packmol inputs -> {outdir}/")
    print(f"  {args.salt} @ {args.molality} mol/kg : {n} pairs + {box['n_water']} waters")
    print(f"  ions packed in z = {ion_lo:.1f} .. {ion_hi:.1f} A  (z0 = {lz_a/2:.1f})")
    print(f"  k_HP = {row.k_HP}   delta_z_FBP = {row.delta_z_FBP}   wall k = {row.k_FBP_wall}")
    print(f"  outputs: {written[0][1]} .. {written[-1][1]}")
    print(f"\n  cd {outdir} && for f in build_*.inp; do packmol < $f; done && bash add_cryst1.sh")
    print(f"  then move the .pdb files into  structures/")
    print(f"\n  expected per file: {{'HOH': {box['n_water']}, "
          f"'{args.salt[:2].upper()}': {n}, ...}}")


def _mi(mol):
    mi1 = f"{mol:.1f}"
    return str(int(mol)) if mol % 1 == 0 else mi1.replace(".", "")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build", help="generate salt_data_<tag>.py for a box")
    b.add_argument("--reference", dest="legacy", default="salt_reference.py",
                   help="experimental data file (default: salt_reference.py)")
    b.add_argument("--box", nargs=3, type=float, required=True,
                   metavar=("LX", "LY", "LZ"))
    b.add_argument("--tag", default=None)
    b.add_argument("--out", default=None)
    b.add_argument("--temperature", type=float, default=298.15)
    b.add_argument("--edge-tol", dest="edge_tol", type=float, default=1e-3)
    b.add_argument("--water-molarity", dest="water_molarity", type=float, default=55.5)
    b.add_argument("--min-pairs", dest="min_pairs", type=int, default=10)
    b.add_argument("--allow-temperature-mismatch", dest="allow_temperature_mismatch",
                   action="store_true",
                   help="build anyway when simulation T differs from the reference T")
    b.add_argument("--clone-from", dest="clone_from", default=None,
                   help="ONE-OFF controlled experiment: copy k_HP, delta_z_FBP and "
                        "num_particles verbatim from another salt_data_*.py so that "
                        "ONLY the box length changes. Do not use for normal builds.")
    b.set_defaults(func=build)

    k = sub.add_parser("packmol", help="emit packmol inputs from a salt_data file")
    k.add_argument("--salt-data", dest="salt_data", required=True)
    k.add_argument("--salt", required=True, help="e.g. CsBr")
    k.add_argument("--molality", type=float, required=True)
    k.add_argument("--replicates", type=int, default=6)
    k.add_argument("--suffix", default="", help="pdb_suffix, e.g. 'd' for the double box")
    k.add_argument("--outdir", default="packmol_inputs")
    k.add_argument("--water-pdb", dest="water_pdb", default="../water.pdb")
    k.add_argument("--cation-pdb", dest="cation_pdb", default=None)
    k.add_argument("--anion-pdb", dest="anion_pdb", default=None)
    k.add_argument("--tolerance", type=float, default=2.0)
    k.add_argument("--margin", type=float, default=0.5, help="wall margin in A")
    k.add_argument("--seed0", type=int, default=1000)
    k.set_defaults(func=packmol)

    args = p.parse_args()
    if getattr(args, "cmd", None) == "packmol":
        import re as _re
        parts = _re.findall(r"[A-Z][a-z0-9]*", args.salt)
        if args.cation_pdb is None:
            args.cation_pdb = f"../{parts[0].lower()}.pdb"
        if args.anion_pdb is None:
            args.anion_pdb = f"../{''.join(parts[1:]).lower()}.pdb"
    args.func(args)


if __name__ == "__main__":
    main()