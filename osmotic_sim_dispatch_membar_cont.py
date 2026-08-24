"""
openMM code to calculate osmotic pressure and osmotic coefficients from harmonic or flat-bottom
potentials, using a membrane barostat with XY isotropic scaling, fixed Z, and zero surface tension.

Important:
- Lz (total Z box length) is specified by --lz_nm (default 28.8).
- z_center (restraint center) is specified by --z_center_nm (default Lz/2).
  For backward compatibility, --z_center is an alias of --z_center_nm.

"""

# =========================
# Standard imports & config
# =========================
import logging
import warnings
warnings.filterwarnings('ignore')

from typing import Optional, List
import argparse
from argparse import Namespace

import os
import os.path
import pickle
from pathlib import Path
import numpy as np
import math
from tqdm import tqdm
import subprocess
from builtins import sum
from tqdm.notebook import tqdm as tqdm_nb  # keep if you sometimes use notebooks

# =========================
# OpenMM / OpenFF imports
# =========================
import openmm
from openmm import CustomExternalForce
from openmm.app import Topology as OMMTopology
from openmm import MonteCarloMembraneBarostat, Vec3  # Membrane barostat (XY iso, Z fixed)

from openff.toolkit import Molecule, Topology
from openff.toolkit import ForceField
from openff.interchange import Interchange

# Your helper packages
from polymerist.genutils.fileutils.pathutils import assemble_path
from polymerist.genutils.decorators.functional import allow_string_paths

from polymerist.mdtools.openfftools import topology
from polymerist.mdtools.openfftools import boxvectors
from polymerist.mdtools.openfftools import TKWRAPPERS, GTR
from polymerist.mdtools.openfftools.partialcharge.molchargers import NAGLCharger

from polymerist.mdtools.openfftools.unitsys import openff_to_openmm
from polymerist.mdtools.openfftools.solvation.solvents import water_TIP3P

from polymerist.mdtools.openmmtools.execution import run_simulation_schedule
from polymerist.mdtools.openmmtools.parameters import (
    SimulationParameters, ThermoParameters, IntegratorParameters, ReporterParameters, ThermostatParameters
)
from polymerist.mdtools.openmmtools.evaluation import get_context_positions

from openff.units import UnitRegistry
ureg = UnitRegistry()

# Sci/plot (kept as in your original imports)
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from scipy.optimize import least_squares, minimize
from scipy.integrate import simpson, quad, trapezoid
from timeit import default_timer as timer

# Units
from openmm.unit import (
    bar, mole, litre, kelvin, kilojoule_per_mole, nanometer, angstrom,
    kilocalorie_per_mole, kilogram, molar, atmosphere, nanosecond, picosecond,
    femtoseconds, Quantity, Unit, AVOGADRO_CONSTANT_NA, BOLTZMANN_CONSTANT_kB
)

# Additional imports already used above
import importlib.util

def load_salt_data(path):
    spec = importlib.util.spec_from_file_location("salt_data_mod", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

# =================================
# Logging (console + file)
# =================================
log_path = Path("failed_replicates.log")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(log_path, mode="a"),
    ],
)

LOGGER = logging.getLogger(__name__)

# ================================================
# Utilities for Z harmonization & membrane barostat
# ================================================
def _vec_length(v: Vec3) -> float:
    return float((v.x*v.x + v.y*v.y + v.z*v.z) ** 0.5)

def _scale_vec_to_length(v: Vec3, new_len: float) -> Vec3:
    curr = _vec_length(v)
    if curr == 0:
        raise ValueError("Periodic box vector has zero length; cannot scale.")
    s = new_len / curr
    return Vec3(v.x * s, v.y * s, v.z * s)

class ReplicateZRegistry:
    """
    Target Z length per concentration key; enforce on Systems and Simulations so
    all replicates share the same Lz.
    """
    def __init__(self):
        self._target_by_key: dict[str, float] = {}

    def set_target_for_key(self, key: str, z_length_nm: float):
        """Explicitly set desired Lz (nm) for a concentration key."""
        self._target_by_key[key] = float(z_length_nm)

    def _ensure_target_from_system(self, system: openmm.System, key: str):
        if key not in self._target_by_key:
            a, b, c = system.getDefaultPeriodicBoxVectors()
            self._target_by_key[key] = _vec_length(c)

    def enforce_on_system(self, system: openmm.System, key: str):
        """
        Scale System's default c-vector length to the registered target.
        If no target is registered yet, infer from this system.
        """
        self._ensure_target_from_system(system, key)
        a, b, c = system.getDefaultPeriodicBoxVectors()
        c_new = _scale_vec_to_length(c, self._target_by_key[key])
        system.setDefaultPeriodicBoxVectors(a, b, c_new)

    def enforce_on_simulation(self, simulation, key: str):
        """
        Scale the Context c-vector length to the registered target.
        If no target is registered yet, infer from this Simulation.
        """
        if key not in self._target_by_key:
            a, b, c = simulation.context.getState(getPositions=False).getPeriodicBoxVectors()
            self._target_by_key[key] = _vec_length(c)
        a, b, c = simulation.context.getState(getPositions=False).getPeriodicBoxVectors()
        c_new = _scale_vec_to_length(c, self._target_by_key[key])
        simulation.context.setPeriodicBoxVectors(a, b, c_new)

class MembraneBarostatManager:
    """
    Attach MonteCarloMembraneBarostat with XY isotropic, Z fixed, configurable pressure/surface tension,
    and attempt frequency. Reinitialize the Context so the Simulation sees the new Force.
    """
    def __init__(self,
                 pressure_bar: float = 1.01325,
                 temperature_K: float = 300.0,
                 surface_tension_bar_nm: float = 0.0,
                 frequency: int = 100):
        self.pressure_bar = pressure_bar
        self.temperature_K = temperature_K
        self.surface_tension_bar_nm = surface_tension_bar_nm
        self.frequency = frequency

    def add_to_simulation(self, simulation):
        barostat = MonteCarloMembraneBarostat(
            self.pressure_bar * bar,
            self.surface_tension_bar_nm * bar * nanometer,
            self.temperature_K * kelvin,
            MonteCarloMembraneBarostat.XYIsotropic,
            MonteCarloMembraneBarostat.ZFixed,
            self.frequency
        )
        simulation.system.addForce(barostat)
        # Required after modifying the System so the Simulation sees the change.
        simulation.context.reinitialize(preserveState=True)
        return barostat


# =========================
# Salt info loader (as-is)
# =========================

def load_salt_info(SD, ion1, ion2):
    salt = ion1 + ion2
    print(f"Salt to be analyzed: {salt}")
    entries = [SD.SaltData(**e) if isinstance(e, dict) else e for e in SD.salt_infos]
    filtered = {
        f"Molality {e.molality} mol/kg": {
            "Molality": e.molality,
            "Molarity": round(e.molarity, 4),        # derived property now
            "Number of Particles": e.num_particles,  # no math.ceil - already int
            "Osmotic Coefficient": e.osmotic_coefficient,
            "Density": e.density,
            "k_HP": e.k_HP,
            "delta_z_FBP": e.delta_z_FBP,
            "k_FBP_wall": e.k_FBP_wall,
        }
        for e in entries if e.salt == salt
    }
    return {salt: filtered} if filtered else {"Error": f"No data found for {salt}"}


# ======================================
# System builders + restraint application
# ======================================
def modify_omm_flat_bottom(concentration, repnum, input_picklefile, ion1, ion2, wdir,
                           center_atom1, center_atom2, k, delta_z, z_center, expected_pairs=None):
    with open(input_picklefile, "rb") as f:
        omm_build = pickle.load(f)

    omm_objects_mod={}
    
    for r in tqdm(range(repnum), desc=f"Modifying {concentration}m systems", leave=False):
        rep_key = f'r{r}'
        print(rep_key)

        omm_top = omm_build[rep_key]["topology"]
        omm_pos = omm_build[rep_key]["positions"]
        omm_sys = omm_build[rep_key]["system"]

        # identify central atoms
        poly_atoms = []
        for residue in omm_top.residues():
            if residue.name == ion1.upper():
                for atom in residue.atoms():
                    if atom.element.symbol == f"{center_atom1}":
                        poly_atoms.append(atom.index)
            if residue.name == ion2.upper():
                for atom in residue.atoms():
                    if atom.element.symbol == f"{center_atom2}":
                        poly_atoms.append(atom.index)

        print("CENTER ATOMS", len(poly_atoms), poly_atoms)

        if expected_pairs is not None and len(poly_atoms) != 2 * expected_pairs:
            raise RuntimeError(
                f"{rep_key}: PDB has {len(poly_atoms)} restrained ions, "
                f"salt_data expects {2 * expected_pairs} "
                f"({expected_pairs} pairs). Structure file does not match "
                f"the selected salt_data.")

        # flat-bottom potential along z
        fb_force = CustomExternalForce('0.5*k*(max(0, abs(z-z0)-rbf)^2)')
        fb_force.addGlobalParameter('k', k)
        fb_force.addGlobalParameter('rbf', delta_z)
        fb_force.addGlobalParameter('z0', z_center)   # <-- z_center is the center point (e.g., 14.4 nm)

        for atom_index in poly_atoms:
            fb_force.addParticle(atom_index, [])

        omm_sys.addForce(fb_force)

        omm_objects_mod[rep_key] = {
            'topology': omm_top,
            'positions': omm_pos,
            'system': omm_sys
        }
        
    output_picklefile = f"{wdir}/omm_modified_{concentration}.pkl"
    with open(output_picklefile, "wb") as f:
        pickle.dump(omm_objects_mod, f)
    print(f"✅ System dictionary for {concentration}m saved as pickle file")

    return omm_objects_mod


def build_system_FBP(concentration, repnum, wdir, ion1, ion2, ff, water, salt_dict,
                     center_atom1, center_atom2, k, delta_z, z_center, lz_total, expected_pairs=None, pdb_suffix=''):
    """
    Construct Interchange with a tetragonal box whose Z length is 'lz_total' (e.g., 28.8 nm).
    The restraint center 'z_center' (e.g., 14.4 nm) is NOT used here; it's applied by modify_*.
    """
    sdf_path1=f'structures/{ion1.lower()}.sdf'
    sdf_path2=f'structures/{ion2.lower()}.sdf'

    if concentration % 1 == 0:
        mi1 = f"{concentration:.1f}"
        mi = str(int(concentration))
    else:
        mi = f"{concentration:.1f}".replace('.', '')
        mi1 = f"{concentration:.1f}"

    print(f'Molality of maximum concentration = {mi} mol/kg',
            salt_dict[f'{ion1}{ion2}'][f'Molality {mi1} mol/kg'])
    
    input_picklefile = f"{wdir}/omm_build_{mi}{pdb_suffix}.pkl"

    if os.path.exists(input_picklefile):
        print("File already exists. Skipping code.")
        with open(input_picklefile, "rb") as f:
            omm_build = pickle.load(f)
        return modify_omm_flat_bottom(mi, repnum, input_picklefile, ion1, ion2, wdir,
                                      center_atom1, center_atom2, k, delta_z, z_center, expected_pairs=expected_pairs)
   
    else:
        print("File not found. Running code...")

        omm_builds = {}

        for r in tqdm(range(repnum), desc=f"Building {mi}m systems", leave=False):
            pdb_path = f'structures/{ion1.lower()}{ion2.lower()}_{mi}m_{pdb_suffix}r{r}.pdb'
            POL1 = Molecule.from_file(sdf_path1)
            POL2 = Molecule.from_file(sdf_path2)
            
            off_top = Topology.from_pdb(pdb_path, unique_molecules=[POL1, POL2])

            if water == 'TIP3P':
                inc = ff.create_interchange(
                    topology=off_top,
                    toolkit_registry=GTR,
                    charge_from_molecules=[water_TIP3P]
                )
            else:
                inc = ff.create_interchange(
                    topology=off_top,
                    toolkit_registry=GTR
                )

            def make_bbox_tetragonal_with_Lz(bbox, Lz_total_openmm):
                """
                bbox: BoxVectorsQuantity (PintQuantity 3x3) from get_topology_bbox(off_top)
                Returns: BoxVectorsQuantity tetragonal reduced-form with Z length = Lz_total_openmm.
                """
                # Use your helper to normalize/handle flexible boxes (Pint)
                bbox = boxvectors.box_vectors_flexible(bbox)

                u = bbox.units
                B = bbox.m_as(u)  # numeric 3x3 in unit u

                Lx, Ly, _ = np.linalg.norm(B, axis=1) * u
                Lxy = max(Lx, Ly)

                # Convert target Lz (OpenMM Quantity) into bbox's Pint units
                # one_u_openmm = openff_to_openmm(1.0 * u)
                q = openff_to_openmm(1.0 * u)
                one_u_openmm = q.unit
                Lz_in_u = Lz_total_openmm.value_in_unit(one_u_openmm)
                Lz_target = Lz_in_u * u

                xyz = np.array([Lxy.m_as(u), Lxy.m_as(u), Lz_target.m_as(u)]) * u
                return boxvectors.xyz_to_box_vectors(xyz)

            # starting box from topology and set Lz to lz_total
            box = boxvectors.get_topology_bbox(off_top)
            inc.box = make_bbox_tetragonal_with_Lz(box, lz_total)

            # Verify in Interchange (Pint)
            u = inc.box.units
            B = inc.box.m_as(u)
            lens = np.linalg.norm(B, axis=1) * u
            lens_nm = [openff_to_openmm(L).value_in_unit(nanometer) for L in lens]

            z_target_nm = lz_total.value_in_unit(nanometer)
            print("inc.box lengths (nm):", lens_nm)
            assert abs(lens_nm[2] - z_target_nm) < 1e-6, f"Lz is {lens_nm[2]} nm, expected {z_target_nm} nm"
            assert abs(lens_nm[0] - lens_nm[1]) < 1e-6, f"Lx != Ly ({lens_nm[0]} vs {lens_nm[1]})"

            # Export to OpenMM
            omm_top = inc.to_openmm_topology(collate=True)
            omm_pos = openff_to_openmm(inc.get_positions(include_virtual_sites=True))
            omm_sys = inc.to_openmm_system(combine_nonbonded_forces=False, add_constrained_forces=True)

            omm_builds[f"r{r}"] = {
                'topology': omm_top,
                'positions': omm_pos,
                'system': omm_sys
            }

        with open(input_picklefile, "wb") as f:
            pickle.dump(omm_builds, f)
        print(f"✅ System dictionary for {mi1}m saved as pickle file")
        
        return modify_omm_flat_bottom(mi, repnum, input_picklefile, ion1, ion2, wdir,
                                      center_atom1, center_atom2, k, delta_z, z_center, expected_pairs=expected_pairs)


def modify_omm_harmonic(concentration, repnum, input_picklefile, ion1, ion2, wdir,
                        center_atom1, center_atom2, k, delta_z, z_center, expected_pairs=None):
    with open(input_picklefile, "rb") as f:
        omm_build = pickle.load(f)

    omm_objects_mod={}
    
    for r in tqdm(range(repnum), desc=f"Modifying {concentration}m systems", leave=False):
        rep_key = f'r{r}'
        print(rep_key)

        omm_top = omm_build[rep_key]["topology"]
        omm_pos = omm_build[rep_key]["positions"]
        omm_sys = omm_build[rep_key]["system"]

        poly_atoms = []
        for residue in omm_top.residues():
            if residue.name == ion1.upper():
                for atom in residue.atoms():
                    if atom.element.symbol == f"{center_atom1}":
                        poly_atoms.append(atom.index)
            if residue.name == ion2.upper():
                for atom in residue.atoms():
                    if atom.element.symbol == f"{center_atom2}":
                        poly_atoms.append(atom.index)

        print("CENTER ATOMS", len(poly_atoms), poly_atoms)

        if expected_pairs is not None and len(poly_atoms) != 2 * expected_pairs:
            raise RuntimeError(
                f"{rep_key}: PDB has {len(poly_atoms)} restrained ions, "
                f"salt_data expects {2 * expected_pairs} "
                f"({expected_pairs} pairs). Structure file does not match "
                f"the selected salt_data.")

        # Harmonic restraint centered at z0
        fb_force = CustomExternalForce('0.5*k*((z-z0)^2)')
        fb_force.addGlobalParameter('k', k)
        fb_force.addGlobalParameter('z0', z_center)  # <-- center point (e.g., 14.4 nm)

        for atom_index in poly_atoms:
            fb_force.addParticle(atom_index, [])

        omm_sys.addForce(fb_force)

        omm_objects_mod[rep_key] = {
            'topology': omm_top,
            'positions': omm_pos,
            'system': omm_sys
        }
        
    output_picklefile = f"{wdir}/omm_modified_{concentration}.pkl"
    with open(output_picklefile, "wb") as f:
        pickle.dump(omm_objects_mod, f)
    print(f"✅ System dictionary for {concentration}m saved as pickle file")

    return omm_objects_mod


def build_system_HP(concentration, repnum, ion1, ion2, wdir, ff, water, salt_dict,
                    center_atom1, center_atom2, k, delta_z, z_center, lz_total,expected_pairs=None, pdb_suffix=''):
    """
    HP builder mirrors the FBP builder but does not change topology box shape beyond using get_topology_bbox.
    The harmonic restraint itself is applied in modify_omm_harmonic() with center z0 = z_center.
    """
    sdf_path1=f'structures/{ion1.lower()}.sdf'
    sdf_path2=f'structures/{ion2.lower()}.sdf'

    if concentration % 1 == 0:
        mi1 = f"{concentration:.1f}"
        mi = str(int(concentration))
    else:
        mi = f"{concentration:.1f}".replace('.', '')
        mi1 = f"{concentration:.1f}"

    print(f'Molality of maximum concentration = {mi} mol/kg',
            salt_dict[f'{ion1}{ion2}'][f'Molality {mi1} mol/kg'])
    
    input_picklefile = f"{wdir}/omm_build_{mi}{pdb_suffix}.pkl"

    if os.path.exists(input_picklefile):
        print("File already exists. Skipping code.")
        with open(input_picklefile, "rb") as f:
            omm_build = pickle.load(f)
        return modify_omm_harmonic(mi, repnum, input_picklefile, ion1, ion2, wdir,
                                   center_atom1, center_atom2, k, delta_z, z_center, expected_pairs=expected_pairs)
   
    else:
        print("File not found. Running code...")

        omm_builds = {}

        for r in tqdm(range(repnum), desc=f"Building {mi}m systems", leave=False):
            pdb_path = f'structures/{ion1.lower()}{ion2.lower()}_{mi}m_{pdb_suffix}r{r}.pdb'
            POL1 = Molecule.from_file(sdf_path1)
            POL2 = Molecule.from_file(sdf_path2)
            
            off_top = Topology.from_pdb(pdb_path, unique_molecules=[POL1, POL2])

            if water == 'TIP3P':
                inc = ff.create_interchange(
                    topology=off_top,
                    toolkit_registry=GTR,
                    charge_from_molecules=[water_TIP3P]
                )
            else:
                inc = ff.create_interchange(
                    topology=off_top,
                    toolkit_registry=GTR
                )

            # Leave box from topology (HP path); we still enforce final Lz via registry before Simulations
            inc.box = boxvectors.get_topology_bbox(off_top)

            omm_top = inc.to_openmm_topology()
            omm_pos = openff_to_openmm(inc.get_positions(include_virtual_sites=True))
            omm_sys = inc.to_openmm_system(combine_nonbonded_forces=False, add_constrained_forces=True)

            omm_builds[f"r{r}"] = {
                'topology': omm_top,
                'positions': omm_pos,
                'system': omm_sys
            }

        with open(input_picklefile, "wb") as f:
            pickle.dump(omm_builds, f)
        print(f"✅ System dictionary for {mi1}m saved as pickle file")
        
        return modify_omm_harmonic(mi, repnum, input_picklefile, ion1, ion2, wdir,
                                   center_atom1, center_atom2, k, delta_z, z_center, expected_pairs=expected_pairs)


# ===================================
# Filesystem helpers
# ===================================
def simulation_exists(base_dir, postfix):
    """Check if both equil and prod simulation folders already exist."""
    base_dir = str(base_dir)
    equil_NVT_dir = os.path.join(base_dir, f"equil_sim_NVT_{postfix}")
    equil_NPT_dir = os.path.join(base_dir, f"equil_sim_NPT_{postfix}")
    prod_dir = os.path.join(base_dir, f"prod_sim_{postfix}")
    return os.path.exists(equil_NVT_dir) and os.path.exists(equil_NPT_dir) and os.path.exists(prod_dir)


# ---- Copy/convert helpers (simplified & robust) ----
import shutil

def copy_file(source_path, destination_path):
    """
    Copy a file preserving metadata (portable), creating parent dirs as needed.
    """
    src = Path(source_path)
    dst = Path(destination_path)
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        shutil.copy2(src, dst)
        print(f"✅ File copied from {src} to {dst}")
    except Exception as e:
        print(f"❌ Error copying file: {e}")


def convert_dcd_to_xtc(input_file, output_file, topology=None, stride=None, force=True):
    """
    Converts a DCD trajectory file to XTC format using mdconvert (MDTraj CLI).
    """
    input_file = str(input_file)
    output_file = str(output_file)
    cmd = ['mdconvert', '-o', output_file, input_file]
    if force:
        cmd.insert(1, '-f')  # add -f before -o for mdconvert
    if stride is not None:
        cmd.extend(['-s', str(stride)])
    if topology is not None:
        cmd.extend(['-t', str(topology)])

    if not force and os.path.exists(output_file):
        os.remove(output_file)

    try:
        res = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"🔁 Successfully converted {input_file} to {output_file}")
        if res.stderr:
            print(res.stderr)
    except subprocess.CalledProcessError as e:
        print(f"❌ Error during conversion: {e}")
        if e.stderr:
            print(e.stderr)


def _format_molality(mol: float):
    """
    Returns (mi, mi1):
      - mi : '2' for 2.0, '25' for 2.5, '64' for 6.4
      - mi1: '2.0', '2.5', '6.4' (for logging)
    """
    if mol % 1 == 0:
        mi1 = f"{mol:.1f}"
        mi = str(int(mol))
    else:
        mi1 = f"{mol:.1f}"
        mi = mi1.replace('.', '')
    return mi, mi1


def copy_simulation_outputs(
    molalities,
    N_replicates,
    wdir,
    rdir,
    *,
    ktags,                   # <- explicit list of ktags (e.g., ['k0p64','k0p74'])
    topology_from='pdb',
    stride=None,
    force_convert=True,
):
    """
    Collect PDBs and convert DCD→XTC for each (molality, ktag, replicate).

    Source (must exist):
      wdir/
        prod_sim_{mi}m_{ktag}_r{i}/
          prod_sim_{mi}m_{ktag}_r{i}_topology.pdb
          prod_sim_{mi}m_{ktag}_r{i}_trajectory.dcd

    Destination:
      rdir/{ktag}_{mi}m/
        md{mi}m_r{i}.pdb
        md{mi}m_r{i}.xtc
    """
    wdir = Path(wdir).resolve()
    rdir = Path(rdir)
    rdir.mkdir(exist_ok=True, parents=True)

    for mol in molalities:
        mi, mi1 = _format_molality(mol)

        for ktag in ktags:
            dest_dir = rdir / f'{ktag}_{mi}m'
            dest_dir.mkdir(parents=True, exist_ok=True)

            for i in range(N_replicates):
                spost = f'{mi}m_{ktag}_r{i}'
                source_dir = wdir / f'prod_sim_{spost}'
                source_pdb = source_dir / f'prod_sim_{spost}_topology.pdb'
                input_dcd  = source_dir / f'prod_sim_{spost}_trajectory.dcd'

                pdb_name = f'md{mi}m_r{i}.pdb'
                xtc_name = f'md{mi}m_r{i}.xtc'
                dest_pdb = dest_dir / pdb_name
                dest_xtc = dest_dir / xtc_name

                print(f"\n📂 Checking: {source_pdb}")
                if source_pdb.exists():
                    copy_file(source_pdb, dest_pdb)
                else:
                    print(f"❌ File '{source_pdb}' does not exist.")

                print(f"🔍 Looking for DCD file: {input_dcd}")
                if input_dcd.exists():
                    topo_for_mdconvert = str(dest_pdb) if (topology_from == 'pdb' and dest_pdb.exists()) else None
                    convert_dcd_to_xtc(
                        str(input_dcd),
                        str(dest_xtc),
                        topology=topo_for_mdconvert,
                        stride=stride,
                        force=force_convert
                    )
                else:
                    print(f"❌ DCD file '{input_dcd}' not found.")



# ===================================
# Restraint registry
# ===================================
RESTRAINT_TYPES = {
    'FBP' : build_system_FBP,
    'HP' : build_system_HP,
}
RESTRAINT_TYPE_ALIASES = {
    'FBP' : 'flat-bottomed potential',
    'HP' : 'harmonic potential',
}


# =========================
# CLI argument parsing
# =========================
def parse_args() -> Namespace:
    """Read user inputs and preprocess quantities."""
    parser = argparse.ArgumentParser(
        description='Run OpenMM simulation schedules according to presets of simulation parameters'
    )

    parser.add_argument('-n', '--name', required=True,
                        help='Tag to use to refer to working directory and internal files')
    parser.add_argument('-s', '--salt', nargs='+', required=True,
                        help='Two ion names, e.g., "NA CL"')
    parser.add_argument('-ci', '--centerion', nargs='+', required=True,
                        help='Two center atom symbols for the polyatomic ions, e.g., "P B"')
    parser.add_argument('-w', '--water', required=True,
                        help='Water model to be used, e.g., TIP3P')
    parser.add_argument('-rn', '--repnumber', required=True,
                        help='Number of replicates to run')
    parser.add_argument('-m', '--concentrations', nargs='+', required=True,
                        help='List of concentrations to be analyzed, e.g., 2.0 3.0')

    parser.add_argument('-c', '--cwd', type=Path, default=Path.cwd(),
                        help='Working directory base path')
    parser.add_argument('-spp', '--sim_param_paths', type=Path, nargs='+', required=False,
                        help='Optional serialized simulation presets (unused, preserved)')
    parser.add_argument('-ff', '--ff_files', nargs='+', required=True,
                        help='Force field OFFXML files to be used (paths or registry-relative)')

    parser.add_argument('-r', '--restraint_type', choices=['HP', 'FBP'], required=True,
                        help='Potential restraint type to apply to ions')

    parser.add_argument('--pdb_suffix', type=str, default='',
                    help='Suffix inserted before replicate number in PDB filename. '
                         'E.g., "d" gives nacl_2m_dr0.pdb (double-box). Default: empty.')
            
    parser.add_argument('--salt_data', required=True,
                    help='per-box salt_data file, e.g. salt_data_48x48x144.py')

    parser.add_argument('-k', help='Single spring constant for the restraint (kJ/mol/nm^2) — optional')

    # ---- Z geometry controls ----
    parser.add_argument('--lz_nm', type=float, default=None,
                        help='Total Z box length Lz in nm from --salt_data BOX.')
    parser.add_argument('--z_center_nm', type=float, default=None,
                        help='Center point z0 in nm for restraints (default Lz/2).')
    # Backward-compatible alias:
    parser.add_argument('--z_center', dest='z_center_nm', type=float, default=None,
                        help='(Alias) Center point z0 in nm for restraints (default Lz/2).')

    parser.add_argument('--delta_z', type=float, default=None,
                        help='Half-width for flat-bottom potential (nm) from --salt_data row')

    parser.add_argument('-du', '--distance_unit', choices=['angstrom', 'angstroms', 'nanometer', 'nanometers'],
                        default='nanometer', help='Unit convention for distance inputs')

    # Barostat and thermal controls
    parser.add_argument('--pressure_bar', type=float, default=1.01325,
                        help='Target lateral pressure (bar) for the membrane barostat (XY).')
    parser.add_argument('--surface_tension_bar_nm', type=float, default=0.0,
                        help='Surface tension gamma (bar*nm). Use 0 for zero-tension runs.')
    parser.add_argument('--barostat_freq', type=int, default=100,
                        help='Attempt frequency (timesteps) for Monte Carlo volume changes.')
    parser.add_argument('--temperature_K', type=float, default=None,
                        help='Simulation temperature in Kelvin; used for integrators and barostat. from --salt_data BOX')

    # Schedules and sampling
    parser.add_argument('--time_step_fs', type=float, default=2.0, help='MD timestep (fs).')
    parser.add_argument('--equil_nvt_ps', type=float, default=50.0, help='NVT equilibration length (ps).')
    parser.add_argument('--equil_npt_ns', type=float, default=3.0, help='Second-phase equilibration length (ns).')
    parser.add_argument('--prod_ns', type=float, default=20.0, help='Production length (ns).')
    parser.add_argument('--nvt_samples', type=int, default=100, help='# samples in NVT equil.')
    parser.add_argument('--prod_samples', type=int, default=1000, help='# samples in production.')

    # Optional extra FF
    parser.add_argument('--extra_ff', type=Path, default=None,
                        help='Optional extra OFFXML to append; if omitted, none is added.')

    args = parser.parse_args()

    # Environment fallbacks (kept from your original)
    salt_env = os.getenv('SALT')
    if not args.salt or args.salt == [""]:
        if salt_env:
            args.salt = salt_env.split()
        else:
            raise ValueError("Salt must be provided via '--salt' or SALT env var.")

    ci_env = os.getenv('CENTERIONS')
    if not args.centerion or args.centerion == [""]:
        if ci_env:
            args.centerion = ci_env.split()
        else:
            raise ValueError("Center ions must be provided via '--centerion' or CENTERIONS env var.")
        
    concs_env = os.getenv('CONCENTRATIONS')
    if not args.concentrations or args.concentrations == [""]:
        if concs_env:
            args.concentrations = concs_env.split()
        else:
            raise ValueError("List of concentrations must be provided via '--concentrations' or CONCENTRATIONS env var.")

    ff_env = os.getenv('FF_FILES')    
    if not args.ff_files or args.ff_files == [""]:
        if ff_env:
            args.ff_files = ff_env.split()
        else:
            raise ValueError("Force field files must be provided via '--ff_files' or FF_FILES env var.")

    # ---- resolve the per-box data file -------------------------------
    SD = load_salt_data(args.salt_data)
    args.SD = SD
    args.box = SD.BOX
    ion1, ion2 = args.salt[0], args.salt[1]
    args.row = SD.lookup(f"{ion1}{ion2}", float(args.concentrations[0]))

    if args.lz_nm is None:
        args.lz_nm = args.box["lz"]
    if args.z_center_nm is None:
        args.z_center_nm = args.box["z_center"]
    if args.temperature_K is None:
        args.temperature_K = args.box["temperature"]
    if args.delta_z is None:
        args.delta_z = args.row.delta_z_FBP

    # derived paths
    args.working_dir = args.cwd / args.name

    # restraint selection
    args.restraint_fn = RESTRAINT_TYPES[args.restraint_type]
    logging.info(f'Using restraint type {args.restraint_type} ({RESTRAINT_TYPE_ALIASES[args.restraint_type]})')

    # units
    args.distance_unit = getattr(openmm.unit, args.distance_unit)
    args.delta_z = args.delta_z * args.distance_unit

    # Convert schedule knobs to quantities
    args.temperature = args.temperature_K * kelvin
    args.time_step = args.time_step_fs * femtoseconds
    args.equil_nvt_time = args.equil_nvt_ps * picosecond
    args.equil_npt_time = args.equil_npt_ns * nanosecond
    args.prod_time = args.prod_ns * nanosecond

    # Z geometry: Lz total and z_center for restraints
    args.lz = args.lz_nm * args.distance_unit
    if args.z_center_nm is None:
        args.z_center = 0.5 * args.lz  # default: center of the box
    else:
        args.z_center = args.z_center_nm * args.distance_unit

    # Force constant handling:
    FORCE_CONST_UNIT = (kilojoule_per_mole / nanometer ** 2)
    if args.k is not None:
        args.k = float(args.k)                      # explicit override
    else:
        args.k = (args.row.k_HP if args.restraint_type == 'HP'
                  else args.row.k_FBP_wall)         # from salt_data
    args.k_quantity = args.k * FORCE_CONST_UNIT
    args.ktag = _k_to_tag(args.k)                   # e.g. 0.6606 -> 'k0p661'

    logging.info(
        f"[salt_data] {args.salt_data}: N={args.row.num_particles} pairs, "
        f"k={args.k} ({args.ktag}, {args.restraint_type}), "
        f"delta_z={args.delta_z}, Lz={args.lz_nm}, z0={args.z_center_nm}, "
        f"T={args.temperature_K}")

    return args


# =========================
# Main
# =========================

def _k_to_tag(k) -> str:
    """Turn a numeric k into a short tag.
    
    Integers -> 'k1'
    Floats   -> 'k0p64'
    """
    # Treat integer-valued floats as integers
    if isinstance(k, int) or (isinstance(k, float) and k.is_integer()):
        return f"k{int(k)}"

    # Non-integer float formatting
    s = f"{k:.3g}".replace('.', 'p').replace('-', 'm')
    return f"k{s}"

def main():
    args = parse_args()
    ion1 = args.salt[0]
    ion2 = args.salt[1]
    water = args.water

    SD = load_salt_data(args.salt_data)
    BOX = SD.BOX
    pdb_suffix=args.pdb_suffix

    concentration_list = [float(x) for x in args.concentrations]
    print("Concentrations:", concentration_list)
    print(f"k (kJ/mol/nm^2): {args.k}   tag: {args.ktag}")
    print(f"Lz (nm): {args.lz.value_in_unit(nanometer):.3f}, z_center (nm): {args.z_center.value_in_unit(nanometer):.3f}")

    repnumber = int(args.repnumber)
    salt_dict = load_salt_info(args.SD, ion1, ion2)

    wdir = args.working_dir
    wdir.mkdir(exist_ok=True)
    rdir = Path(f'{wdir}/result_files')
    
    # Resolve FF paths
    ff_specifiers = []
    for ff_file in args.ff_files:
        ff_path = Path(ff_file)
        if ff_path.exists():
            ff_specifiers.append(ff_path)
        else:
            raise FileNotFoundError(f"Missing force field file: {ff_path}")

    if args.extra_ff is not None:
        if not args.extra_ff.exists():
            raise FileNotFoundError(f"Missing extra force field file: {args.extra_ff}")
        ff_specifiers.append(args.extra_ff)

    for ff_path in ff_specifiers:
        if not ff_path.exists():
            raise FileNotFoundError(f"Missing force field file: {ff_path}")

    ff = ForceField(*ff_specifiers)

    # Helpers
    z_registry = ReplicateZRegistry()
    mb_manager = MembraneBarostatManager(
        pressure_bar=args.pressure_bar,
        temperature_K=args.temperature_K,
        surface_tension_bar_nm=args.surface_tension_bar_nm,
        frequency=args.barostat_freq
    )

    # Ensure each concentration key uses the same Lz target you requested
    # (set once before we touch any System/Simulation)
    # conc_key format: 'X.Ym' (e.g., '2.0m')
    lz_nm_value = args.lz.value_in_unit(nanometer)

    omm_modobj_all = {}

    for mol in concentration_list:
        # ---- consistent molality strings ----
        if mol % 1 == 0:
            mi1 = f"{mol:.1f}"
            mi = str(int(mol))
        else:
            mi  = f"{mol:.1f}".replace('.', '')
            mi1 = f"{mol:.1f}"

        conc_key = f'{mi1}m'
        z_registry.set_target_for_key(conc_key, lz_nm_value)

        # Build/load systems for this concentration (independent of k)
        placeholder_k = args.k_quantity

        if args.restraint_type == 'FBP':
            builder = build_system_FBP
        else:
            builder = build_system_HP

        omm_modobj = builder(
            concentration=mol, 
            repnum=repnumber, 
            ion1=ion1,
            ion2=ion2,
            wdir=wdir,
            ff=ff,
            water=water,
            salt_dict=salt_dict,
            center_atom1=args.centerion[0],
            center_atom2=args.centerion[1],
            k=placeholder_k,                 # placeholder; we reset 'k' later if needed
            delta_z=args.delta_z, 
            z_center=args.z_center,          # center point used by modify_*
            lz_total=args.lz,                 # total Z length
            expected_pairs=args.row.num_particles,
            pdb_suffix=pdb_suffix,
        )
        # Store by conc_key for consistent access later
        omm_modobj_all[conc_key] = omm_modobj

        # Now loop over k values
        for r in range(repnumber):
            postfix = f'{mi}m_{args.ktag}_r{r}'

            if simulation_exists(wdir, postfix):
                LOGGER.info(f"Skipping {postfix}, simulation folders already exist.")
                continue

            try:
                LOGGER.info(f"▶ Starting simulation {postfix}")

                # ---- pull base objects ----
                omm_top = omm_modobj_all[conc_key][f'r{r}']['topology']
                omm_sys = omm_modobj_all[conc_key][f'r{r}']['system']
                omm_pos = omm_modobj_all[conc_key][f'r{r}']['positions']

                # ---- update restraint parameters ----
                for i_force in range(omm_sys.getNumForces()):
                    f = omm_sys.getForce(i_force)
                    if isinstance(f, CustomExternalForce):
                        for p in range(f.getNumGlobalParameters()):
                            pname = f.getGlobalParameterName(p)
                            if pname == 'k':
                                f.setGlobalParameterDefaultValue(p, args.k_quantity)
                            elif pname == 'z0':
                                f.setGlobalParameterDefaultValue(p, args.z_center)

                # ---- enforce Z ----
                z_registry.enforce_on_system(omm_sys, conc_key)

                # ---------- Simulation schedules ----------

                schedule1: dict[str, SimulationParameters] = {
                    f"equil_sim_NVT_{postfix}": SimulationParameters(
                        integ_params=IntegratorParameters(
                            time_step=1.0 * femtoseconds,
                            total_time=args.equil_nvt_time,
                            num_samples=args.nvt_samples,
                        ),
                        thermo_params=ThermoParameters(
                            thermostat_params=ThermostatParameters(
                            temperature=args.temperature,
                            )
                        ),
                        reporter_params=ReporterParameters(traj_ext="dcd"),
                    )
                }

                schedule2: dict[str, SimulationParameters] = {
                    f"equil_sim_NPT_{postfix}": SimulationParameters(
                        integ_params=IntegratorParameters(
                            time_step=args.time_step,
                            total_time=args.equil_npt_time,
                            num_samples=args.nvt_samples,
                        ),
                        thermo_params=ThermoParameters(
                            thermostat_params=ThermostatParameters(
                            temperature=args.temperature,
                            )
                        ),
                        reporter_params=ReporterParameters(traj_ext="dcd"),
                    ),
                    f"prod_sim_{postfix}": SimulationParameters(
                        integ_params=IntegratorParameters(
                            time_step=args.time_step,
                            total_time=args.prod_time,
                            num_samples=args.prod_samples,
                        ),
                        thermo_params=ThermoParameters(
                            thermostat_params=ThermostatParameters(
                            temperature=args.temperature,
                            )
                        ),
                        reporter_params=ReporterParameters(),
                    ),
                }

                # ---------- NVT ----------
                history1 = run_simulation_schedule(
                    working_dir=wdir,
                    schedule=schedule1,
                    init_top=omm_top,
                    init_sys=omm_sys,
                    init_pos=omm_pos,
                    return_history=True
                )

                nvt_sim = history1[f'equil_sim_NVT_{postfix}']['simulation']
                nvt_sys = nvt_sim.system

                mb_manager.add_to_simulation(nvt_sim)
                z_registry.enforce_on_simulation(nvt_sim, conc_key)

                logging.info(nvt_sys.getForces())
                nvt_top = nvt_sim.topology
                nvt_pos = get_context_positions(nvt_sim.context)

                
                nvt_top = nvt_sim.topology
                nvt_pos = get_context_positions(nvt_sim.context) 

                z_registry.enforce_on_system(nvt_sys, conc_key)

                # ---------- NPT + PROD ----------
                history2 = run_simulation_schedule(
                    working_dir=wdir,
                    schedule=schedule2,
                    init_top=nvt_top,
                    init_sys=nvt_sys,
                    init_pos=nvt_pos,
                    return_history=True
                )

                LOGGER.info(f"✅ Completed {postfix}")

            except openmm.OpenMMException as e:
                LOGGER.error(
                    f"❌ OpenMM failure (NaN likely) — skipped replicate | "
                    f"molality={mi1}m, k={args.ktag}, r={r}\n{e}"
                )
                continue

            except Exception as e:
                LOGGER.exception(
                    f"🔥 Unexpected error — skipped replicate | "
                    f"molality={mi1}m, k={args.ktag}, r={r}"
                )
                continue

    # --------- Collect outputs for all concentrations & ktags ---------
    copy_simulation_outputs(
        molalities=concentration_list,
        N_replicates=repnumber,
        wdir=wdir,
        rdir=rdir,
        ktags=[args.ktag],
        topology_from='pdb',
        stride=None,
        force_convert=True
    )


if __name__ == '__main__':
    main()