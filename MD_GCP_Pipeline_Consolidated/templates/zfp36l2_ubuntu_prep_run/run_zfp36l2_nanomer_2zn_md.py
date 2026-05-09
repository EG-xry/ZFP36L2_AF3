#!/usr/bin/env python3
"""
MD simulation script for ZFP36L2 (WT) + RNA 9-mer + 2 Zn ions.

Structure: AF3 best model PDB for WT + RNA(9-mer) + 2 Zn
Entities: 1 protein chain (A), 1 RNA chain (D), 2 Zn ions (B, C)

Zn Coordination Sites (HARD-CODED):
- ZF1: Zn1 (chain B) - Cys9 SG, Cys18 SG, Cys24 SG, His28 NE2
- ZF2: Zn2 (chain C) - Cys47 SG, Cys56 SG, Cys62 SG, His66 NE2

Force Field: AMBER14 protein + RNA OL3 + TIP3P water
Zn distance restraints: RAMPED from weak to strong to prevent explosions

Workflow:
A) Energy minimization with WEAK Zn restraints (100,000 max iterations)
B) NVT equilibration (100 ps) with MEDIUM Zn restraints
C) NPT equilibration (500 ps) with MEDIUM Zn restraints
D) Restrained NPT settle (1 ns) with FINAL Zn restraints
E) Production NPT (configurable) with FINAL Zn restraints

Output:
- Single log file with stage banners
- Zn distance log with stage column
- Trajectory: DCD every 10 ps
- Checkpoint: every 5 minutes
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from openmm import (
    CustomBondForce,
    CustomExternalForce,
    LangevinMiddleIntegrator,
    MonteCarloBarostat,
    Platform,
    unit,
)
from openmm.app import (
    CheckpointReporter,
    DCDReporter,
    Element,
    ForceField,
    HBonds,
    Modeller,
    PDBFile,
    PME,
    Simulation,
    StateDataReporter,
)


# ============================================================================
# Configuration
# ============================================================================

BASE_DIR = Path(__file__).resolve().parent
CLEANED_DIR = BASE_DIR / "cleaned_pdbs"
DEFAULT_PDB = CLEANED_DIR / "fold_zfp36l2_wt_nanomer_9_wt_2zn_model_0_fixed.pdb"
OUTPUT_ROOT = BASE_DIR / "zfp36l2_nanomer_2zn_md"

# HARD-CODED Zn coordination sites - residue numbers and atom names
# ZF1 Zn ligands: Cys9 SG, Cys18 SG, Cys24 SG, His28 NE2
ZF1_LIGANDS = [
    ("CYS", 9, "SG", 0.230),    # Zn-Cys: 2.30 Å = 0.230 nm
    ("CYS", 18, "SG", 0.230),
    ("CYS", 24, "SG", 0.230),
    ("HIS", 28, "NE2", 0.210),  # Zn-His: 2.10 Å = 0.210 nm
]

# ZF2 Zn ligands: Cys47 SG, Cys56 SG, Cys62 SG, His66 NE2
ZF2_LIGANDS = [
    ("CYS", 47, "SG", 0.230),
    ("CYS", 56, "SG", 0.230),
    ("CYS", 62, "SG", 0.230),
    ("HIS", 66, "NE2", 0.210),
]

# Zn restraint ramping schedule (kJ/mol/nm²)
ZN_K_WEAK = 100.0       # For minimization - gentle pull
ZN_K_MEDIUM = 500.0     # For early equilibration
ZN_K_FINAL = 1000.0     # For late equilibration and production

RNA_BASES = {"A", "C", "G", "U"}


# ============================================================================
# Utility Functions
# ============================================================================

def time_to_steps(duration_ns: float, timestep_fs: float) -> int:
    """Convert nanoseconds to steps."""
    total_ps = duration_ns * 1000.0
    step_ps = timestep_fs / 1000.0
    return int(total_ps / step_ps)


def ps_to_steps(interval_ps: float, timestep_fs: float) -> int:
    """Convert picoseconds to steps."""
    step_ps = timestep_fs / 1000.0
    return max(1, int(interval_ps / step_ps))


def create_forcefield() -> ForceField:
    """Load Amber14 protein + RNA OL3 + TIP3P water forcefields."""
    combos = [
        (
            "amber14/protein.ff14SB.xml",
            "amber14/RNA.OL3.xml",
            "amber14/tip3p.xml",
        ),
        ("amber14-all.xml", "amber14/tip3p.xml"),
        ("amber14-all.xml", "tip3p.xml"),
    ]

    last_error: Optional[ValueError] = None
    for files in combos:
        try:
            print(f"[ForceField] Trying: {files}")
            ff = ForceField(*files)
            print(f"[ForceField] Successfully loaded: {files}")
            return ff
        except ValueError as exc:
            last_error = exc
            print(f"[ForceField] Failed: {exc}")
    raise last_error


def _canonical_rna_base(name: str) -> Optional[str]:
    """Extract canonical RNA base letter."""
    upper = name.strip().upper()
    if not upper:
        return None
    letter = upper[-1]
    return letter if letter in RNA_BASES else None


def identify_terminal_template_ids(topology) -> Dict[Tuple[str, str], str]:
    """Determine which residues should use 5'/3' RNA templates."""
    template_ids: Dict[Tuple[str, str], str] = {}
    for chain in topology.chains():
        residues = list(chain.residues())
        if not residues:
            continue
        rna_residues = [res for res in residues if _canonical_rna_base(res.name)]
        if len(rna_residues) == 0 or len(rna_residues) != len(residues):
            continue
        first = rna_residues[0]
        last = rna_residues[-1]
        first_base = _canonical_rna_base(first.name)
        last_base = _canonical_rna_base(last.name)
        if first_base:
            template_ids[(chain.id, first.id)] = f"{first_base}5"
        if last_base:
            template_ids[(chain.id, last.id)] = f"{last_base}3"
    return template_ids


def build_residue_template_map(topology, template_ids: Dict[Tuple[str, str], str]):
    """Build residue template mapping."""
    mapping = {}
    for residue in topology.residues():
        key = (residue.chain.id, residue.id)
        template = template_ids.get(key)
        if template:
            mapping[residue] = template
    return mapping


def write_banner(log_path: Path, stage_name: str, append: bool = True):
    """Write a banner line to the log file."""
    mode = "a" if append else "w"
    with open(log_path, mode) as f:
        f.write(f"\n{'='*60}\n")
        f.write(f"===== STAGE: {stage_name} =====\n")
        f.write(f"{'='*60}\n\n")


# ============================================================================
# Position Restraints
# ============================================================================

def add_position_restraints(
    system,
    topology,
    reference_positions,
    force_constant,
    heavy_only: bool = True,
    name: str = "PositionRestraints"
) -> int:
    """Add position restraints to the system. Returns the force index."""
    force = CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addGlobalParameter("k", force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")
    force.setName(name)

    positions_nm = reference_positions.value_in_unit(unit.nanometer)
    n_restrained = 0
    for atom in topology.atoms():
        if heavy_only and atom.element is not None and atom.element.symbol == "H":
            continue
        pos = positions_nm[atom.index]
        force.addParticle(atom.index, pos)
        n_restrained += 1

    force_idx = system.addForce(force)
    print(f"[Restraints] Added {n_restrained} position restraints (k={force_constant})")
    return force_idx


def update_restraint_force_constant(context, force_constant, param_name: str = "k"):
    """Update the force constant for position restraints."""
    context.setParameter(param_name, force_constant)
    print(f"[Restraints] Updated {param_name} to {force_constant}")


# ============================================================================
# Zn Coordination Restraints (HARD-CODED)
# ============================================================================

def find_atom_index(topology, chain_id: str, res_name: str, res_num: int, atom_name: str) -> Optional[int]:
    """Find atom index by chain, residue name/number, and atom name."""
    for chain in topology.chains():
        if chain.id != chain_id:
            continue
        for residue in chain.residues():
            res_idx = int(residue.id) if residue.id.isdigit() else None
            if res_idx != res_num:
                continue
            if residue.name != res_name:
                continue
            for atom in residue.atoms():
                if atom.name == atom_name:
                    return atom.index
    return None


def find_zn_indices(topology) -> List[Tuple[int, str]]:
    """Find all Zn atom indices with their residue IDs."""
    zn_indices = []
    for atom in topology.atoms():
        if atom.element == Element.getBySymbol("Zn"):
            zn_indices.append((atom.index, atom.residue.id))
    return zn_indices


def assign_zn_to_ligands_hardcoded(
    topology,
    positions,
    zf1_ligands: List[Tuple[str, int, str, float]],
    zf2_ligands: List[Tuple[str, int, str, float]],
) -> List[Tuple[int, int, float, str]]:
    """
    HARD-CODE the Zn-ligand assignments.
    
    We find which Zn is closest to which zinc finger and assign ALL 4 ligands
    per Zn, regardless of current distances.
    
    Returns list of (zn_idx, ligand_idx, r0_nm, description).
    """
    zn_data = find_zn_indices(topology)
    if len(zn_data) != 2:
        raise ValueError(f"Expected 2 Zn ions, found {len(zn_data)}")
    
    positions_nm = positions.value_in_unit(unit.nanometer)
    
    # Find all ligand atom indices first
    zf1_atom_indices = []
    for res_name, res_num, atom_name, r0 in zf1_ligands:
        idx = find_atom_index(topology, "A", res_name, res_num, atom_name)
        if idx is None:
            raise ValueError(f"Could not find ZF1 ligand: {res_name}{res_num}:{atom_name}")
        zf1_atom_indices.append((idx, res_name, res_num, atom_name, r0))
    
    zf2_atom_indices = []
    for res_name, res_num, atom_name, r0 in zf2_ligands:
        idx = find_atom_index(topology, "A", res_name, res_num, atom_name)
        if idx is None:
            raise ValueError(f"Could not find ZF2 ligand: {res_name}{res_num}:{atom_name}")
        zf2_atom_indices.append((idx, res_name, res_num, atom_name, r0))
    
    # Calculate center of mass for each zinc finger ligand set
    zf1_center = np.mean([positions_nm[idx] for idx, *_ in zf1_atom_indices], axis=0)
    zf2_center = np.mean([positions_nm[idx] for idx, *_ in zf2_atom_indices], axis=0)
    
    # Assign each Zn to the closer zinc finger
    zn_assignments = {}  # zn_idx -> "ZF1" or "ZF2"
    
    for zn_idx, zn_res_id in zn_data:
        zn_pos = np.array(positions_nm[zn_idx])
        dist_to_zf1 = np.linalg.norm(zn_pos - zf1_center)
        dist_to_zf2 = np.linalg.norm(zn_pos - zf2_center)
        
        if dist_to_zf1 < dist_to_zf2:
            zn_assignments[zn_idx] = "ZF1"
            print(f"[Zn] Zn atom {zn_idx} (res {zn_res_id}) → ZF1 (dist: {dist_to_zf1:.3f} nm)")
        else:
            zn_assignments[zn_idx] = "ZF2"
            print(f"[Zn] Zn atom {zn_idx} (res {zn_res_id}) → ZF2 (dist: {dist_to_zf2:.3f} nm)")
    
    # Check that we have one Zn for each finger
    assigned_fingers = set(zn_assignments.values())
    if assigned_fingers != {"ZF1", "ZF2"}:
        # Force assignment if both Zn went to same finger
        zn_list = list(zn_assignments.keys())
        print(f"[Zn] WARNING: Both Zn assigned to same finger, forcing separation")
        zn_assignments[zn_list[0]] = "ZF1"
        zn_assignments[zn_list[1]] = "ZF2"
    
    # Build ALL 8 restraint pairs (4 per Zn)
    pairs = []
    for zn_idx, finger in zn_assignments.items():
        if finger == "ZF1":
            ligands = zf1_atom_indices
        else:
            ligands = zf2_atom_indices
        
        for lig_idx, res_name, res_num, atom_name, r0 in ligands:
            desc = f"{finger}:Zn{zn_idx}-{res_name}{res_num}:{atom_name}"
            current_dist = np.linalg.norm(np.array(positions_nm[zn_idx]) - np.array(positions_nm[lig_idx]))
            print(f"[Zn] Restraint: {desc}, target={r0:.3f} nm, current={current_dist:.3f} nm")
            pairs.append((zn_idx, lig_idx, r0, desc))
    
    if len(pairs) != 8:
        raise ValueError(f"Expected 8 Zn-ligand restraints, got {len(pairs)}")
    
    return pairs


def add_zn_distance_restraints(
    system,
    topology,
    positions,
    force_constant: float,
    zn_ligand_pairs: List[Tuple[int, int, float]],  # (zn_idx, ligand_idx, r0_nm)
    force_name: str = "ZnDistanceRestraints",
) -> int:
    """
    Add Zn-ligand distance restraints using a CustomBondForce.
    
    Uses a separate parameter name so we can update Zn restraints independently.
    Returns force_index.
    """
    zn_force = CustomBondForce("0.5*k_zn*(r-r0)^2")
    zn_force.addGlobalParameter("k_zn", force_constant)
    zn_force.addPerBondParameter("r0")
    zn_force.setName(force_name)

    bonds_added = 0
    for zn_idx, ligand_idx, r0_nm in zn_ligand_pairs:
        zn_force.addBond(zn_idx, ligand_idx, [r0_nm])
        bonds_added += 1

    force_idx = system.addForce(zn_force)
    print(f"[Zn] Added {bonds_added} Zn-ligand restraints (k={force_constant} kJ/mol/nm²)")
    return force_idx


def compute_zn_ligand_distances(context, zn_ligand_pairs: List[Tuple[int, int, float, str]]) -> Dict[str, float]:
    """Compute current Zn-ligand distances."""
    state = context.getState(getPositions=True)
    positions_nm = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    
    distances = {}
    for zn_idx, lig_idx, r0, desc in zn_ligand_pairs:
        zn_pos = positions_nm[zn_idx]
        lig_pos = positions_nm[lig_idx]
        dist = np.linalg.norm(zn_pos - lig_pos)
        distances[desc] = dist
    return distances


def log_zn_distances(
    context, 
    zn_ligand_pairs: List[Tuple[int, int, float, str]], 
    zn_log_path: Path,
    stage_name: str,
    step: int,
    append: bool = True,
    write_header: bool = False,
):
    """Log Zn-ligand distances with stage column."""
    mode = "a" if append else "w"
    with open(zn_log_path, mode) as f:
        if write_header:
            header = "stage\tstep\t" + "\t".join(desc for _, _, _, desc in zn_ligand_pairs)
            f.write(header + "\n")
        
        dists = compute_zn_ligand_distances(context, zn_ligand_pairs)
        line = f"{stage_name}\t{step}\t" + "\t".join(f"{dists[desc]:.4f}" for _, _, _, desc in zn_ligand_pairs)
        f.write(line + "\n")


# ============================================================================
# System Setup
# ============================================================================

def build_modeller(
    pdb_path: Path,
    forcefield: ForceField,
    padding_nm: float,
    ionic_strength: float,
    neutralize: bool = True,
) -> Tuple[Modeller, Dict[Tuple[str, str], str]]:
    """Load PDB, add hydrogens, solvate, and add ions."""
    print(f"[Setup] Loading PDB: {pdb_path}")
    pdb = PDBFile(str(pdb_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    
    # Identify RNA terminals
    template_ids = identify_terminal_template_ids(modeller.topology)
    residue_templates = build_residue_template_map(modeller.topology, template_ids)
    
    # Count entities
    n_protein = sum(1 for r in modeller.topology.residues() if r.name not in RNA_BASES and r.name != "ZN")
    n_rna = sum(1 for r in modeller.topology.residues() if r.name in RNA_BASES)
    n_zn = len([a for a in modeller.topology.atoms() if a.element == Element.getBySymbol("Zn")])
    print(f"[Setup] Entities: {n_protein} protein residues, {n_rna} RNA residues, {n_zn} Zn ions")
    
    # Add hydrogens
    print("[Setup] Adding hydrogens (pH 7.0)...")
    modeller.addHydrogens(forcefield, pH=7.0, residueTemplates=residue_templates)
    
    # Update residue templates
    residue_templates = build_residue_template_map(modeller.topology, template_ids)
    
    # Solvate
    print(f"[Setup] Adding solvent (TIP3P, padding={padding_nm} nm, ions={ionic_strength} M)...")
    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=padding_nm * unit.nanometer,
        ionicStrength=ionic_strength * unit.molar,
        neutralize=neutralize,
        residueTemplates=residue_templates,
    )
    
    n_atoms = modeller.topology.getNumAtoms()
    n_residues = modeller.topology.getNumResidues()
    print(f"[Setup] System: {n_atoms} atoms, {n_residues} residues")
    
    return modeller, template_ids


def create_system(
    forcefield: ForceField,
    topology,
    residue_templates: Optional[Dict],
    with_barostat: bool,
    temperature,
    pressure,
    barostat_interval: int = 25,
    nonbonded_cutoff_nm: float = 1.2,
):
    """Create OpenMM System with specified parameters."""
    print(f"[System] Creating system (PME, cutoff={nonbonded_cutoff_nm} nm, barostat={with_barostat})...")
    
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=PME,
        nonbondedCutoff=nonbonded_cutoff_nm * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
        removeCMMotion=True,
        residueTemplates=residue_templates or {},
    )
    
    if with_barostat:
        barostat = MonteCarloBarostat(pressure, temperature, barostat_interval)
        system.addForce(barostat)
        print(f"[System] Added MonteCarloBarostat (P={pressure}, T={temperature}, interval={barostat_interval})")
    
    return system


def get_platform(name: Optional[str]):
    """Get specified platform or None for automatic selection."""
    if name is None:
        return None
    return Platform.getPlatformByName(name)


# ============================================================================
# Simulation Stages
# ============================================================================

def run_minimization(
    simulation: Simulation, 
    max_iterations: int,
    stage_name: str,
    log_path: Path,
    append_log: bool,
    zn_pairs: Optional[List] = None,
    zn_log_path: Optional[Path] = None,
):
    """Run energy minimization with logging."""
    write_banner(log_path, stage_name, append=append_log)
    
    print(f"[{stage_name}] Starting energy minimization (max {max_iterations} iterations)...")
    start_time = time.time()
    
    state_before = simulation.context.getState(getEnergy=True)
    e_before = state_before.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    
    # Log initial Zn distances
    if zn_pairs and zn_log_path:
        log_zn_distances(
            simulation.context, zn_pairs, zn_log_path, 
            stage_name, 0, append=append_log, write_header=not append_log
        )
    
    simulation.minimizeEnergy(maxIterations=max_iterations)
    
    state_after = simulation.context.getState(getEnergy=True)
    e_after = state_after.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    
    elapsed = time.time() - start_time
    print(f"[{stage_name}] Done in {elapsed:.1f}s: E = {e_before:.1f} → {e_after:.1f} kJ/mol")
    
    # Log final Zn distances
    if zn_pairs and zn_log_path:
        log_zn_distances(
            simulation.context, zn_pairs, zn_log_path,
            stage_name, max_iterations, append=True
        )
    
    # Write minimization result to log
    with open(log_path, "a") as f:
        f.write(f"# Minimization: {max_iterations} iterations, E: {e_before:.1f} → {e_after:.1f} kJ/mol\n")
        f.write(f"# Time: {elapsed:.1f}s\n\n")
    
    return elapsed


def run_stage(
    stage_name: str,
    simulation: Simulation,
    steps: int,
    log_path: Path,
    log_interval_steps: int,
    trajectory_reporter: Optional[DCDReporter] = None,
    checkpoint_path: Optional[Path] = None,
    checkpoint_interval_steps: Optional[int] = None,
    zn_ligand_pairs: Optional[List] = None,
    zn_log_path: Optional[Path] = None,
    zn_log_interval_steps: Optional[int] = None,
):
    """Run a simulation stage with stage banner and reporters."""
    
    # Write banner
    write_banner(log_path, stage_name, append=True)
    
    # State data reporter (append mode)
    reporters = [
        StateDataReporter(
            str(log_path),
            log_interval_steps,
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            volume=True,
            density=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=steps,
            separator="\t",
            append=True,  # Always append since we write banner first
        )
    ]
    
    if trajectory_reporter is not None:
        reporters.append(trajectory_reporter)
    
    if checkpoint_path is not None and checkpoint_interval_steps is not None:
        reporters.append(CheckpointReporter(str(checkpoint_path), checkpoint_interval_steps))
    
    simulation.reporters.extend(reporters)
    
    print(f"[{stage_name}] Running {steps} steps...")
    start_time = time.time()
    
    # Run with periodic Zn distance logging
    if zn_ligand_pairs and zn_log_path and zn_log_interval_steps:
        steps_done = 0
        while steps_done < steps:
            chunk = min(zn_log_interval_steps, steps - steps_done)
            simulation.step(chunk)
            steps_done += chunk
            
            # Log Zn distances with stage name
            current_step = simulation.currentStep
            log_zn_distances(
                simulation.context, zn_ligand_pairs, zn_log_path,
                stage_name, current_step, append=True
            )
    else:
        simulation.step(steps)
    
    elapsed = time.time() - start_time
    ns_per_day = (steps * 0.002 / 1000) / (elapsed / 86400) if elapsed > 0 else 0
    print(f"[{stage_name}] Done in {elapsed:.1f}s ({ns_per_day:.1f} ns/day)")
    
    simulation.reporters.clear()
    return elapsed


def write_pdb(topology, positions, path: Path):
    """Write positions to PDB file."""
    with path.open("w") as f:
        PDBFile.writeFile(topology, positions, f, keepIds=True)
    print(f"[Output] Wrote {path}")


# ============================================================================
# Main Pipeline
# ============================================================================

def run_pipeline(
    pdb_path: Path,
    output_root: Path,
    forcefield: ForceField,
    platform_name: Optional[str],
    # Solvation
    padding_nm: float = 1.0,
    ionic_strength: float = 0.15,
    # System
    nonbonded_cutoff_nm: float = 1.2,
    temperature_k: float = 300.0,
    pressure_bar: float = 1.0,
    friction_per_ps: float = 1.0,
    timestep_fs: float = 2.0,
    barostat_interval: int = 25,
    # Minimization
    minimize_max_iters: int = 100000,
    # Equilibration
    nvt_ps: float = 100.0,
    npt_eq_ps: float = 500.0,
    npt_settle_ns: float = 1.0,
    heavy_restraint_k: float = 1000.0,
    heavy_restraint_k_reduced: float = 200.0,
    # Production
    production_ns: float = 5.0,
    # Zn restraints (ramped)
    zn_k_weak: float = ZN_K_WEAK,
    zn_k_medium: float = ZN_K_MEDIUM,
    zn_k_final: float = ZN_K_FINAL,
    # Output
    log_interval_ps: float = 1.0,
    traj_interval_ps: float = 10.0,
    checkpoint_interval_min: float = 5.0,
    zn_log_interval_ps: float = 10.0,
):
    """Run the complete MD pipeline with Zn restraint ramping."""
    
    start_time = time.time()
    
    # Setup output directory
    pdb_stem = pdb_path.stem
    out_dir = output_root / pdb_stem
    out_dir.mkdir(parents=True, exist_ok=True)
    
    log_path = out_dir / "md.log"
    traj_path = out_dir / "production.dcd"
    checkpoint_path = out_dir / "checkpoint.chk"
    final_pdb_path = out_dir / "final_frame.pdb"
    minimized_pdb_path = out_dir / "minimized.pdb"
    equilibrated_pdb_path = out_dir / "equilibrated.pdb"
    zn_dist_log_path = out_dir / "zn_distances.tsv"
    
    print("=" * 70)
    print(f"ZFP36L2 + RNA 9-mer + 2 Zn MD Simulation")
    print(f"Input: {pdb_path}")
    print(f"Output: {out_dir}")
    print(f"Zn restraint ramping: weak={zn_k_weak} → medium={zn_k_medium} → final={zn_k_final} kJ/mol/nm²")
    print("=" * 70)
    
    # Build system
    modeller, template_ids = build_modeller(
        pdb_path, forcefield, padding_nm, ionic_strength
    )
    topology = modeller.topology
    positions = modeller.positions
    residue_templates = build_residue_template_map(topology, template_ids)
    
    # Find Zn-ligand pairs for restraints (HARD-CODED assignments)
    print("\n[Zn] Setting up HARD-CODED Zn-ligand restraints...")
    zn_pairs = assign_zn_to_ligands_hardcoded(topology, positions, ZF1_LIGANDS, ZF2_LIGANDS)
    zn_pair_tuples = [(zn, lig, r0) for zn, lig, r0, _ in zn_pairs]
    
    # Platform and integrator setup
    platform = get_platform(platform_name)
    timestep = timestep_fs * unit.femtoseconds
    temperature = temperature_k * unit.kelvin
    pressure = pressure_bar * unit.bar
    friction = friction_per_ps / unit.picoseconds
    
    def new_integrator():
        return LangevinMiddleIntegrator(temperature, friction, timestep)
    
    # Convert time intervals to steps
    nvt_steps = ps_to_steps(nvt_ps, timestep_fs)
    npt_eq_steps = ps_to_steps(npt_eq_ps, timestep_fs)
    npt_settle_steps = time_to_steps(npt_settle_ns, timestep_fs)
    production_steps = time_to_steps(production_ns, timestep_fs)
    log_interval_steps = ps_to_steps(log_interval_ps, timestep_fs)
    traj_interval_steps = ps_to_steps(traj_interval_ps, timestep_fs)
    checkpoint_interval_steps = ps_to_steps(checkpoint_interval_min * 60 * 1000, timestep_fs)
    zn_log_steps = ps_to_steps(zn_log_interval_ps, timestep_fs)
    
    heavy_k = heavy_restraint_k * unit.kilojoule_per_mole / unit.nanometer**2
    heavy_k_reduced = heavy_restraint_k_reduced * unit.kilojoule_per_mole / unit.nanometer**2
    
    # ========================================================================
    # Stage 1: Energy Minimization with WEAK Zn restraints
    # ========================================================================
    print("\n" + "=" * 70)
    print(f"Stage 1: Energy Minimization (WEAK Zn k={zn_k_weak})")
    print("=" * 70)
    
    system_min = create_system(
        forcefield, topology, residue_templates,
        with_barostat=False, temperature=temperature, pressure=pressure,
        nonbonded_cutoff_nm=nonbonded_cutoff_nm,
    )
    add_position_restraints(system_min, topology, positions, heavy_k, heavy_only=True)
    add_zn_distance_restraints(system_min, topology, positions, zn_k_weak, zn_pair_tuples)
    
    integrator_min = new_integrator()
    sim_min = Simulation(topology, system_min, integrator_min, platform)
    sim_min.context.setPositions(positions)
    
    run_minimization(
        sim_min, minimize_max_iters,
        f"Minimize (Zn k={zn_k_weak})",
        log_path, append_log=False,
        zn_pairs=zn_pairs,
        zn_log_path=zn_dist_log_path,
    )
    state_min = sim_min.context.getState(getPositions=True)
    pos_min = state_min.getPositions()
    write_pdb(topology, pos_min, minimized_pdb_path)
    
    # ========================================================================
    # Stage 2: NVT Equilibration with MEDIUM Zn restraints
    # ========================================================================
    print("\n" + "=" * 70)
    print(f"Stage 2: NVT Equilibration ({nvt_ps} ps, MEDIUM Zn k={zn_k_medium})")
    print("=" * 70)
    
    system_nvt = create_system(
        forcefield, topology, residue_templates,
        with_barostat=False, temperature=temperature, pressure=pressure,
        nonbonded_cutoff_nm=nonbonded_cutoff_nm,
    )
    add_position_restraints(system_nvt, topology, pos_min, heavy_k, heavy_only=True)
    add_zn_distance_restraints(system_nvt, topology, pos_min, zn_k_medium, zn_pair_tuples)
    
    integrator_nvt = new_integrator()
    sim_nvt = Simulation(topology, system_nvt, integrator_nvt, platform)
    sim_nvt.context.setPositions(pos_min)
    sim_nvt.context.setVelocitiesToTemperature(temperature)
    
    run_stage(
        f"NVT ({nvt_ps} ps)",
        sim_nvt, nvt_steps, log_path, log_interval_steps,
        zn_ligand_pairs=zn_pairs,
        zn_log_path=zn_dist_log_path,
        zn_log_interval_steps=zn_log_steps,
    )
    state_nvt = sim_nvt.context.getState(getPositions=True)
    pos_nvt = state_nvt.getPositions()
    
    # ========================================================================
    # Stage 3: NPT Equilibration with MEDIUM Zn restraints
    # ========================================================================
    print("\n" + "=" * 70)
    print(f"Stage 3: NPT Equilibration ({npt_eq_ps} ps, MEDIUM Zn k={zn_k_medium})")
    print("=" * 70)
    
    system_npt_eq = create_system(
        forcefield, topology, residue_templates,
        with_barostat=True, temperature=temperature, pressure=pressure,
        barostat_interval=barostat_interval,
        nonbonded_cutoff_nm=nonbonded_cutoff_nm,
    )
    add_position_restraints(system_npt_eq, topology, pos_nvt, heavy_k, heavy_only=True)
    add_zn_distance_restraints(system_npt_eq, topology, pos_nvt, zn_k_medium, zn_pair_tuples)
    
    integrator_npt_eq = new_integrator()
    sim_npt_eq = Simulation(topology, system_npt_eq, integrator_npt_eq, platform)
    sim_npt_eq.context.setPositions(pos_nvt)
    sim_npt_eq.context.setVelocitiesToTemperature(temperature)
    
    run_stage(
        f"NPT-eq ({npt_eq_ps} ps)",
        sim_npt_eq, npt_eq_steps, log_path, log_interval_steps,
        zn_ligand_pairs=zn_pairs,
        zn_log_path=zn_dist_log_path,
        zn_log_interval_steps=zn_log_steps,
    )
    state_npt_eq = sim_npt_eq.context.getState(getPositions=True)
    pos_npt_eq = state_npt_eq.getPositions()
    
    # ========================================================================
    # Stage 4: Restrained NPT Settle with FINAL Zn restraints
    # ========================================================================
    print("\n" + "=" * 70)
    print(f"Stage 4: Restrained NPT Settle ({npt_settle_ns} ns, FINAL Zn k={zn_k_final})")
    print("=" * 70)
    
    system_settle = create_system(
        forcefield, topology, residue_templates,
        with_barostat=True, temperature=temperature, pressure=pressure,
        barostat_interval=barostat_interval,
        nonbonded_cutoff_nm=nonbonded_cutoff_nm,
    )
    add_position_restraints(system_settle, topology, pos_npt_eq, heavy_k, heavy_only=True)
    add_zn_distance_restraints(system_settle, topology, pos_npt_eq, zn_k_final, zn_pair_tuples)
    
    integrator_settle = new_integrator()
    sim_settle = Simulation(topology, system_settle, integrator_settle, platform)
    sim_settle.context.setPositions(pos_npt_eq)
    sim_settle.context.setVelocitiesToTemperature(temperature)
    
    # First half with full heavy-atom restraints
    half_steps = npt_settle_steps // 2
    run_stage(
        f"NPT-settle-1 ({npt_settle_ns/2} ns, heavy-k={heavy_restraint_k})",
        sim_settle, half_steps, log_path, log_interval_steps,
        zn_ligand_pairs=zn_pairs,
        zn_log_path=zn_dist_log_path,
        zn_log_interval_steps=zn_log_steps,
    )
    
    # Reduce heavy-atom restraints for second half
    update_restraint_force_constant(sim_settle.context, heavy_k_reduced, "k")
    
    run_stage(
        f"NPT-settle-2 ({npt_settle_ns/2} ns, heavy-k={heavy_restraint_k_reduced})",
        sim_settle, npt_settle_steps - half_steps, log_path, log_interval_steps,
        zn_ligand_pairs=zn_pairs,
        zn_log_path=zn_dist_log_path,
        zn_log_interval_steps=zn_log_steps,
    )
    
    state_settle = sim_settle.context.getState(getPositions=True)
    pos_settle = state_settle.getPositions()
    write_pdb(topology, pos_settle, equilibrated_pdb_path)
    
    # ========================================================================
    # Stage 5: Production NPT (no heavy-atom restraints, FINAL Zn restraints)
    # ========================================================================
    print("\n" + "=" * 70)
    print(f"Stage 5: Production NPT ({production_ns} ns, Zn k={zn_k_final}, no heavy-atom restraints)")
    print("=" * 70)
    
    system_prod = create_system(
        forcefield, topology, residue_templates,
        with_barostat=True, temperature=temperature, pressure=pressure,
        barostat_interval=barostat_interval,
        nonbonded_cutoff_nm=nonbonded_cutoff_nm,
    )
    # Zn restraints ON, heavy-atom restraints OFF
    add_zn_distance_restraints(system_prod, topology, pos_settle, zn_k_final, zn_pair_tuples)
    
    integrator_prod = new_integrator()
    sim_prod = Simulation(topology, system_prod, integrator_prod, platform)
    sim_prod.context.setPositions(pos_settle)
    sim_prod.context.setVelocitiesToTemperature(temperature)
    
    traj_reporter = DCDReporter(str(traj_path), traj_interval_steps)
    
    prod_elapsed = run_stage(
        f"Production ({production_ns} ns)",
        sim_prod, production_steps, log_path, log_interval_steps,
        trajectory_reporter=traj_reporter,
        checkpoint_path=checkpoint_path,
        checkpoint_interval_steps=checkpoint_interval_steps,
        zn_ligand_pairs=zn_pairs,
        zn_log_path=zn_dist_log_path,
        zn_log_interval_steps=zn_log_steps,
    )
    
    # Final snapshot
    final_state = sim_prod.context.getState(getPositions=True)
    write_pdb(topology, final_state.getPositions(), final_pdb_path)
    
    # ========================================================================
    # Summary
    # ========================================================================
    total_elapsed = time.time() - start_time
    total_hours = total_elapsed / 3600
    
    print("\n" + "=" * 70)
    print("SIMULATION COMPLETE")
    print("=" * 70)
    print(f"Total runtime: {total_hours:.2f} hours ({total_elapsed:.1f} seconds)")
    print(f"Production time: {prod_elapsed:.1f} seconds")
    print(f"Output directory: {out_dir}")
    print(f"  - Log: md.log (single file with stage banners)")
    print(f"  - Trajectory: production.dcd")
    print(f"  - Zn distances: zn_distances.tsv (with stage column)")
    print(f"  - Final structure: final_frame.pdb")
    print("=" * 70)
    
    # Write summary file
    summary_path = out_dir / "simulation_summary.txt"
    with open(summary_path, "w") as f:
        f.write(f"ZFP36L2 + RNA 9-mer + 2 Zn MD Simulation Summary\n")
        f.write(f"Date: {datetime.now().isoformat()}\n")
        f.write(f"Input: {pdb_path}\n")
        f.write(f"\nParameters:\n")
        f.write(f"  Temperature: {temperature_k} K\n")
        f.write(f"  Pressure: {pressure_bar} bar\n")
        f.write(f"  Timestep: {timestep_fs} fs\n")
        f.write(f"  Padding: {padding_nm} nm\n")
        f.write(f"  Ionic strength: {ionic_strength} M\n")
        f.write(f"  Nonbonded cutoff: {nonbonded_cutoff_nm} nm\n")
        f.write(f"  NVT: {nvt_ps} ps\n")
        f.write(f"  NPT equilibration: {npt_eq_ps} ps\n")
        f.write(f"  NPT settle: {npt_settle_ns} ns\n")
        f.write(f"  Production: {production_ns} ns\n")
        f.write(f"\nZn Restraint Ramping:\n")
        f.write(f"  Minimization: k = {zn_k_weak} kJ/mol/nm² (WEAK)\n")
        f.write(f"  NVT + NPT-eq: k = {zn_k_medium} kJ/mol/nm² (MEDIUM)\n")
        f.write(f"  NPT-settle + Production: k = {zn_k_final} kJ/mol/nm² (FINAL)\n")
        f.write(f"\nRuntime:\n")
        f.write(f"  Total: {total_hours:.2f} hours\n")
        f.write(f"  Production: {prod_elapsed/3600:.2f} hours\n")
    
    return out_dir


# ============================================================================
# CLI
# ============================================================================

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="MD simulation for ZFP36L2 + RNA 9-mer + 2 Zn (with Zn restraint ramping)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    
    parser.add_argument(
        "--pdb", type=Path, default=DEFAULT_PDB,
        help="Input PDB file",
    )
    parser.add_argument(
        "--output-root", type=Path, default=OUTPUT_ROOT,
        help="Output directory root",
    )
    parser.add_argument(
        "--platform", type=str, default=None,
        help="OpenMM platform (CUDA, OpenCL, CPU, or None for auto)",
    )
    
    # Simulation parameters
    parser.add_argument("--padding-nm", type=float, default=1.0, help="Solvation padding (nm)")
    parser.add_argument("--ionic-strength", type=float, default=0.15, help="Ionic strength (M)")
    parser.add_argument("--temperature", type=float, default=300.0, help="Temperature (K)")
    parser.add_argument("--pressure", type=float, default=1.0, help="Pressure (bar)")
    parser.add_argument("--friction", type=float, default=1.0, help="Friction coefficient (1/ps)")
    parser.add_argument("--timestep", type=float, default=2.0, help="Timestep (fs)")
    parser.add_argument("--cutoff", type=float, default=1.2, help="Nonbonded cutoff (nm)")
    
    # Stage durations
    parser.add_argument("--minimize-iters", type=int, default=100000, help="Max minimization iterations")
    parser.add_argument("--nvt-ps", type=float, default=100.0, help="NVT equilibration (ps)")
    parser.add_argument("--npt-eq-ps", type=float, default=500.0, help="NPT equilibration (ps)")
    parser.add_argument("--npt-settle-ns", type=float, default=1.0, help="Restrained NPT settle (ns)")
    parser.add_argument("--production-ns", type=float, default=5.0, help="Production run (ns)")
    
    # Restraints
    parser.add_argument("--heavy-k", type=float, default=1000.0, help="Heavy-atom restraint k (kJ/mol/nm²)")
    parser.add_argument("--heavy-k-reduced", type=float, default=200.0, help="Reduced restraint k (kJ/mol/nm²)")
    parser.add_argument("--zn-k-weak", type=float, default=ZN_K_WEAK, help="Zn restraint k (WEAK, minimization)")
    parser.add_argument("--zn-k-medium", type=float, default=ZN_K_MEDIUM, help="Zn restraint k (MEDIUM, early equil)")
    parser.add_argument("--zn-k-final", type=float, default=ZN_K_FINAL, help="Zn restraint k (FINAL, settle+prod)")
    
    # Output intervals
    parser.add_argument("--log-ps", type=float, default=1.0, help="Log interval (ps)")
    parser.add_argument("--traj-ps", type=float, default=10.0, help="Trajectory save interval (ps)")
    parser.add_argument("--checkpoint-min", type=float, default=5.0, help="Checkpoint interval (minutes)")
    
    # Test mode
    parser.add_argument("--test", action="store_true", help="Run a quick test (short durations)")
    
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    
    # Test mode: short durations
    if args.test:
        print("[TEST MODE] Using short simulation times for testing")
        args.minimize_iters = 1000
        args.nvt_ps = 1.0
        args.npt_eq_ps = 1.0
        args.npt_settle_ns = 0.001
        args.production_ns = 0.002
    
    forcefield = create_forcefield()
    
    run_pipeline(
        pdb_path=args.pdb,
        output_root=args.output_root,
        forcefield=forcefield,
        platform_name=args.platform,
        padding_nm=args.padding_nm,
        ionic_strength=args.ionic_strength,
        temperature_k=args.temperature,
        pressure_bar=args.pressure,
        friction_per_ps=args.friction,
        timestep_fs=args.timestep,
        nonbonded_cutoff_nm=args.cutoff,
        minimize_max_iters=args.minimize_iters,
        nvt_ps=args.nvt_ps,
        npt_eq_ps=args.npt_eq_ps,
        npt_settle_ns=args.npt_settle_ns,
        production_ns=args.production_ns,
        heavy_restraint_k=args.heavy_k,
        heavy_restraint_k_reduced=args.heavy_k_reduced,
        zn_k_weak=args.zn_k_weak,
        zn_k_medium=args.zn_k_medium,
        zn_k_final=args.zn_k_final,
        log_interval_ps=args.log_ps,
        traj_interval_ps=args.traj_ps,
        checkpoint_interval_min=args.checkpoint_min,
    )


if __name__ == "__main__":
    main()
