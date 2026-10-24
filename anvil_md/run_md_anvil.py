#!/usr/bin/env python3
"""
ZFP36L2 + RNA 9-mer + 2 Zn — HPC GPU MD simulation script.

Designed for Anvil A100 / Expanse V100 with checkpoint/restart safety for
wall-time limits. Workflow matches the proven GCP run_minimal_md.py pipeline.

Force field: AMBER14 protein + RNA OL3 + TIP3P
Temperature: 310 K | Pressure: 1 atm | Timestep: 2 fs | Cutoff: 1.0 nm

Workflow:
  1. Solvate + ions (1.2 nm padding, 0.15 M)
  2. Energy minimization (10k iterations max)
  3. NVT equilibration (0.2 ns, 1000 kJ/mol/nm² heavy-atom restraints)
  4. NPT equilibration (0.5 ns, 100 kJ/mol/nm² restraints)
  5. Production NPT (target ns, Zn–Cys restraints only, checkpointed)

Usage:
  python run_md_anvil.py --pdb input.pdb --variant WT --production-ns 100
  python run_md_anvil.py --pdb input.pdb --variant WT --production-ns 100 --resume
"""

from __future__ import annotations

import argparse
import math
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from openmm import (
    CustomBondForce,
    CustomExternalForce,
    LangevinIntegrator,
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

RNA_BASES = {"A", "C", "G", "U"}
WALL_TIME_SAFETY_MARGIN_S = 30 * 60
ZN_RAMP_PS = 50.0  # ramp Zn restraint k over first 50 ps of production
ZN_K_RAMP_START = 50.0  # kJ/mol/nm² at production start


class _TimestampedThroughputReporter:
    """Append a wall-clock-timestamped throughput line every reportInterval steps.

    Writes a TSV: ISO-datetime, step, sim_ns, wall_seconds_since_start,
    interval_ns_per_day. Unlike OpenMM's "Speed" column (an instantaneous
    estimate that ignores I/O/checkpoint overhead), this records the REAL
    wall-clock pace between consecutive reports, so any future slow phase is
    visible after the fact with the exact time it occurred. This is the data we
    were missing when trying to diagnose why some runs were slow.
    """

    def __init__(self, file_path: str, report_interval: int, timestep_fs: float):
        self._interval = max(1, int(report_interval))
        self._timestep_fs = timestep_fs
        self._t0 = None
        self._t_prev = None
        self._step_prev = None
        new = not Path(file_path).exists()
        self._fh = open(file_path, "a")
        if new:
            self._fh.write("datetime\tstep\tsim_ns\twall_s\tinterval_ns_per_day\n")
            self._fh.flush()

    def describeNextReport(self, simulation):
        steps = self._interval - simulation.currentStep % self._interval
        # (steps, needPositions, needVelocities, needForces, needEnergy)
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        now = time.time()
        step = simulation.currentStep
        sim_ns = step * self._timestep_fs / 1e6
        if self._t0 is None:
            self._t0 = now
        rate = 0.0
        if self._t_prev is not None and now > self._t_prev:
            d_ns = (step - self._step_prev) * self._timestep_fs / 1e6
            rate = d_ns / ((now - self._t_prev) / 86400.0)
        self._fh.write(
            f"{datetime.now().isoformat()}\t{step}\t{sim_ns:.4f}\t"
            f"{now - self._t0:.1f}\t{rate:.1f}\n"
        )
        self._fh.flush()
        self._t_prev = now
        self._step_prev = step

    def __del__(self):
        try:
            self._fh.close()
        except Exception:
            pass



def time_to_steps(duration_ns: float, timestep_fs: float) -> int:
    total_ps = duration_ns * 1000.0
    step_ps = timestep_fs / 1000.0
    return int(total_ps / step_ps)


def ps_to_steps(interval_ps: float, timestep_fs: float) -> int:
    step_ps = timestep_fs / 1000.0
    return max(1, int(interval_ps / step_ps))


def steps_to_ns(steps: int, timestep_fs: float) -> float:
    return steps * timestep_fs / 1e6


def create_forcefield() -> ForceField:
    combos = [
        ("amber14/protein.ff14SB.xml", "amber14/RNA.OL3.xml", "amber14/tip3p.xml"),
        ("amber14-all.xml", "amber14/tip3p.xml"),
        ("amber14-all.xml", "tip3p.xml"),
    ]
    last_error = None
    for files in combos:
        try:
            ff = ForceField(*files)
            print(f"[ForceField] Loaded: {files}")
            return ff
        except ValueError as exc:
            last_error = exc
    raise last_error


def _canonical_rna_base(name: str) -> Optional[str]:
    upper = name.strip().upper()
    if not upper:
        return None
    letter = upper[-1]
    return letter if letter in RNA_BASES else None


def identify_terminal_template_ids(topology) -> Dict[Tuple[str, str], str]:
    template_ids = {}
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


def build_residue_template_map(topology, template_ids):
    mapping = {}
    for residue in topology.residues():
        key = (residue.chain.id, residue.id)
        template = template_ids.get(key)
        if template:
            mapping[residue] = template
    return mapping


def add_position_restraints(system, topology, reference_positions, force_constant, heavy_only=True):
    # PBC-safe restraints required for CUDA under periodic boundary conditions
    force = CustomExternalForce("0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
    force.addGlobalParameter("k", force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    positions_nm = reference_positions.value_in_unit(unit.nanometer)
    n = 0
    for atom in topology.atoms():
        if heavy_only and atom.element is not None and atom.element.symbol == "H":
            continue
        pos = positions_nm[atom.index]
        force.addParticle(atom.index, pos)
        n += 1
    system.addForce(force)
    print(f"[Restraints] {n} position restraints added (k={force_constant})")
    return n


def add_zinc_restraints(system, topology, positions, force_constant, cutoff_nm=0.35):
    atoms = list(topology.atoms())
    zn_indices = [a.index for a in atoms if a.element == Element.getBySymbol("Zn")]
    cys_sg = {a.index: a for a in atoms if a.residue.name == "CYS" and a.name == "SG"}
    if not zn_indices or not cys_sg:
        print("[Zn] No Zn or Cys SG found; skipping Zn restraints")
        return 0

    positions_nm = positions.value_in_unit(unit.nanometer)
    zn_force = CustomBondForce("0.5*k*(r-r0)^2")
    zn_force.addGlobalParameter("k", force_constant)
    zn_force.addPerBondParameter("r0")

    bonds_added = 0
    for zn_idx in zn_indices:
        zn_pos = positions_nm[zn_idx]
        for sg_idx in cys_sg:
            sg_pos = positions_nm[sg_idx]
            dist = math.sqrt(
                (zn_pos[0] - sg_pos[0]) ** 2
                + (zn_pos[1] - sg_pos[1]) ** 2
                + (zn_pos[2] - sg_pos[2]) ** 2
            )
            if dist <= cutoff_nm:
                zn_force.addBond(zn_idx, sg_idx, [dist])
                bonds_added += 1
    if bonds_added:
        system.addForce(zn_force)
        print(f"[Zn] {bonds_added} Zn–Cys restraints added (k={force_constant})")
    return bonds_added


def build_modeller(pdb_path, forcefield, padding_nm, ionic_strength):
    print(f"[Setup] Loading PDB: {pdb_path}")
    pdb = PDBFile(str(pdb_path))
    modeller = Modeller(pdb.topology, pdb.positions)

    template_ids = identify_terminal_template_ids(modeller.topology)
    res_templates = build_residue_template_map(modeller.topology, template_ids)

    try:
        modeller.addHydrogens(forcefield, pH=7.0, residueTemplates=res_templates)
    except Exception as exc:
        if "User-supplied template" not in str(exc):
            raise
        modeller.addHydrogens(forcefield, pH=7.0)
    res_templates = build_residue_template_map(modeller.topology, template_ids)

    try:
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=padding_nm * unit.nanometer,
            ionicStrength=ionic_strength * unit.molar,
            residueTemplates=res_templates,
        )
    except Exception as exc:
        if "User-supplied template" not in str(exc):
            raise
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=padding_nm * unit.nanometer,
            ionicStrength=ionic_strength * unit.molar,
        )

    n_atoms = modeller.topology.getNumAtoms()
    print(f"[Setup] Solvated system: {n_atoms} atoms")
    return modeller, template_ids


def create_system(forcefield, topology, res_templates, with_barostat, temperature, pressure):
    kwargs = dict(
        topology=topology,
        nonbondedMethod=PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=HBonds,
        rigidWater=True,
        ewaldErrorTolerance=0.0005,
        removeCMMotion=True,
    )
    try:
        system = forcefield.createSystem(
            residueTemplates=res_templates or {},
            **kwargs,
        )
    except Exception as exc:
        if not res_templates or "User-supplied template" not in str(exc):
            raise
        system = forcefield.createSystem(**kwargs)
    if with_barostat:
        system.addForce(MonteCarloBarostat(pressure, temperature, 25))
    return system


def get_platform_and_properties(platform_name: Optional[str], precision: str):
    if platform_name:
        platform = Platform.getPlatformByName(platform_name)
    else:
        platform = None
        for name in ("CUDA", "OpenCL", "CPU"):
            try:
                platform = Platform.getPlatformByName(name)
                break
            except Exception:
                continue

    properties = None
    if platform and platform.getName() in ("CUDA", "OpenCL"):
        properties = {"Precision": precision}
        if platform.getName() == "CUDA":
            properties["DeviceIndex"] = "0"
        print(f"[Platform] {platform.getName()} (precision={precision})")
    elif platform:
        print(f"[Platform] {platform.getName()}")
    else:
        print("[Platform] auto (OpenMM default)")
    return platform, properties


def _is_nan_coordinate_error(exc: Exception) -> bool:
    message = str(exc).lower()
    return "particle coordinate is nan" in message or "particle position is nan" in message


def is_valid_dcd(path: Path) -> bool:
    """Return True if path looks like a non-empty DCD with a valid header."""
    if not path.exists():
        return False
    if path.stat().st_size < 100:
        return False
    with path.open("rb") as handle:
        return len(handle.read(4)) == 4


def prepare_dcd_path(path: Path, append: bool) -> bool:
    """
    Prepare trajectory file for DCDReporter.

    Returns the append flag to pass to DCDReporter.

    IMPORTANT (fix for the 2026-05/06 truncation bug):
    If a valid DCD already exists we ALWAYS append to it rather than deleting
    it — even when the caller did not explicitly request append.  This means a
    resubmitted job that forgot the ``--resume`` flag will *extend* the
    existing trajectory instead of silently overwriting the first part of it.
    Only a genuinely fresh start (no existing DCD, or an explicit
    ``fresh_start`` NaN-retry) creates a new file.
    """
    if is_valid_dcd(path):
        # A real trajectory is on disk: preserve it and append, regardless of
        # the requested append flag. Overwriting would lose completed frames.
        if not append:
            print(
                f"[DCD] Existing trajectory found at {path.name}; appending to "
                f"it instead of overwriting (trajectory-preservation safeguard)."
            )
        return True
    if path.exists():
        # Exists but is not a valid DCD (e.g. truncated/zero-byte) -> reset it.
        path.unlink()
    return False


def clear_production_restart_artifacts(traj_path: Path, prod_checkpoint: Path):
    """Remove partial production outputs before a clean production restart."""
    if traj_path.exists() and not is_valid_dcd(traj_path):
        traj_path.unlink()
    if prod_checkpoint.exists():
        prod_checkpoint.unlink()


def load_equilibrated_state(equil_pdb: Path, equil_checkpoint: Path):
    """Load equilibrated coordinates (and velocities if checkpoint exists)."""
    pdb = PDBFile(str(equil_pdb))
    velocities = None
    if equil_checkpoint.exists():
        print(f"[Skip-equil] Found checkpoint: {equil_checkpoint}")
    return pdb.topology, pdb.positions, velocities


def set_zn_restraint_k(context, k_value):
    """Update Zn–Cys restraint strength if the force is present."""
    try:
        context.setParameter("k", k_value)
    except Exception:
        pass


def run_pipeline(args):
    variant = args.variant
    pdb_path = Path(args.pdb).resolve()
    out_dir = Path(args.output_dir) / variant
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / "md.log"
    timing_log_path = out_dir / "throughput.tsv"
    traj_path = out_dir / "production.dcd"
    prod_checkpoint = out_dir / "production.chk"

    equil_checkpoint = out_dir / "equilibrated.chk"
    equil_pdb = out_dir / "equilibrated.pdb"
    final_pdb = out_dir / "final_frame.pdb"
    progress_file = out_dir / "progress.txt"

    timestep_fs = args.timestep_fs
    temperature = args.temperature * unit.kelvin
    pressure = args.pressure * unit.atmospheres
    friction = 1.0 / unit.picoseconds
    timestep = timestep_fs * unit.femtoseconds
    production_ns = args.production_ns
    production_steps = time_to_steps(production_ns, timestep_fs)

    log_interval = ps_to_steps(args.log_interval_ps, timestep_fs)
    traj_interval = ps_to_steps(args.traj_interval_ps, timestep_fs)
    checkpoint_interval = ps_to_steps(15 * 60 * 1000, timestep_fs)

    nvt_steps = time_to_steps(args.stage1_ns, timestep_fs)
    npt_steps = time_to_steps(args.stage2_ns, timestep_fs)

    wall_limit_s = args.wall_hours * 3600
    job_start_time = time.time()

    def time_remaining():
        return wall_limit_s - (time.time() - job_start_time) - WALL_TIME_SAFETY_MARGIN_S

    print("=" * 70)
    print(f"ZFP36L2 variant={variant}  production={production_ns} ns  T={args.temperature}K")
    print(f"Padding={args.padding} nm  Ionic={args.ionic_strength} M  P={args.pressure} atm")
    print(f"Wall limit={args.wall_hours}h  Safety margin={WALL_TIME_SAFETY_MARGIN_S/60:.0f}min")
    print(f"Output: {out_dir}")
    print("=" * 70)

    forcefield = create_forcefield()
    platform, properties = get_platform_and_properties(args.platform, args.precision)

    resuming = args.resume and prod_checkpoint.exists() and equil_checkpoint.exists()
    skip_equil = args.skip_equil and equil_pdb.exists() and not resuming
    if resuming:
        print("\n[RESUME] Found production checkpoint. Resuming production run...")
    elif skip_equil:
        print(f"\n[Skip-equil] Using existing equilibrated structure: {equil_pdb}")

    if skip_equil or resuming:
        if not equil_pdb.exists():
            raise FileNotFoundError(
                f"Equilibrated PDB required for resume/skip-equil: {equil_pdb}"
            )
        topology, equil_positions, _ = load_equilibrated_state(equil_pdb, equil_checkpoint)
        positions = equil_positions
        template_ids = identify_terminal_template_ids(topology)
        res_templates = build_residue_template_map(topology, template_ids)
        equil_velocities = None
        if resuming:
            print(f"[RESUME] Loaded solvated topology from {equil_pdb} (skipping re-solvation)")
    else:
        modeller, template_ids = build_modeller(pdb_path, forcefield, args.padding, args.ionic_strength)
        topology = modeller.topology
        positions = modeller.positions
        res_templates = build_residue_template_map(topology, template_ids)
        equil_positions = None
        equil_velocities = None

    def new_integrator(local_timestep_fs: float = timestep_fs):
        return LangevinIntegrator(
            temperature, friction, local_timestep_fs * unit.femtoseconds
        )

    def make_simulation(system, local_timestep_fs: float = timestep_fs):
        return Simulation(
            topology, system, new_integrator(local_timestep_fs), platform, properties
        )

    if not resuming and not skip_equil:
        # Stage 1: Minimization + NVT (same context, no Zn restraints)
        print("\n" + "=" * 70)
        print(f"Stage 1: Minimization + NVT ({args.stage1_ns} ns)")
        print("=" * 70)
        sys_stage1 = create_system(forcefield, topology, res_templates, False, temperature, pressure)
        add_position_restraints(
            sys_stage1, topology, positions,
            1000.0 * unit.kilojoule_per_mole / unit.nanometer**2,
        )
        sim_stage1 = make_simulation(sys_stage1)
        sim_stage1.context.setPositions(positions)
        print(f"[Minimize] Running ({args.minimize_max_steps} max iterations)...")
        sim_stage1.minimizeEnergy(maxIterations=args.minimize_max_steps)
        sim_stage1.context.setVelocitiesToTemperature(temperature)
        sim_stage1.reporters.append(StateDataReporter(
            str(log_path), log_interval, step=True, potentialEnergy=True,
            kineticEnergy=True, temperature=True, speed=True, totalSteps=nvt_steps,
            separator="\t", append=False,
        ))
        print(f"[NVT] Running {nvt_steps} steps...")
        sim_stage1.step(nvt_steps)
        state1 = sim_stage1.context.getState(getPositions=True, getVelocities=True)
        pos_nvt = state1.getPositions()
        sim_stage1.reporters.clear()
        del sim_stage1

        # Stage 2: NPT equilibration (no Zn restraints)
        print("\n" + "=" * 70)
        print(f"Stage 2: NPT Equilibration ({args.stage2_ns} ns)")
        print("=" * 70)
        sys_npt = create_system(forcefield, topology, res_templates, True, temperature, pressure)
        add_position_restraints(
            sys_npt, topology, pos_nvt,
            100.0 * unit.kilojoule_per_mole / unit.nanometer**2,
        )
        sim_npt = make_simulation(sys_npt)
        sim_npt.context.setPositions(pos_nvt)
        sim_npt.context.setVelocitiesToTemperature(temperature)
        sim_npt.reporters.append(StateDataReporter(
            str(log_path), log_interval, step=True, potentialEnergy=True,
            temperature=True, volume=True, density=True, speed=True,
            totalSteps=npt_steps, separator="\t", append=True,
        ))
        print(f"[NPT-eq] Running {npt_steps} steps...")
        sim_npt.step(npt_steps)
        state_npt = sim_npt.context.getState(getPositions=True, getVelocities=True)
        pos_npt = state_npt.getPositions()
        equil_velocities = state_npt.getVelocities()
        sim_npt.reporters.clear()
        del sim_npt

        with open(equil_pdb, "w") as f:
            PDBFile.writeFile(topology, pos_npt, f, keepIds=True)
        print(f"[Output] Equilibrated PDB: {equil_pdb}")

        equil_positions = pos_npt

    # Stage 3: Production NPT
    print("\n" + "=" * 70)
    print(f"Stage 3: Production NPT ({production_ns} ns)")
    print("=" * 70)

    zinc_k = args.zinc_k * unit.kilojoule_per_mole / unit.nanometer**2
    zn_ramp_steps = ps_to_steps(ZN_RAMP_PS, timestep_fs)
    zn_k_ramp_start = ZN_K_RAMP_START * unit.kilojoule_per_mole / unit.nanometer**2

    def setup_production_sim(
        local_timestep_fs: float,
        extra_heavy_k: Optional[float],
        append_traj: bool,
        local_zinc_k=None,
        fresh_start: bool = False,
    ):
        if fresh_start:
            clear_production_restart_artifacts(traj_path, prod_checkpoint)

        sys_prod = create_system(forcefield, topology, res_templates, True, temperature, pressure)
        ref_positions = equil_positions if equil_positions is not None else positions
        effective_zinc_k = local_zinc_k if local_zinc_k is not None else zinc_k
        zn_added = 0
        if not args.no_zinc_restraints:
            zn_added = add_zinc_restraints(sys_prod, topology, ref_positions, effective_zinc_k)
        if extra_heavy_k is not None:
            add_position_restraints(
                sys_prod, topology, ref_positions,
                extra_heavy_k * unit.kilojoule_per_mole / unit.nanometer**2,
            )
        sim = make_simulation(sys_prod, local_timestep_fs)
        if resuming and not fresh_start:
            sim.loadCheckpoint(str(prod_checkpoint))
        else:
            sim.context.setPositions(equil_positions)
            if equil_velocities is not None and not fresh_start:
                sim.context.setVelocities(equil_velocities)
            else:
                sim.context.setVelocitiesToTemperature(temperature)
            if not skip_equil or not equil_checkpoint.exists():
                sim.saveCheckpoint(str(equil_checkpoint))
                print(f"[Checkpoint] Equilibration state saved: {equil_checkpoint}")
            elif skip_equil:
                print(f"[Skip-equil] Reusing equilibration checkpoint: {equil_checkpoint}")
        local_steps = time_to_steps(production_ns, local_timestep_fs)
        local_log_interval = ps_to_steps(args.log_interval_ps, local_timestep_fs)
        local_traj_interval = ps_to_steps(args.traj_interval_ps, local_timestep_fs)
        dcd_append = prepare_dcd_path(traj_path, append_traj)
        sim.reporters.append(DCDReporter(str(traj_path), local_traj_interval, append=dcd_append))
        sim.reporters.append(StateDataReporter(
            str(log_path), local_log_interval, step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, volume=True, density=True,
            speed=True, remainingTime=True, progress=True,
            elapsedTime=True,
            totalSteps=local_steps, separator="\t", append=True,
        ))
        # Wall-clock timestamped throughput log so any future slow phase is
        # diagnosable after the fact (records real datetime + step + ns each
        # reporting interval). This is independent of OpenMM's "Speed" column,
        # which is only an instantaneous estimate.
        sim.reporters.append(
            _TimestampedThroughputReporter(str(timing_log_path), local_log_interval, timestep_fs)
        )

        sim.reporters.append(CheckpointReporter(str(prod_checkpoint), checkpoint_interval))
        if zn_added:
            set_zn_restraint_k(sim.context, zn_k_ramp_start)
            print(f"[Zn] Ramping restraint k: {ZN_K_RAMP_START} → {args.zinc_k} kJ/mol/nm² "
                  f"over {ZN_RAMP_PS} ps")
        return sim, local_steps, zn_added

    def run_production_chunks(sim, steps_remaining, steps_done_local, ramp_zn: bool):
        nonlocal prod_start
        chunk_steps = ps_to_steps(500.0, timestep_fs)
        total_remaining = steps_remaining
        steps_at_start = steps_done_local
        production_steps_local = steps_at_start + total_remaining
        ramp_done = not ramp_zn

        while total_remaining > 0:
            tr = time_remaining()
            if tr <= 0:
                print("\n[WALL TIME] Approaching wall limit. Saving checkpoint and exiting.")
                sim.saveCheckpoint(str(prod_checkpoint))
                steps_completed = production_steps_local - total_remaining
                ns_completed = steps_to_ns(steps_completed, timestep_fs)
                with open(progress_file, "w") as f:
                    f.write(f"variant={variant}\n")
                    f.write(f"steps_completed={steps_completed}\n")
                    f.write(f"ns_completed={ns_completed:.4f}\n")
                    f.write(f"target_ns={production_ns}\n")
                    f.write(f"timestamp={datetime.now().isoformat()}\n")
                    f.write(f"status=checkpointed\n")
                print(f"[WALL TIME] Saved at {ns_completed:.2f} / {production_ns} ns")
                final_state = sim.context.getState(getPositions=True)
                with open(final_pdb, "w") as f:
                    PDBFile.writeFile(topology, final_state.getPositions(), f, keepIds=True)
                return steps_completed, True

            if ramp_zn and not ramp_done:
                steps_completed_so_far = production_steps_local - total_remaining
                if steps_completed_so_far >= zn_ramp_steps:
                    set_zn_restraint_k(sim.context, zinc_k)
                    ramp_done = True
                    print(f"[Zn] Restraint k reached target: {args.zinc_k} kJ/mol/nm²")
                else:
                    frac = steps_completed_so_far / zn_ramp_steps
                    k_now = ZN_K_RAMP_START + frac * (args.zinc_k - ZN_K_RAMP_START)
                    set_zn_restraint_k(
                        sim.context,
                        k_now * unit.kilojoule_per_mole / unit.nanometer**2,
                    )

            this_chunk = min(chunk_steps, total_remaining)
            sim.step(this_chunk)
            total_remaining -= this_chunk

            steps_completed = production_steps_local - total_remaining
            ns_completed = steps_to_ns(steps_completed, timestep_fs)
            elapsed = time.time() - prod_start
            if elapsed > 0:
                ns_per_day = (ns_completed - steps_to_ns(steps_at_start, timestep_fs)) / (elapsed / 86400)
            else:
                ns_per_day = 0

            if total_remaining % (chunk_steps * 10) == 0 or total_remaining == 0:
                print(f"[Production] {ns_completed:.2f}/{production_ns} ns "
                      f"({100*steps_completed/production_steps_local:.1f}%) "
                      f"Speed: {ns_per_day:.1f} ns/day  "
                      f"Wall remaining: {time_remaining()/3600:.2f}h")

        return production_steps_local, False

    if resuming:
        sim_prod, production_steps, zn_added = setup_production_sim(
            timestep_fs, None, append_traj=True, fresh_start=False,
        )
        steps_done = sim_prod.currentStep
        steps_remaining = production_steps - steps_done
        print(f"[RESUME] Completed {steps_to_ns(steps_done, timestep_fs):.2f} ns ({steps_done} steps)")
        print(f"[RESUME] Remaining: {steps_to_ns(steps_remaining, timestep_fs):.2f} ns")
        if steps_remaining <= 0:
            print("[RESUME] Production already complete!")
            final_state = sim_prod.context.getState(getPositions=True)
            with open(final_pdb, "w") as f:
                PDBFile.writeFile(topology, final_state.getPositions(), f, keepIds=True)
            return
        resume_extend = True
    else:
        # Fresh (non-resume) start. Only clear an *invalid/partial* DCD — a
        # valid existing trajectory is preserved & appended to by
        # prepare_dcd_path() below. This prevents the truncation bug where a
        # resubmitted job without --resume overwrote completed frames.
        if not is_valid_dcd(traj_path):
            clear_production_restart_artifacts(traj_path, prod_checkpoint)
        sim_prod, production_steps, zn_added = setup_production_sim(
            timestep_fs, None, append_traj=False, fresh_start=False,
        )
        steps_done = 0
        steps_remaining = production_steps
        resume_extend = False
        if is_valid_dcd(traj_path):
            print(
                "[WARN] A valid production.dcd already exists but --resume was "
                "NOT given. Appending new frames to preserve existing data. "
                "If you intended a clean restart, delete production.dcd first."
            )

    print(f"[Production] Starting: {steps_to_ns(steps_remaining, timestep_fs):.2f} ns remaining")
    print(f"[Production] Wall time remaining: {time_remaining()/3600:.2f} hours")

    prod_start = time.time()
    production_failed = False
    nan_retry_mode = None

    try:
        steps_completed, wall_exit = run_production_chunks(
            sim_prod, steps_remaining, steps_done,
            ramp_zn=bool(zn_added) and not resume_extend,
        )
        if wall_exit:
            return
        production_steps = steps_completed
    except Exception as exc:
        if not _is_nan_coordinate_error(exc):
            raise
        production_failed = True
        nan_retry_mode = "1fs_25kJ_zn100"
        print(
            "[retry] NaN during production. Restarting with 1.0 fs timestep, "
            "25 kJ heavy-atom restraints, and Zn k=100."
        )
        sim_prod.reporters.clear()
        del sim_prod

        retry_zinc_k = 100.0 * unit.kilojoule_per_mole / unit.nanometer**2
        sim_prod, production_steps, zn_added = setup_production_sim(
            1.0,
            25.0,
            append_traj=False,
            local_zinc_k=retry_zinc_k,
            fresh_start=True,
        )
        prod_start = time.time()
        steps_completed, wall_exit = run_production_chunks(
            sim_prod, production_steps, 0, ramp_zn=False,
        )
        if wall_exit:
            return

    prod_elapsed = time.time() - prod_start
    sim_prod.saveCheckpoint(str(prod_checkpoint))

    final_state = sim_prod.context.getState(getPositions=True)
    with open(final_pdb, "w") as f:
        PDBFile.writeFile(topology, final_state.getPositions(), f, keepIds=True)

    with open(progress_file, "w") as f:
        f.write(f"variant={variant}\n")
        f.write(f"steps_completed={production_steps}\n")
        f.write(f"ns_completed={production_ns:.4f}\n")
        f.write(f"target_ns={production_ns}\n")
        f.write(f"timestamp={datetime.now().isoformat()}\n")
        f.write(f"status=completed\n")
        f.write(f"production_hours={prod_elapsed/3600:.2f}\n")
        if production_failed:
            f.write(f"nan_retry={nan_retry_mode}\n")

    print("\n" + "=" * 70)
    print("SIMULATION COMPLETE")
    print(f"Production: {production_ns} ns in {prod_elapsed/3600:.2f} hours")
    print(f"Output: {out_dir}")
    print("=" * 70)


def parse_args():
    parser = argparse.ArgumentParser(
        description="ZFP36L2 MD simulation for HPC GPU (Anvil/Expanse)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--pdb", type=str, required=True, help="Input cleaned PDB file")
    parser.add_argument("--variant", type=str, required=True, help="Variant name (e.g., WT, Y154H)")
    parser.add_argument("--output-dir", type=str, default="./md_output", help="Output directory")
    parser.add_argument("--production-ns", type=float, default=100.0, help="Production run length (ns)")
    parser.add_argument("--padding", type=float, default=1.2, help="Solvation padding (nm)")
    parser.add_argument("--ionic-strength", type=float, default=0.15, help="Ionic strength (M)")
    parser.add_argument("--temperature", type=float, default=310.0, help="Temperature (K)")
    parser.add_argument("--pressure", type=float, default=1.0, help="Pressure (atm)")
    parser.add_argument("--timestep-fs", type=float, default=2.0, help="Integrator timestep (fs)")
    parser.add_argument("--stage1-ns", type=float, default=0.2, help="NVT equilibration duration (ns)")
    parser.add_argument("--stage2-ns", type=float, default=0.5, help="NPT equilibration duration (ns)")
    parser.add_argument("--minimize-max-steps", type=int, default=10000, help="Max minimization iterations")
    parser.add_argument("--zinc-k", type=float, default=200.0, help="Zn–Cys restraint constant (kJ/mol/nm²)")
    parser.add_argument("--no-zinc-restraints", action="store_true", help="Disable Zn restraints in production")
    parser.add_argument("--log-interval-ps", type=float, default=10.0, help="Log interval (ps)")
    parser.add_argument("--traj-interval-ps", type=float, default=100.0, help="Trajectory interval (ps)")
    parser.add_argument("--platform", type=str, default="CUDA", help="OpenMM platform (CUDA, CPU, OpenCL)")
    parser.add_argument("--precision", type=str, default="single", choices=["single", "double", "mixed"],
                        help="CUDA/OpenCL precision (never use mixed for this system)")
    parser.add_argument("--wall-hours", type=float, default=47.5, help="Wall-time limit (hours)")
    parser.add_argument("--resume", action="store_true", help="Resume from checkpoint if available")
    parser.add_argument(
        "--skip-equil",
        action="store_true",
        help="Skip solvation/equilibration; start production from equilibrated.pdb",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_pipeline(args)
