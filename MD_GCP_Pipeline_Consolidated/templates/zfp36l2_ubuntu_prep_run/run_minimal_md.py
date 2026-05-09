#!/usr/bin/env python3
"""
Minimal MD driver for ZFP36L2–GM-CSF RNA complexes.

Workflow per complex:
- Solvate the cleaned protein/RNA PDB (chains A/B, Zn included)
- Energy minimization (steepest descent, 10k iterations max)
- NVT equilibration (0.2 ns, 1000 kJ/mol/nm^2 heavy-atom restraints)
- NPT equilibration (0.5 ns, 100 kJ/mol/nm^2 restraints, 1 atm)
- NPT production (10 ns, Zn–Cys restraints only, 2 ps trajectory stride)

Outputs are written under prep_run/minimal_md/<pdb_stem>/:
- md.log (aggregated reporter spanning all stages)
- production.dcd (2 ps stride)
- final_frame.pdb (last snapshot)
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
from openmm import CustomBondForce, CustomExternalForce, LangevinIntegrator, MonteCarloBarostat, Platform, unit
from openmm.app import (
    DCDReporter,
    Element,
    ForceField,
    Modeller,
    PDBFile,
    PME,
    Simulation,
    StateDataReporter,
    HBonds,
)


BASE_DIR = Path(__file__).resolve().parent
CLEAN_DIR = BASE_DIR / "cleaned_pdbs"
DEFAULT_PDBS = [
    CLEAN_DIR / "fold_zfp36l2_wt_gm_csf_wt_model_0_cleaned.pdb",
    CLEAN_DIR / "fold_zfp36l2_c206f_gm_csf_wt_model_0_cleaned.pdb",
    CLEAN_DIR / "fold_zfp36l2_d210e_gm_csf_wt_model_0_cleaned.pdb",
    CLEAN_DIR / "fold_zfp36l2_wt_gm_csf_gc_mutant_model_0_cleaned.pdb",
]

RNA_BASES = {"A", "C", "G", "U"}


def time_to_steps(duration_ns: float, timestep_fs: float) -> int:
    total_ps = duration_ns * 1000.0
    step_ps = timestep_fs / 1000.0
    return int(total_ps / step_ps)


def ps_to_steps(interval_ps: float, timestep_fs: float) -> int:
    step_ps = timestep_fs / 1000.0
    return max(1, int(interval_ps / step_ps))


def create_forcefield() -> ForceField:
    """
    Load Amber14 protein/RNA + TIP3P water forcefields.

    We prefer the standalone RNA OL3 file (older releases used different file
    names), but if it is unavailable fall back to the aggregated Amber14 file.
    """

    combos = [
        (
            "amber14/protein.ff14SB.xml",
            "amber14/RNA.OL3.xml",
            "amber14/tip3p.xml",
        ),
        ("amber14-all.xml", "amber14/tip3p.xml"),
        ("amber14-all.xml", "tip3p.xml"),
    ]

    last_error: ValueError | None = None
    for files in combos:
        try:
            return ForceField(*files)
        except ValueError as exc:
            last_error = exc
    raise last_error


def _canonical_rna_base(name: str) -> str | None:
    upper = name.strip().upper()
    if not upper:
        return None
    letter = upper[-1]
    return letter if letter in RNA_BASES else None


def identify_terminal_template_ids(topology) -> dict[Tuple[str, str], str]:
    """
    Determine which residues should use 5'/3' RNA templates based on chain content.

    Returns
    -------
    dict mapping (chain_id, residue_id) -> template name (e.g., "A5", "U3")
    """

    template_ids: dict[Tuple[str, str], str] = {}
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


def build_residue_template_map(topology, template_ids: dict[Tuple[str, str], str]):
    mapping = {}
    for residue in topology.residues():
        key = (residue.chain.id, residue.id)
        template = template_ids.get(key)
        if template:
            mapping[residue] = template
    return mapping


def add_position_restraints(system, topology, reference_positions, force_constant, heavy_only=True):
    force = CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    force.addGlobalParameter("k", force_constant)
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    positions_nm = reference_positions.value_in_unit(unit.nanometer)
    for atom in topology.atoms():
        if heavy_only and atom.element is not None and atom.element.symbol == "H":
            continue
        pos = positions_nm[atom.index]
        force.addParticle(atom.index, pos)
    system.addForce(force)


def add_zinc_restraints(system, topology, positions, force_constant, cutoff_nm=0.35):
    atoms = list(topology.atoms())
    zn_indices = [a.index for a in atoms if a.element == Element.getBySymbol("Zn")]
    cys_sg = {
        a.index: a for a in atoms if a.residue.name == "CYS" and a.name == "SG"
    }
    if not zn_indices or not cys_sg:
        return

    positions_nm = positions.value_in_unit(unit.nanometer)
    zn_force = CustomBondForce("0.5*k*(r-r0)^2")
    zn_force.addGlobalParameter("k", force_constant)
    zn_force.addPerBondParameter("r0")

    cutoff = cutoff_nm
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
            if dist <= cutoff:
                zn_force.addBond(zn_idx, sg_idx, [dist])
                bonds_added += 1
    if bonds_added:
        system.addForce(zn_force)


def build_modeller(pdb_path: Path, forcefield: ForceField, padding_nm: float, ionic_strength: float):
    pdb = PDBFile(str(pdb_path))
    modeller = Modeller(pdb.topology, pdb.positions)
    template_ids = identify_terminal_template_ids(modeller.topology)
    residue_templates = build_residue_template_map(modeller.topology, template_ids)
    try:
        modeller.addHydrogens(forcefield, pH=7.0, residueTemplates=residue_templates)
    except Exception as exc:
        # Some AF3-derived RNA residues don't match forced terminal templates.
        # Fall back to automatic template matching for broader compatibility.
        if "User-supplied template" not in str(exc):
            raise
        modeller.addHydrogens(forcefield, pH=7.0)
    residue_templates = build_residue_template_map(modeller.topology, template_ids)
    try:
        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=padding_nm * unit.nanometer,
            ionicStrength=ionic_strength * unit.molar,
            residueTemplates=residue_templates,
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
    return modeller, template_ids


def create_system(
    forcefield: ForceField,
    topology,
    residue_templates: dict | None,
    with_barostat: bool,
    temperature,
    pressure,
):
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
            residueTemplates=residue_templates or {},
            **kwargs,
        )
    except Exception as exc:
        if not residue_templates or "User-supplied template" not in str(exc):
            raise
        system = forcefield.createSystem(**kwargs)
    if with_barostat:
        system.addForce(MonteCarloBarostat(pressure, temperature, 25))
    return system


def get_platform(name: str | None):
    if name is None:
        return None
    return Platform.getPlatformByName(name)


def run_stage(
    stage_name: str,
    simulation: Simulation,
    steps: int,
    log_path: Path,
    log_interval_steps: int,
    append_log: bool,
    trajectory_reporter: DCDReporter | None = None,
):
    reporters: List = [
        StateDataReporter(
            str(log_path),
            log_interval_steps,
            step=True,
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
            append=append_log,
        )
    ]
    if trajectory_reporter is not None:
        reporters.append(trajectory_reporter)
    simulation.reporters.extend(reporters)
    print(f"[{stage_name}] Running {steps} steps...")
    simulation.step(steps)
    simulation.reporters.clear()


def write_final_snapshot(topology, positions, path: Path):
    with path.open("w") as handle:
        PDBFile.writeFile(topology, positions, handle, keepIds=True)


def _is_nan_coordinate_error(exc: Exception) -> bool:
    message = str(exc).lower()
    return "particle coordinate is nan" in message or "particle position is nan" in message


def run_pipeline(
    pdb_path: Path,
    output_root: Path,
    forcefield: ForceField,
    timestep_fs: float,
    platform_name: str | None,
    padding_nm: float,
    ionic_strength: float,
    temperature: float,
    pressure: float,
    restrain_zinc: bool,
    zinc_k: float,
    stage1_ns: float,
    stage2_ns: float,
    production_ns: float,
    log_interval_ps: float,
    traj_interval_ps: float,
    minimize_max_steps: int,
):
    pdb_path = pdb_path.resolve()
    if not pdb_path.exists():
        raise FileNotFoundError(pdb_path)
    system_name = pdb_path.stem
    out_dir = output_root / system_name
    out_dir.mkdir(parents=True, exist_ok=True)

    log_path = out_dir / "md.log"
    traj_path = out_dir / "production.dcd"
    final_pdb_path = out_dir / "final_frame.pdb"

    print(f"==> Preparing {system_name}")
    modeller, template_ids = build_modeller(pdb_path, forcefield, padding_nm, ionic_strength)
    topology = modeller.topology
    positions = modeller.positions
    residue_templates = build_residue_template_map(topology, template_ids)

    platform = get_platform(platform_name)
    timestep = timestep_fs * unit.femtoseconds
    temperature_k = temperature * unit.kelvin
    pressure_atm = pressure * unit.atmospheres

    stage1_steps = time_to_steps(stage1_ns, timestep_fs)
    stage2_steps = time_to_steps(stage2_ns, timestep_fs)
    production_steps = time_to_steps(production_ns, timestep_fs)
    log_interval_steps = ps_to_steps(log_interval_ps, timestep_fs)
    traj_interval_steps = ps_to_steps(traj_interval_ps, timestep_fs)

    # Shared kwargs
    def new_integrator():
        return LangevinIntegrator(temperature_k, 1.0 / unit.picosecond, timestep)

    # Minimization + NVT
    system_stage1 = create_system(
        forcefield,
        topology,
        residue_templates,
        with_barostat=False,
        temperature=temperature_k,
        pressure=pressure_atm,
    )
    add_position_restraints(
        system_stage1,
        topology,
        positions,
        1000.0 * unit.kilojoule_per_mole / unit.nanometer**2,
        heavy_only=True,
    )
    integrator1 = new_integrator()
    simulation1 = Simulation(topology, system_stage1, integrator1, platform)
    simulation1.context.setPositions(positions)

    print("[minimize] Performing energy minimization (10000 max iterations)")
    simulation1.minimizeEnergy(maxIterations=minimize_max_steps)
    simulation1.context.setVelocitiesToTemperature(temperature_k)
    run_stage(
        f"NVT ({stage1_ns} ns)",
        simulation1,
        stage1_steps,
        log_path,
        log_interval_steps,
        append_log=False,
    )
    state1 = simulation1.context.getState(getPositions=True, getVelocities=True)
    positions_stage1 = state1.getPositions()

    # NPT equilibration
    system_stage2 = create_system(
        forcefield,
        topology,
        residue_templates,
        with_barostat=True,
        temperature=temperature_k,
        pressure=pressure_atm,
    )
    add_position_restraints(
        system_stage2,
        topology,
        positions_stage1,
        100.0 * unit.kilojoule_per_mole / unit.nanometer**2,
        heavy_only=True,
    )
    integrator2 = new_integrator()
    simulation2 = Simulation(topology, system_stage2, integrator2, platform)
    simulation2.context.setPositions(positions_stage1)
    simulation2.context.setVelocitiesToTemperature(temperature_k)
    run_stage(
        f"NPT ({stage2_ns} ns)",
        simulation2,
        stage2_steps,
        log_path,
        log_interval_steps,
        append_log=True,
    )
    state2 = simulation2.context.getState(getPositions=True, getVelocities=True)
    positions_stage2 = state2.getPositions()

    def run_production_once(
        *,
        timestep_fs_local: float,
        extra_heavy_atom_restraint_k: float | None,
        append_log: bool,
    ):
        local_production_steps = time_to_steps(production_ns, timestep_fs_local)
        local_log_interval_steps = ps_to_steps(log_interval_ps, timestep_fs_local)
        local_traj_interval_steps = ps_to_steps(traj_interval_ps, timestep_fs_local)
        system_prod = create_system(
            forcefield,
            topology,
            residue_templates,
            with_barostat=True,
            temperature=temperature_k,
            pressure=pressure_atm,
        )
        if restrain_zinc:
            add_zinc_restraints(
                system_prod,
                topology,
                positions_stage2,
                zinc_k * unit.kilojoule_per_mole / unit.nanometer**2,
            )
        if extra_heavy_atom_restraint_k is not None:
            add_position_restraints(
                system_prod,
                topology,
                positions_stage2,
                extra_heavy_atom_restraint_k * unit.kilojoule_per_mole / unit.nanometer**2,
                heavy_only=True,
            )
        integrator = LangevinIntegrator(
            temperature_k,
            1.0 / unit.picosecond,
            timestep_fs_local * unit.femtoseconds,
        )
        simulation = Simulation(topology, system_prod, integrator, platform)
        simulation.context.setPositions(positions_stage2)
        simulation.context.setVelocitiesToTemperature(temperature_k)
        traj_reporter = DCDReporter(str(traj_path), local_traj_interval_steps)
        run_stage(
            f"Production ({production_ns} ns)",
            simulation,
            local_production_steps,
            log_path,
            local_log_interval_steps,
            append_log=append_log,
            trajectory_reporter=traj_reporter,
        )
        return simulation.context.getState(getPositions=True)

    # Production with one conservative retry for NaN blowups.
    # This keeps long jobs from fully failing on isolated unstable variants.
    try:
        final_state = run_production_once(
            timestep_fs_local=timestep_fs,
            extra_heavy_atom_restraint_k=None,
            append_log=True,
        )
    except Exception as exc:
        if not _is_nan_coordinate_error(exc):
            raise
        print(
            "[retry] NaN coordinates detected during production. "
            "Retrying with 1.0 fs timestep and light heavy-atom restraints."
        )
        final_state = run_production_once(
            timestep_fs_local=1.0,
            extra_heavy_atom_restraint_k=25.0,
            append_log=True,
        )

    write_final_snapshot(topology, final_state.getPositions(), final_pdb_path)
    print(f"[Done] Outputs in {out_dir}")


def parse_args(argv: Sequence[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pdbs",
        nargs="+",
        type=Path,
        default=None,
        help="Paths to cleaned PDBs (defaults to the four GM-CSF complexes).",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=BASE_DIR / "minimal_md",
        help="Directory to store simulation outputs.",
    )
    parser.add_argument("--platform", type=str, default=None, help="OpenMM platform (CUDA, OpenCL, CPU, etc.)")
    parser.add_argument("--padding-nm", type=float, default=1.2, help="Solvation padding (nm)")
    parser.add_argument("--ionic-strength", type=float, default=0.15, help="Ionic strength for addSolvent (mol/L)")
    parser.add_argument("--temperature", type=float, default=300.0, help="Temperature in Kelvin")
    parser.add_argument("--pressure", type=float, default=1.0, help="Pressure in atm for NPT stages")
    parser.add_argument("--timestep-fs", type=float, default=2.0, help="Integrator timestep in femtoseconds")
    parser.add_argument("--log-interval-ps", type=float, default=1.0, help="StateDataReporter stride in picoseconds")
    parser.add_argument("--traj-interval-ps", type=float, default=2.0, help="Trajectory stride in picoseconds (production)")
    parser.add_argument("--stage1-ns", type=float, default=0.2, help="NVT duration (ns)")
    parser.add_argument("--stage2-ns", type=float, default=0.5, help="NPT equilibration duration (ns)")
    parser.add_argument("--production-ns", type=float, default=10.0, help="Production duration (ns)")
    parser.add_argument("--minimize-max-steps", type=int, default=10000, help="Maximum minimization iterations")
    parser.add_argument(
        "--no-zinc-restraints",
        action="store_true",
        help="Disable Zn–Cys harmonic restraints during production.",
    )
    parser.add_argument(
        "--zinc-k",
        type=float,
        default=200.0,
        help="Zn–Cys restraint constant (kJ/mol/nm^2) when enabled.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None):
    args = parse_args(argv)
    pdbs = args.pdbs or DEFAULT_PDBS
    forcefield = create_forcefield()
    for pdb in pdbs:
        run_pipeline(
            pdb_path=pdb,
            output_root=args.output_root,
            forcefield=forcefield,
            timestep_fs=args.timestep_fs,
            platform_name=args.platform,
            padding_nm=args.padding_nm,
            ionic_strength=args.ionic_strength,
            temperature=args.temperature,
            pressure=args.pressure,
            restrain_zinc=not args.no_zinc_restraints,
            zinc_k=args.zinc_k,
            stage1_ns=args.stage1_ns,
            stage2_ns=args.stage2_ns,
            production_ns=args.production_ns,
            log_interval_ps=args.log_interval_ps,
            traj_interval_ps=args.traj_interval_ps,
            minimize_max_steps=args.minimize_max_steps,
        )


if __name__ == "__main__":
    main()
