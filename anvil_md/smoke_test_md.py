#!/usr/bin/env python3
"""Preflight smoke test: clean PDB, minimize, short NVT on CPU/CUDA."""

from __future__ import annotations

import argparse
import sys
import tempfile
from pathlib import Path

import numpy as np

from openmm import LangevinIntegrator, Platform, unit
from openmm.app import Simulation

# Import shared helpers from the main driver
sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_md_anvil as md
from clean_af3_structures import clean_pdb


def resolve_pdb(pdb_path: Path) -> Path:
    """Return a cleaned PDB path, cleaning raw AF3 input if needed."""
    if "_cleaned" in pdb_path.stem:
        return pdb_path.resolve()
    cleaned = Path(tempfile.mkdtemp()) / f"{pdb_path.stem}_cleaned.pdb"
    clean_pdb(pdb_path.resolve(), cleaned)
    return cleaned


def check_zn_distances(topology, positions, cutoff_nm=0.35):
    positions_nm = positions.value_in_unit(unit.nanometer)
    zn_indices = [
        a.index for a in topology.atoms()
        if a.element and a.element.symbol == "Zn"
    ]
    cys_sg = [
        a.index for a in topology.atoms()
        if a.residue.name == "CYS" and a.name == "SG"
    ]
    print(f"[Zn] {len(zn_indices)} Zn ions, {len(cys_sg)} Cys SG atoms")
    for zn_idx in zn_indices:
        zn_pos = np.array(positions_nm[zn_idx])
        nearby = []
        for sg_idx in cys_sg:
            sg_pos = np.array(positions_nm[sg_idx])
            dist = np.linalg.norm(zn_pos - sg_pos)
            if dist <= cutoff_nm:
                nearby.append(dist)
        print(f"  Zn atom {zn_idx}: {len(nearby)} Cys SG within {cutoff_nm} nm "
              f"(distances: {[f'{d:.3f}' for d in sorted(nearby)]})")


def run_smoke(pdb_path: Path, platform_name: str, precision: str, nvt_steps: int):
    cleaned = resolve_pdb(pdb_path)
    print(f"[Input] {cleaned}")

    ff = md.create_forcefield()
    modeller, template_ids = md.build_modeller(cleaned, ff, padding_nm=1.2, ionic_strength=0.15)
    topology = modeller.topology
    positions = modeller.positions
    res_templates = md.build_residue_template_map(topology, template_ids)

    check_zn_distances(topology, positions)

    temperature = 310.0 * unit.kelvin
    pressure = 1.0 * unit.atmospheres
    timestep = 2.0 * unit.femtoseconds

    system = md.create_system(ff, topology, res_templates, False, temperature, pressure)
    md.add_position_restraints(
        system, topology, positions,
        1000.0 * unit.kilojoule_per_mole / unit.nanometer**2,
    )

    platform, properties = md.get_platform_and_properties(platform_name, precision)
    integrator = LangevinIntegrator(temperature, 1.0 / unit.picosecond, timestep)
    sim = Simulation(topology, system, integrator, platform, properties)
    sim.context.setPositions(positions)

    print("[Minimize] 200 iterations...")
    sim.minimizeEnergy(maxIterations=200)
    e_min = sim.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"  Energy after minimize: {e_min}")

    sim.context.setVelocitiesToTemperature(temperature)
    print(f"[NVT] {nvt_steps} steps on {platform_name} ({precision})...")
    sim.step(nvt_steps)

    state = sim.context.getState(getPositions=True, getEnergy=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    has_nan = bool(np.isnan(pos).any())
    e = state.getPotentialEnergy()
    temp = sim.context.getState(getEnergy=True).getKineticEnergy()

    print(f"  Potential energy: {e}")
    print(f"  NaN in coordinates: {has_nan}")

    if has_nan:
        print("\nSMOKE TEST FAILED: NaN coordinates detected")
        sys.exit(1)

    print("\nSMOKE TEST PASSED")
    return 0


def main():
    parser = argparse.ArgumentParser(description="MD preflight smoke test")
    parser.add_argument("--pdb", type=Path, required=True, help="Raw or cleaned input PDB")
    parser.add_argument("--platform", type=str, default="CPU", help="OpenMM platform")
    parser.add_argument("--precision", type=str, default="single",
                        choices=["single", "double", "mixed"])
    parser.add_argument("--nvt-steps", type=int, default=500, help="NVT steps to run")
    args = parser.parse_args()

    if not args.pdb.exists():
        print(f"ERROR: PDB not found: {args.pdb}")
        sys.exit(2)

    run_smoke(args.pdb, args.platform, args.precision, args.nvt_steps)


if __name__ == "__main__":
    main()
