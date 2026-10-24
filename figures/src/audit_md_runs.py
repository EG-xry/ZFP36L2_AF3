#!/usr/bin/env python3
"""
Audit what every MD run actually simulated.

For each production run this reads the run's own topology (equilibrated.pdb),
compares the protein sequence with the wild-type construct, and records the
closest MD input model (anvil_md/best_models/) by Ca RMSD. The sequence check
is decisive; the closest-model column is only supporting evidence, because the
AF3 models of all variants lie within ~1 Å of each other. The audit exists
because a folder name is not proof of the system: every Anvil replicate in a
variant folder turned out to contain the wild-type protein.

The topologies live in the trajectory folders, which are too large for the
repository, so this runs once locally and writes a small CSV that
make_figures.py reads.

Output: figures/data/md_run_inventory.csv
    source, folder, protein, closest_input_model, closest_rmsd_A, rna, zn

Env: /Users/eric/micromamba/envs/omm/bin/python  (needs mdtraj)
"""
from __future__ import annotations

from pathlib import Path

import mdtraj as md
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures" / "data" / "md_run_inventory.csv"
SOURCES = [
    ("expanse_orig", ROOT / "expanse_md_results"),
    ("anvil_rep1", ROOT / "anvil_md_results" / "rep1"),
    ("anvil_rep2", ROOT / "anvil_md_results" / "rep2"),
    ("expanse_rep2", ROOT / "expanse_md_results_rep2"),
]
WT_SEQ = "STRYKTELCRPFEESGTCKYGEKCQFAHGFHELRSLTRHPKYKTELCRTFHTIGFCPYGPRCHFIHNADE"
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
RNA_NAMES = {"A", "C", "G", "U", "A5", "A3", "C5", "C3", "G5", "G3", "U5", "U3"}


def ca_only(traj):
    return traj.atom_slice(traj.topology.select("protein and name CA"))


def protein_label(top) -> str:
    seq = "".join(AA3.get(r.name, "?") for r in top.residues if r.is_protein)
    diffs = [f"{w}{i + 151}{m}" for i, (w, m) in enumerate(zip(WT_SEQ, seq)) if w != m]
    return "+".join(diffs) if diffs else "WT"


def main() -> None:
    models = {p.name.split("_")[0]: ca_only(md.load(str(p)))
              for p in sorted((ROOT / "anvil_md" / "best_models").glob("*_best_model.pdb"))}
    rows = []
    for source, base in SOURCES:
        for pdb in sorted(base.glob("*/equilibrated.pdb")):
            traj = md.load(str(pdb))
            top = traj.topology
            rmsd = {k: float(md.rmsd(ca_only(traj), m, 0)[0] * 10) for k, m in models.items()}
            start = min(rmsd, key=rmsd.get)
            rows.append({
                "source": source,
                "folder": pdb.parent.name,
                "protein": protein_label(top),
                "closest_input_model": start,
                "closest_rmsd_A": round(rmsd[start], 2),
                "rna": "".join(r.name.rstrip("35")[-1] for r in top.residues
                               if r.name.upper() in RNA_NAMES),
                "zn": sum(1 for a in top.atoms if a.element is not None
                          and a.element.symbol == "Zn"),
            })
            print(f"  {source:13s} {pdb.parent.name:6s} protein={rows[-1]['protein']:7s} "
                  f"start={start} ({rmsd[start]:.2f} Å)", flush=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
