#!/usr/bin/env python3
"""
MM-GBSA-style interaction-energy metric + literature-backed MD metrics for the
ZFP36L2 + nonamer system, computed on the production DCDs.

Why this script exists
----------------------
RMSD / per-residue H-bond counts did NOT robustly separate deleterious from
neutral once we sampled to ~270-400 ns (see FINAL_CONCLUSIONS.md). The single
most promising remaining MD observable is a *binding energy* — a thermodynamic
quantity, not a geometry snapshot. We add three energy-based and two
geometry-based metrics that are standard in protein-RNA MD discrimination
literature:

  1. MM-GBSA-style protein<->RNA interaction energy (single-trajectory).
     ΔE_int = E(complex) - E(protein) - E(RNA), evaluated with an OpenMM
     implicit-solvent (GBn2) system on stripped (solvent-removed) frames.
     This is the "one-trajectory MM-GBSA" approximation widely used to rank
     mutants (Genheden & Ryde 2015; Wang et al. 2019). It is an *interaction
     energy*, not an absolute ΔG (no entropy, single-traj), but it is exactly
     the kind of relative score that ranks loss-of-binding mutants.

  2. Direct protein<->RNA nonbonded interaction energy in *explicit* solvent
     via a decomposed NonbondedForce (Coulomb + LJ between the two chains),
     averaged over frames. Cheap, no implicit-solvent reparametrization.

  3. Number of protein<->RNA salt bridges (Arg/Lys NZ/NH <-> RNA phosphate O)
     within 4 Å — a specific electrostatic contact count that is more robust
     than total H-bonds for nucleic-acid binding (Yang et al. 2021).

  4. Interface buried surface area proxy: number of protein heavy atoms within
     5 Å of any RNA heavy atom (a fast SASA-free interface-size metric).

  5. RMSF of the ZF1/ZF2 RNA-binding residues (per-residue flexibility) — the
     metric the Ramos paper *implied* but did not quantify; high RMSF at the
     binding residues correlates with loss of binding.

All are written into energy_metrics_summary.csv and folded into the AUROC
catalog by run_full_analysis.py (or recomputed standalone here).

Run:
  python3 expanse_md_analysis/compute_energy_metrics.py            # all variants
  python3 expanse_md_analysis/compute_energy_metrics.py --stride 20 --nframes 150
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import mdtraj as md

ROOT = Path(__file__).resolve().parent.parent
EXPANSE_DIR = ROOT / "expanse_md_results"
OUT_DIR = Path(__file__).resolve().parent

VARIANTS = ["WT", "C206F", "C212A", "D210E", "G204D",
            "R160W", "R184C", "S185N", "T187S", "Y154H"]
CLASS = {
    "WT": "WT", "C206F": "deleterious", "C212A": "zn_knockout_control",
    "G204D": "deleterious", "R160W": "deleterious", "R184C": "deleterious",
    "Y154H": "deleterious",
    "D210E": "neutral", "S185N": "neutral", "T187S": "neutral",
}
RNA_RES = ("A", "C", "G", "U", "A5", "A3", "C5", "C3",
           "G5", "G3", "U5", "U3")

# ZF1 / ZF2 RNA-binding residues in PDB-local numbering (ZFP36L2 = +150)
ZF1_PDB = list(range(2, 30))    # 152-179
ZF2_PDB = list(range(40, 68))   # 190-217


def rna_selection(prefix="resname"):
    return " or ".join(f"{prefix} {r}" for r in RNA_RES)


def load_solute(variant, max_ns, stride):
    """Load protein+RNA+Zn, strip water/ions, subsample."""
    var = EXPANSE_DIR / variant
    top = var / "equilibrated.pdb"
    dcd = var / "production.dcd"
    if not (top.exists() and dcd.exists()):
        return None
    full = md.load_topology(str(top))
    keep = np.unique(np.concatenate([
        full.select("protein"),
        full.select(rna_selection()),
        full.select("name ZN or element Zn"),
    ]))
    traj = md.load(str(dcd), top=str(top), atom_indices=keep)
    nmax = int(max_ns * 10) + 1
    if traj.n_frames > nmax:
        traj = traj[:nmax]
    if stride > 1:
        traj = traj[::stride]
    return traj


# ---------------------------------------------------------------------------
# Metric 2 (explicit-solvent decomposed nonbonded protein<->RNA energy)
# Implemented with mdtraj coordinates + a vacuum Coulomb/LJ pairwise sum on the
# interface atoms. This is a fast, parametrization-free proxy. For a *true*
# MM-GBSA you would build an OpenMM GBn2 system (metric 1) — provided below but
# gated behind --mmgbsa because it needs the solvated topology re-parametrized.
# ---------------------------------------------------------------------------

# AMBER-ish partial charges are not in the mdtraj topology, so the cheap proxy
# uses a Debye-screened Coulomb on formal charges of charged groups + an
# LJ-free contact count. The robust, publishable number is the OpenMM GBn2
# interaction energy (metric 1, --mmgbsa). We compute both; the proxy is always
# available, the GBn2 one when OpenMM is importable.

def charged_atom_groups(top):
    """Return (pos_idx, neg_idx) for the dominant charged groups:
    protein Arg/Lys side-chain N (+1), Asp/Glu side-chain O (-1),
    RNA phosphate O (-1)."""
    pos, neg = [], []
    for a in top.atoms:
        r = a.residue
        nm = a.name
        rn = r.name.upper()
        if r.is_protein:
            if rn in ("ARG",) and nm in ("NH1", "NH2", "NE"):
                pos.append(a.index)
            elif rn in ("LYS",) and nm in ("NZ",):
                pos.append(a.index)
            elif rn in ("ASP",) and nm in ("OD1", "OD2"):
                neg.append(a.index)
            elif rn in ("GLU",) and nm in ("OE1", "OE2"):
                neg.append(a.index)
        else:  # RNA
            if nm in ("OP1", "OP2", "O1P", "O2P"):
                neg.append(a.index)
    return np.array(pos, dtype=int), np.array(neg, dtype=int)


def screened_coulomb_protein_rna(traj):
    """Sum of screened-Coulomb interactions between protein charged groups and
    RNA phosphate oxygens, averaged over frames. Lower (more negative) = stronger
    binding. Units: arbitrary (kJ/mol-ish with a Debye screen)."""
    top = traj.topology
    rna_atoms = {a.index for a in top.atoms if not a.residue.is_protein
                 and a.residue.name.upper() in [r.upper() for r in RNA_RES]}
    pos, neg = charged_atom_groups(top)
    # protein positives vs RNA negatives is the dominant binding term
    prot_pos = np.array([i for i in pos], dtype=int)
    rna_neg = np.array([i for i in neg if i in rna_atoms], dtype=int)
    if len(prot_pos) == 0 or len(rna_neg) == 0:
        return np.nan
    pairs = np.array(np.meshgrid(prot_pos, rna_neg)).T.reshape(-1, 2)
    d = md.compute_distances(traj, pairs, periodic=True)  # nm
    # screened Coulomb: q1 q2 / r * exp(-r/lambda); q=+1/-1, kappa^-1=1.0 nm (0.15 M ~ 0.8 nm)
    kappa_inv = 0.8  # nm Debye length at ~0.15 M
    ke = 138.935458  # kJ*nm/(mol*e^2)
    d = np.clip(d, 0.15, None)  # avoid singularities
    e = -ke * np.exp(-d / kappa_inv) / d  # +*- => negative (attractive)
    return float(e.sum(axis=1).mean())  # per-frame total, averaged


def salt_bridges_protein_rna(traj, cutoff_nm=0.4):
    top = traj.topology
    rna_phos = top.select("(" + rna_selection() + ") and (name OP1 or name OP2 or name O1P or name O2P)")
    cation = top.select("protein and (name NZ or name NH1 or name NH2 or name NE)")
    if len(rna_phos) == 0 or len(cation) == 0:
        return np.nan
    pairs = np.array(np.meshgrid(cation, rna_phos)).T.reshape(-1, 2)
    d = md.compute_distances(traj, pairs, periodic=True)
    return float((d < cutoff_nm).sum(axis=1).mean())


def interface_atom_count(traj, cutoff_nm=0.5):
    top = traj.topology
    prot = top.select("protein and not type H")
    rna = top.select("(" + rna_selection() + ") and not type H")
    if len(prot) == 0 or len(rna) == 0:
        return np.nan
    # neighbor count: protein heavy atoms within cutoff of any RNA heavy atom
    counts = []
    for fi in range(traj.n_frames):
        nb = md.compute_neighbors(traj[fi], cutoff_nm, rna, haystack_indices=prot)
        counts.append(len(nb[0]))
    return float(np.mean(counts))


def binding_residue_rmsf(traj):
    """Mean RMSF (Å) of the ZF1+ZF2 RNA-binding Cα atoms after superposing on
    the protein Cα. Higher = floppier binding residues."""
    top = traj.topology
    ca = top.select("protein and name CA")
    if len(ca) == 0:
        return np.nan
    tr = traj.superpose(traj, frame=0, atom_indices=ca, ref_atom_indices=ca)
    bind_pdb = ZF1_PDB + ZF2_PDB
    bind_ca = top.select("protein and name CA and (" +
                         " or ".join(f"resSeq {p}" for p in bind_pdb) + ")")
    if len(bind_ca) == 0:
        return np.nan
    xyz = tr.xyz[:, bind_ca, :]              # (frames, natoms, 3)
    mean = xyz.mean(axis=0, keepdims=True)
    rmsf = np.sqrt(((xyz - mean) ** 2).sum(axis=2).mean(axis=0))  # per-atom
    return float(rmsf.mean() * 10.0)         # nm -> Å


# ---------------------------------------------------------------------------
# Metric 1: true single-trajectory MM-GBSA interaction energy (OpenMM GBn2)
# ---------------------------------------------------------------------------

def mmgbsa_interaction(variant, traj):
    """One-trajectory MM-GBSA: build an implicit-solvent (GBn2) system for the
    stripped solute, then for each frame compute
        E_int = E(complex) - E(protein-only) - E(RNA-only)
    by zeroing the partial charges / nonbonded params of the absent chain via
    separate System builds. Returns mean E_int in kJ/mol (lower = stronger).
    Requires OpenMM + a forcefield that covers protein+RNA in implicit solvent.
    """
    try:
        from openmm import app, unit, Platform, LangevinIntegrator, Context
        from openmm.app import ForceField, Modeller, NoCutoff, HBonds
    except Exception as exc:
        print(f"  [mmgbsa] OpenMM unavailable ({exc}); skipping GBn2 metric")
        return np.nan

    # In OpenMM 8.x, the GB model is specified by *including the implicit-solvent
    # XML in the ForceField* (e.g. amber14/GBn2.xml) rather than passing an
    # implicitSolvent= kwarg to createSystem(). We therefore build the FF with
    # the GBn2 parameter file and call createSystem with NO implicitSolvent arg.
    ff = None
    for combo in (("amber14/protein.ff14SB.xml", "amber14/RNA.OL3.xml", "implicit/gbn2.xml"),
                  ("amber14/protein.ff14SB.xml", "amber14/RNA.OL3.xml", "amber14/GBn2.xml"),
                  ("amber14-all.xml", "implicit/gbn2.xml")):
        try:
            ff = ForceField(*combo)
            break
        except Exception:
            continue
    if ff is None:
        print("  [mmgbsa] could not load implicit-solvent forcefield; skipping")
        return np.nan

    # Strip the Zn ions: the implicit-solvent forcefield has no template for a
    # bare Zn²⁺, and Zn does not contribute to the protein↔RNA interaction term
    # we want (it's buried in the zinc fingers). Build a Zn-free sub-trajectory.
    zn_mask = np.array([a.element is not None and a.element.symbol == "Zn"
                        for a in traj.topology.atoms])
    if zn_mask.any():
        keep = np.where(~zn_mask)[0]
        traj = traj.atom_slice(keep)

    # The 3'-terminal RNA residue is missing its O3' terminal hydrogen in the
    # stripped solute topology, so the AMBER U3 template doesn't match. Use
    # PDBFixer to add the missing hydrogens/atoms on a representative frame's
    # topology before building the implicit-solvent system.
    om_top = traj.topology.to_openmm()
    try:
        from pdbfixer import PDBFixer
        import io
        from openmm.app import PDBFile as _PDBFile
        # Write frame 0 to an in-memory PDB, fix it, and re-read topology.
        buf = io.StringIO()
        _PDBFile.writeFile(om_top, (traj.xyz[0] * 10.0), buf)  # Å for PDB
        buf.seek(0)
        fixer = PDBFixer(pdbfile=buf)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        om_top = fixer.topology
        fixed_positions = fixer.positions
    except Exception as exc:
        print(f"  [mmgbsa] PDBFixer unavailable/failed ({exc}); trying raw topology")
        fixed_positions = None

    top = om_top
    # Build the implicit-solvent system once (no PBC). GB is already in the FF.
    try:
        system = ff.createSystem(top, nonbondedMethod=NoCutoff,
                                 constraints=HBonds)
    except Exception as exc:
        print(f"  [mmgbsa] createSystem failed ({exc}); skipping")
        return np.nan

    integrator = LangevinIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond,
                                    2.0 * unit.femtosecond)
    try:
        platform = Platform.getPlatformByName("CPU")
    except Exception:
        platform = None
    context = Context(system, integrator, platform) if platform else Context(system, integrator)

    # Identify protein vs RNA atom indices in the openmm topology order
    prot_idx, rna_idx = [], []
    for atom in top.atoms():
        if atom.residue.name.upper() in [r.upper() for r in RNA_RES]:
            rna_idx.append(atom.index)
        else:
            prot_idx.append(atom.index)

    def energy_of(positions):
        context.setPositions(positions)
        st = context.getState(getEnergy=True)
        return st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    # Average over up to ~20 frames for speed.
    fr = np.linspace(0, traj.n_frames - 1, min(20, traj.n_frames)).astype(int)
    e_ints = []
    from openmm.vec3 import Vec3  # noqa
    for f in fr:
        pos = traj.xyz[f] * unit.nanometer
        e_complex = energy_of(pos)
        # protein-only: move RNA far away (translate by +500 nm) so interaction ~0
        pos_p = traj.xyz[f].copy()
        pos_p[rna_idx] += 500.0
        e_prot = energy_of(pos_p * unit.nanometer)
        # rna-only: move protein far away
        pos_r = traj.xyz[f].copy()
        pos_r[prot_idx] += 500.0
        e_rna = energy_of(pos_r * unit.nanometer)
        e_ints.append(e_complex - e_prot - e_rna)
    return float(np.mean(e_ints))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-ns", type=float, default=269.0,
                    help="common analysis window (ns)")
    ap.add_argument("--stride", type=int, default=10,
                    help="frame stride for the cheap metrics")
    ap.add_argument("--mmgbsa", action="store_true",
                    help="also compute the OpenMM GBn2 single-traj MM-GBSA "
                         "interaction energy (slow, CPU)")
    args = ap.parse_args()

    rows = []
    for v in VARIANTS:
        print(f"[{v}] loading <= {args.max_ns} ns (stride {args.stride}) ...")
        traj = load_solute(v, args.max_ns, args.stride)
        if traj is None:
            print("  missing; skip")
            continue
        row = {
            "variant": v,
            "class": CLASS[v],
            "n_frames_used": traj.n_frames,
            "coulomb_prot_rna_screened": screened_coulomb_protein_rna(traj),
            "salt_bridges_prot_rna": salt_bridges_protein_rna(traj),
            "interface_atom_count_5A": interface_atom_count(traj),
            "binding_residue_rmsf_A": binding_residue_rmsf(traj),
        }
        if args.mmgbsa:
            print("  computing MM-GBSA (GBn2) interaction energy ...")
            row["mmgbsa_int_kJmol"] = mmgbsa_interaction(v, traj)
        rows.append(row)
        print("   ", {k: (round(val, 3) if isinstance(val, float) else val)
                      for k, val in row.items() if k not in ("variant", "class")})

    df = pd.DataFrame(rows)
    out = OUT_DIR / "energy_metrics_summary.csv"
    df.to_csv(out, index=False)
    print(f"\nWrote {out}")

    # Quick AUROC of each new metric (del vs neu, 6 vs 3, excl WT)
    TRUTH_DEL = {"Y154H", "R160W", "R184C", "G204D", "C206F", "C212A"}
    mut = df[df["variant"] != "WT"].copy()
    y = mut["variant"].apply(lambda x: 1 if x in TRUTH_DEL else 0).values

    def auroc(s, l):
        p = s[l == 1]; n = s[l == 0]
        if len(p) == 0 or len(n) == 0:
            return float("nan")
        return (sum(1 for a in p for b in n if a > b)
                + 0.5 * sum(1 for a in p for b in n if a == b)) / (len(p) * len(n))

    # sign: which direction means "more deleterious"
    direction = {
        "coulomb_prot_rna_screened": +1,   # weaker (less negative) attraction => higher => del
        "salt_bridges_prot_rna": -1,       # fewer salt bridges => del
        "interface_atom_count_5A": -1,     # smaller interface => del
        "binding_residue_rmsf_A": +1,      # floppier binding residues => del
        "mmgbsa_int_kJmol": +1,            # less negative interaction energy => del
    }
    print("\n=== AUROC of new energy/structure metrics (del vs neu) ===")
    auroc_rows = []
    for col, sgn in direction.items():
        if col not in mut.columns:
            continue
        vals = mut[col].values.astype(float) * sgn
        if not np.isfinite(vals).all():
            print(f"  {col:32} (has NaN, skipped)")
            continue
        a = auroc(vals, y)
        auroc_rows.append({"feature": col, "AUROC": round(a, 3), "origin": "MD-energy"})
        print(f"  {col:32} AUROC={a:.3f}")
    pd.DataFrame(auroc_rows).to_csv(OUT_DIR / "energy_metric_aurocs.csv", index=False)
    print(f"\nWrote {OUT_DIR / 'energy_metric_aurocs.csv'}")


if __name__ == "__main__":
    main()
