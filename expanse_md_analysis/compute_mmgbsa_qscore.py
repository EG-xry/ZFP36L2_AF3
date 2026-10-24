#!/usr/bin/env python3
"""
MM-GBSA interaction energy + native-contact Q-score for all replicates.

1. MM-GBSA (single-trajectory, OpenMM GBn2):
   ΔE_int = E(complex) − E(protein_only) − E(RNA_only)
   Evaluated in implicit solvent (GBn2) on stripped (no water/ion) frames.
   Zn atoms are removed (no GBn2 template); the protein fold is still maintained
   because the trajectory was run WITH Zn restraints — stripping Zn only affects
   the energy evaluation, not the geometry.
   Lower (more negative) ΔE_int = stronger binding.

2. Native-contact Q-score:
   Defines "native contacts" from the WT equilibrated structure (frame 0):
   all protein-heavy ↔ RNA-heavy atom pairs within 4.5 Å at the reference frame.
   Q(t) = fraction of those native contacts present at frame t.
   A deleterious mutation that disrupts the binding pocket should have lower Q.
   This is register-aware because only the canonical WT contacts count.

Run:
  python3 expanse_md_analysis/compute_mmgbsa_qscore.py
"""
from __future__ import annotations
from pathlib import Path
import time
import numpy as np
import pandas as pd
import mdtraj as md

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent

VARIANTS = ["WT", "C206F", "C212A", "D210E", "G204D",
            "R160W", "R184C", "S185N", "T187S", "Y154H"]
TRUTH_DEL = {"Y154H", "R160W", "R184C", "G204D", "C206F", "C212A"}
CLASS = {
    "WT": "WT", "C206F": "deleterious", "C212A": "zn_knockout_control",
    "G204D": "deleterious", "R160W": "deleterious", "R184C": "deleterious",
    "Y154H": "deleterious",
    "D210E": "neutral", "S185N": "neutral", "T187S": "neutral",
}
RNA_RES = {"A", "C", "G", "U", "A5", "A3", "C5", "C3", "G5", "G3", "U5", "U3"}

SOURCES = [
    ("expanse_orig", ROOT / "expanse_md_results"),
    ("anvil_rep1", ROOT / "anvil_md_results" / "rep1"),
    ("anvil_rep2", ROOT / "anvil_md_results" / "rep2"),
    ("expanse_rep2", ROOT / "expanse_md_results_rep2"),
]
MIN_FRAMES = 1800
MAX_NS = 200.0
Q_CUTOFF_NM = 0.45  # 4.5 Å for native contact definition
Q_CUTOFF_EVAL_NM = 0.55  # 5.5 Å for evaluating if contact is "present" (standard 1.2× native cutoff)
MMGBSA_NFRAMES = 50  # subsample for MM-GBSA (CPU-intensive)


def rna_sel():
    return " or ".join(f"resname {r}" for r in RNA_RES)


def load_solute(base_dir, variant):
    """Load protein+RNA+Zn, strip water/ions, cap at MAX_NS."""
    for vname in [variant, "C208A" if variant == "C212A" else variant]:
        d = base_dir / vname
        top = d / "equilibrated.pdb"
        dcd = d / "production.dcd"
        if top.exists() and dcd.exists():
            break
    else:
        return None
    try:
        full_top = md.load_topology(str(top))
        keep = np.unique(np.concatenate([
            full_top.select("protein"),
            full_top.select(rna_sel()),
            full_top.select("name ZN or element Zn"),
        ]))
        traj = md.load(str(dcd), top=str(top), atom_indices=keep)
    except Exception as e:
        print(f"    load error: {e}")
        return None
    nmax = int(MAX_NS * 10) + 1
    if traj.n_frames > nmax:
        traj = traj[:nmax]
    return traj


# ============================================================================
# Q-score: fraction of native protein-RNA contacts
# ============================================================================

def compute_native_contacts(ref_traj):
    """Define native contacts as (protein_resSeq, protein_atomName, rna_resSeq, rna_atomName)
    tuples from frame 0 of the WT reference trajectory. This is topology-independent
    so it can be evaluated on any variant's trajectory regardless of atom count."""
    top = ref_traj.topology
    prot_heavy = top.select("protein and not type H")
    rna_heavy = top.select(f"({rna_sel()}) and not type H")

    pairs = np.array(np.meshgrid(prot_heavy, rna_heavy)).T.reshape(-1, 2)
    d0 = md.compute_distances(ref_traj[0], pairs)[0]
    native_mask = d0 < Q_CUTOFF_NM

    # Store as residue/atom-name tuples (topology-independent)
    atoms = list(top.atoms)
    native_contact_ids = []
    for idx in np.where(native_mask)[0]:
        a1, a2 = atoms[pairs[idx, 0]], atoms[pairs[idx, 1]]
        native_contact_ids.append((
            a1.residue.resSeq, a1.name,
            a2.residue.resSeq, a2.name,
        ))
    print(f"    Q-score: {len(native_contact_ids)} native contacts defined "
          f"(from {len(pairs)} total protein-RNA heavy-atom pairs)")
    return native_contact_ids


def compute_qscore(traj, native_contact_ids, eval_cutoff=Q_CUTOFF_EVAL_NM):
    """Compute Q(t) = fraction of native contacts present at each frame.
    Maps contact IDs (resSeq, atomName) to the current trajectory's atom indices."""
    top = traj.topology
    # Build lookup: (resSeq, atomName) -> atom index
    atom_lookup = {}
    for a in top.atoms:
        atom_lookup[(a.residue.resSeq, a.name)] = a.index

    # Map native contacts to current topology's indices
    pairs = []
    for pres, pname, rres, rname in native_contact_ids:
        idx1 = atom_lookup.get((pres, pname))
        idx2 = atom_lookup.get((rres, rname))
        if idx1 is not None and idx2 is not None:
            pairs.append([idx1, idx2])

    if len(pairs) == 0:
        return np.full(traj.n_frames, np.nan)
    pairs = np.array(pairs)

    # Subsample for speed
    stride = max(1, traj.n_frames // 200)
    sub = traj[::stride]
    d = md.compute_distances(sub, pairs)
    q = (d < eval_cutoff).mean(axis=1)
    return q


# ============================================================================
# MM-GBSA: single-trajectory interaction energy via OpenMM GBn2
# ============================================================================

def compute_mmgbsa(traj, variant, n_frames=MMGBSA_NFRAMES):
    """Single-trajectory MM-GBSA: E(complex) - E(protein) - E(RNA) in GBn2.
    Zn atoms are stripped for the energy evaluation (no GBn2 template).
    Returns array of ΔE_int values (kJ/mol), one per sampled frame."""
    from openmm import unit, Platform, LangevinIntegrator, Context
    from openmm.app import ForceField, NoCutoff, HBonds
    from pdbfixer import PDBFixer
    import io
    from openmm.app import PDBFile as _PDBFile

    # Strip Zn
    zn_mask = np.array([a.element is not None and a.element.symbol == "Zn"
                        for a in traj.topology.atoms])
    if zn_mask.any():
        keep = np.where(~zn_mask)[0]
        traj = traj.atom_slice(keep)

    # Use PDBFixer to add missing terminal H (for RNA 3' end)
    om_top = traj.topology.to_openmm()
    try:
        buf = io.StringIO()
        # Write frame 0 positions in nm (PDBFile expects nm)
        positions_angstrom = traj.xyz[0] * 10.0  # nm → Å for PDB format
        _PDBFile.writeFile(om_top, positions_angstrom, buf)
        buf.seek(0)
        fixer = PDBFixer(pdbfile=buf)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        om_top = fixer.topology
    except Exception as e:
        print(f"    PDBFixer: {e}; using raw topology")

    # Build GBn2 system
    ff = ForceField('amber14/protein.ff14SB.xml', 'amber14/RNA.OL3.xml', 'implicit/gbn2.xml')
    try:
        system = ff.createSystem(om_top, nonbondedMethod=NoCutoff, constraints=HBonds)
    except Exception as e:
        print(f"    MM-GBSA createSystem failed: {e}")
        return np.full(n_frames, np.nan)

    integrator = LangevinIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 2.0 * unit.femtosecond)
    # Use whatever platform is available (Reference is always present)
    platform = None
    for pname in ("CPU", "Reference"):
        try:
            platform = Platform.getPlatformByName(pname)
            break
        except Exception:
            continue
    context = Context(system, integrator, platform) if platform else Context(system, integrator)

    # Identify protein vs RNA atoms in the fixed topology
    prot_idx, rna_idx = [], []
    for atom in om_top.atoms():
        if atom.residue.name.upper() in [r.upper() for r in RNA_RES]:
            rna_idx.append(atom.index)
        else:
            prot_idx.append(atom.index)

    def energy_of(positions):
        context.setPositions(positions)
        return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)

    # Sample frames
    frame_indices = np.linspace(0, traj.n_frames - 1, min(n_frames, traj.n_frames)).astype(int)
    e_ints = []
    for f in frame_indices:
        pos = traj.xyz[f] * unit.nanometer
        try:
            e_complex = energy_of(pos)
            # Protein-only: move RNA far away
            pos_p = traj.xyz[f].copy()
            pos_p[rna_idx] += 500.0
            e_prot = energy_of(pos_p * unit.nanometer)
            # RNA-only: move protein far away
            pos_r = traj.xyz[f].copy()
            pos_r[prot_idx] += 500.0
            e_rna = energy_of(pos_r * unit.nanometer)
            e_ints.append(e_complex - e_prot - e_rna)
        except Exception:
            e_ints.append(np.nan)

    return np.array(e_ints)


# ============================================================================
# Main
# ============================================================================

def auroc(scores, labels):
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    nc = sum(1 for p in pos for n_ in neg if p > n_)
    nt = sum(1 for p in pos for n_ in neg if p == n_)
    return (nc + 0.5 * nt) / (len(pos) * len(neg))


def main():
    # Load WT reference for Q-score native contacts
    print("Loading WT reference for native-contact definition...")
    wt_ref = load_solute(ROOT / "expanse_md_results", "WT")
    if wt_ref is None:
        print("ERROR: cannot load WT reference")
        return
    native_pairs = compute_native_contacts(wt_ref)
    del wt_ref  # free memory

    all_rows = []
    for src_name, base_dir in SOURCES:
        print(f"\n{'='*60}")
        print(f"SOURCE: {src_name}")
        print(f"{'='*60}")
        for v in VARIANTS:
            dcd_path = None
            for vname in [v, "C208A" if v == "C212A" else v]:
                p = base_dir / vname / "production.dcd"
                if p.exists():
                    dcd_path = p
                    break
            if dcd_path is None:
                continue

            print(f"  {v}: loading...", end=" ", flush=True)
            traj = load_solute(base_dir, v)
            if traj is None or traj.n_frames < MIN_FRAMES:
                print(f"skip ({traj.n_frames if traj else 0} frames)")
                continue
            print(f"{traj.n_frames} frames", end="", flush=True)

            row = {"variant": v, "class": CLASS[v], "source": src_name,
                   "n_frames": traj.n_frames}

            # Q-score
            print(" | Q", end="", flush=True)
            q = compute_qscore(traj, native_pairs)
            last20 = max(1, int(0.2 * len(q)))
            row["qscore_mean"] = float(np.nanmean(q))
            row["qscore_last20pct"] = float(np.nanmean(q[-last20:]))

            # MM-GBSA
            print(" | GBSA", end="", flush=True)
            t0 = time.time()
            e_int = compute_mmgbsa(traj, v, n_frames=30)
            dt = time.time() - t0
            row["mmgbsa_mean_kJmol"] = float(np.nanmean(e_int))
            row["mmgbsa_std_kJmol"] = float(np.nanstd(e_int))

            print(f" | {dt:.0f}s | Q={row['qscore_last20pct']:.3f} GBSA={row['mmgbsa_mean_kJmol']:.0f}")
            all_rows.append(row)

    df = pd.DataFrame(all_rows)
    df.to_csv(OUT / "mmgbsa_qscore_all.csv", index=False)
    print(f"\nWrote {OUT / 'mmgbsa_qscore_all.csv'}")

    # Replicate-averaged values per variant
    metric_cols = ["qscore_mean", "qscore_last20pct", "mmgbsa_mean_kJmol"]
    avg = df.groupby("variant")[metric_cols].mean().reset_index()
    avg["class"] = avg["variant"].map(CLASS)
    rep_counts = df.groupby("variant").size().reset_index(name="n_reps")
    avg = avg.merge(rep_counts, on="variant")

    print(f"\n{'='*60}")
    print("REPLICATE-AVERAGED VALUES")
    print(f"{'='*60}")
    for _, r in avg.sort_values("qscore_last20pct").iterrows():
        tag = 'del' if r['variant'] in TRUTH_DEL else ('neu' if r['variant'] in ('D210E','S185N','T187S') else 'WT')
        print(f"  {r['variant']:8s} ({tag}) Q={r['qscore_last20pct']:.3f}  "
              f"GBSA={r['mmgbsa_mean_kJmol']:.0f}  ({int(r['n_reps'])} reps)")

    # AUROC (all 9 variants: 6 del vs 3 neu)
    print(f"\n{'='*60}")
    print("AUROC (6 del vs 3 neu, including C212A)")
    print(f"{'='*60}")
    mut = avg[avg["variant"] != "WT"]
    y = mut["variant"].apply(lambda v: 1 if v in TRUTH_DEL else 0).values
    for feat, sign, label in [
        ("qscore_last20pct", -1, "Q-score last20% (lower=more del)"),
        ("qscore_mean", -1, "Q-score mean (lower=more del)"),
        ("mmgbsa_mean_kJmol", +1, "MM-GBSA (less negative=weaker=more del)"),
    ]:
        vals = mut[feat].values.astype(float) * sign
        if np.isfinite(vals).all():
            a = auroc(vals, y)
            print(f"  {label:45s} AUROC = {a:.3f}")

    # AUROC without C212A (5 del vs 3 neu)
    print(f"\n{'='*60}")
    print("AUROC (5 del vs 3 neu, excluding C212A)")
    print(f"{'='*60}")
    TRUTH_5 = {"Y154H", "R160W", "R184C", "G204D", "C206F"}
    mut5 = avg[(avg["variant"] != "WT") & (avg["variant"] != "C212A")]
    y5 = mut5["variant"].apply(lambda v: 1 if v in TRUTH_5 else 0).values
    for feat, sign, label in [
        ("qscore_last20pct", -1, "Q-score last20% (lower=more del)"),
        ("qscore_mean", -1, "Q-score mean (lower=more del)"),
        ("mmgbsa_mean_kJmol", +1, "MM-GBSA (less negative=weaker=more del)"),
    ]:
        vals = mut5[feat].values.astype(float) * sign
        if np.isfinite(vals).all():
            a = auroc(vals, y5)
            print(f"  {label:45s} AUROC = {a:.3f}")

    avg.to_csv(OUT / "mmgbsa_qscore_averaged.csv", index=False)
    print(f"\nWrote {OUT / 'mmgbsa_qscore_averaged.csv'}")


if __name__ == "__main__":
    main()
