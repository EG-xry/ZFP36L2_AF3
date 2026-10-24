#!/usr/bin/env python3
"""
Compute trajectory-based metrics for each variant's production.dcd:
  - Cα RMSD vs starting (equilibrated.pdb) — full backbone stability
  - Cα RMSD vs starting on ZF1 (res 152–179) and ZF2 (res 190–217)
  - RNA P-atom RMSD vs starting — RNA placement stability
  - Min P-atom distance to nearest Cα — interface contact tracking
  - Number of inter-chain (protein–RNA) heavy-atom contacts < 4 Å

Robust to common AF3 + OpenMM topology where:
  - Chain A = ZFP36L2 protein (~70 aa, includes the zinc finger domain)
  - Chain B = 9-nt RNA nonamer
  - Plus Zn ions, ions, waters

Produces:
  - <variant>_rmsd.csv  per-variant time series
  - rmsd_summary.csv    final-vs-start + average summaries
  - rmsd_summary.md
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import mdtraj as md

EXPANSE_DIR = Path(__file__).resolve().parent.parent / "expanse_md_results"
OUT_DIR = Path(__file__).resolve().parent
OUT_DIR.mkdir(parents=True, exist_ok=True)

VARIANTS = ["WT", "C206F", "C212A", "D210E", "G204D",
            "R160W", "R184C", "S185N", "T187S", "Y154H"]

CLASS_MAP = {
    "WT": "WT",
    "C206F": "deleterious",
    "C212A": "zn_knockout_control",  # was mislabeled C208A; mutation is at C212
    "G204D": "deleterious",
    "R160W": "deleterious",
    "R184C": "deleterious",
    "Y154H": "deleterious",
    "D210E": "neutral",  # NOTE: actual mutation in this run is D219E
    "S185N": "neutral",
    "T187S": "neutral",
}

# AF3 numbering for ZFP36L2 zinc finger 1 / 2 (approximate)
# NOTE: AF3 PDB outputs use 1..70 numbering; ZFP36L2_resnum = PDB_resSeq + 150
# ZF1 (ZFP36L2 152-179) -> PDB 2-29
# ZF2 (ZFP36L2 190-217) -> PDB 40-67
ZF1_RANGE_PDB = (2, 29)
ZF2_RANGE_PDB = (40, 67)


def load_traj(variant: str):
    var_dir = EXPANSE_DIR / variant
    top = var_dir / "equilibrated.pdb"
    dcd = var_dir / "production.dcd"
    if not (top.exists() and dcd.exists()):
        return None, None
    print(f"  loading {dcd.name} ({dcd.stat().st_size / 1e6:.0f} MB)...")
    # Strip waters/ions on the fly to save memory: load with topology then atom_slice
    top_t = md.load_topology(str(top))
    # Identify solute atoms = protein + RNA + Zn
    protein_idx = top_t.select("protein")
    rna_idx = top_t.select("resname A or resname C or resname G or resname U or "
                           "resname A5 or resname A3 or resname C5 or resname C3 or "
                           "resname G5 or resname G3 or resname U5 or resname U3")
    zn_idx = top_t.select("name ZN or element Zn")
    keep = np.unique(np.concatenate([protein_idx, rna_idx, zn_idx]))
    traj = md.load(str(dcd), top=str(top), atom_indices=keep)
    sub_top = traj.topology
    return traj, sub_top


def select_atoms(top, sel: str) -> np.ndarray:
    return top.select(sel)


def compute_rmsd_series(traj, ref_idx, target_idx):
    """RMSD over target_idx after fitting on ref_idx (both per-frame)."""
    # mdtraj.superpose & rmsd are simplest:
    traj_aligned = traj.superpose(traj, frame=0, atom_indices=ref_idx, ref_atom_indices=ref_idx)
    diff = traj_aligned.xyz[:, target_idx, :] - traj_aligned.xyz[0:1, target_idx, :]
    return np.sqrt((diff ** 2).sum(axis=2).mean(axis=1)) * 10.0  # nm -> Å


def interface_contacts(traj, prot_heavy_idx, rna_heavy_idx, cutoff_nm=0.4, stride=10):
    """Number of protein-heavy / RNA-heavy atom pairs within cutoff per frame.
    To save memory & time, subsample frames with `stride`."""
    n_frames = traj.n_frames
    frame_idxs = np.arange(0, n_frames, stride)
    out = np.zeros(len(frame_idxs), dtype=int)
    # Build the cross pair index for distance calculation
    pairs = np.array(np.meshgrid(prot_heavy_idx, rna_heavy_idx)).T.reshape(-1, 2)
    sub = traj[frame_idxs]
    dists = md.compute_distances(sub, pairs, periodic=True)  # nm
    in_cut = dists < cutoff_nm
    out = in_cut.sum(axis=1)
    return frame_idxs, out


def main():
    summary_rows = []
    for variant in VARIANTS:
        print(f"\n[{variant}]")
        traj, top = load_traj(variant)
        if traj is None:
            print("  missing input; skipping")
            continue

        # Time axis (production.dcd traj_interval defaulted to 100 ps in run script)
        # We can derive from traj.time if present, else fall back to 100 ps spacing.
        if traj.time is not None and traj.time[-1] > 0:
            time_ns = traj.time / 1000.0
        else:
            time_ns = np.arange(traj.n_frames) * 0.1  # 100 ps per frame

        # Atom selections
        ca_all = top.select("name CA and protein")
        rna_p = top.select("(resname A or resname C or resname G or resname U or "
                           "resname A5 or resname A3 or resname C5 or resname C3 or "
                           "resname G5 or resname G3 or resname U5 or resname U3) and name P")
        if len(ca_all) == 0:
            print("  WARNING no CA atoms found; skip")
            continue
        if len(rna_p) == 0:
            print("  WARNING no RNA P atoms found; skip")

        # ZF1 / ZF2 selections by residue id (AF3 uses sequence numbering)
        residues = list(top.residues)
        protein_resids = sorted({r.resSeq for r in residues if r.is_protein})
        if not protein_resids:
            print("  WARNING no protein residues; skip")
            continue
        zf1_ca = top.select(f"name CA and protein and resSeq >= {ZF1_RANGE_PDB[0]} and resSeq <= {ZF1_RANGE_PDB[1]}")
        zf2_ca = top.select(f"name CA and protein and resSeq >= {ZF2_RANGE_PDB[0]} and resSeq <= {ZF2_RANGE_PDB[1]}")
        print(f"  CA: {len(ca_all)}, RNA P: {len(rna_p)}, ZF1 CA: {len(zf1_ca)}, ZF2 CA: {len(zf2_ca)}")

        # RMSD series
        rmsd_ca = compute_rmsd_series(traj, ca_all, ca_all)
        rmsd_zf1 = compute_rmsd_series(traj, ca_all, zf1_ca) if len(zf1_ca) > 0 else np.full(traj.n_frames, np.nan)
        rmsd_zf2 = compute_rmsd_series(traj, ca_all, zf2_ca) if len(zf2_ca) > 0 else np.full(traj.n_frames, np.nan)
        rmsd_rna = compute_rmsd_series(traj, ca_all, rna_p) if len(rna_p) > 0 else np.full(traj.n_frames, np.nan)

        # Interface contacts (heavy-atom < 4 Å), subsampled
        prot_heavy = top.select("protein and (not name H and not type H)")
        rna_heavy = top.select("(resname A or resname C or resname G or resname U or "
                               "resname A5 or resname A3 or resname C5 or resname C3 or "
                               "resname G5 or resname G3 or resname U5 or resname U3) and (not name H and not type H)")
        if len(rna_heavy) > 0:
            stride = max(1, traj.n_frames // 200)
            cf_idx, cf_count = interface_contacts(traj, prot_heavy, rna_heavy, cutoff_nm=0.4, stride=stride)
            contact_time = time_ns[cf_idx]
        else:
            cf_idx = np.arange(traj.n_frames)
            cf_count = np.full(traj.n_frames, np.nan)
            contact_time = time_ns

        # Per-variant CSV
        df = pd.DataFrame({
            "time_ns": time_ns,
            "rmsd_ca_A": rmsd_ca,
            "rmsd_zf1_A": rmsd_zf1,
            "rmsd_zf2_A": rmsd_zf2,
            "rmsd_rnaP_A": rmsd_rna,
        })
        df.to_csv(OUT_DIR / f"{variant}_rmsd.csv", index=False)
        cdf = pd.DataFrame({"time_ns": contact_time, "n_pr_contacts_4A": cf_count})
        cdf.to_csv(OUT_DIR / f"{variant}_contacts.csv", index=False)

        # Summary using last 25%
        last_q = df.iloc[int(0.75 * len(df)):]
        last_q_c = cdf.iloc[int(0.75 * len(cdf)):]
        summary_rows.append({
            "variant": variant,
            "class": CLASS_MAP.get(variant, "?"),
            "n_frames": traj.n_frames,
            "t_final_ns": float(time_ns[-1]),
            "rmsd_ca_final_A": float(rmsd_ca[-1]),
            "rmsd_ca_last25_mean_A": float(last_q["rmsd_ca_A"].mean()),
            "rmsd_zf1_last25_mean_A": float(last_q["rmsd_zf1_A"].mean()),
            "rmsd_zf2_last25_mean_A": float(last_q["rmsd_zf2_A"].mean()),
            "rmsd_rnaP_final_A": float(rmsd_rna[-1]) if not np.isnan(rmsd_rna[-1]) else np.nan,
            "rmsd_rnaP_last25_mean_A": float(last_q["rmsd_rnaP_A"].mean()),
            "n_pr_contacts_last25_mean": float(last_q_c["n_pr_contacts_4A"].mean()) if "n_pr_contacts_4A" in last_q_c else np.nan,
        })

    summary = pd.DataFrame(summary_rows)
    # Compute delta vs WT
    if "WT" in summary["variant"].values:
        wt = summary[summary["variant"] == "WT"].iloc[0]
        for col in ["rmsd_ca_last25_mean_A", "rmsd_zf1_last25_mean_A",
                    "rmsd_zf2_last25_mean_A", "rmsd_rnaP_last25_mean_A",
                    "n_pr_contacts_last25_mean"]:
            summary[f"delta_{col}_vs_WT"] = summary[col] - wt[col]

    summary = summary.sort_values(["class", "variant"])
    out_csv = OUT_DIR / "rmsd_summary.csv"
    summary.to_csv(out_csv, index=False)
    cols = list(summary.columns)
    with (OUT_DIR / "rmsd_summary.md").open("w") as f:
        f.write("| " + " | ".join(cols) + " |\n")
        f.write("| " + " | ".join("---" for _ in cols) + " |\n")
        for _, row in summary.iterrows():
            f.write("| " + " | ".join(
                (f"{v:.4f}" if isinstance(v, float) and not np.isnan(v) else str(v))
                for v in row
            ) + " |\n")

    print("\n=== Summary (per-variant trajectory metrics) ===")
    with pd.option_context("display.max_columns", None, "display.width", 240):
        print(summary.to_string(index=False))
    print(f"\nWrote: {out_csv}")


if __name__ == "__main__":
    main()
