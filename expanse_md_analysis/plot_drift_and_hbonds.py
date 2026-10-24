#!/usr/bin/env python3
"""
Generate two publication-style figures for the ZFP36L2 + nonamer MD study:

  Figure 1: MD-drift deltas vs WT (common 0–50 ns window)
            - Cα RMSD, ZF1 RMSD, ZF2 RMSD, RNA-P RMSD, protein–RNA contact count
            - For every variant, Δ = variant − WT, restricted to first 50 ns
            - Grouped by class (deleterious / neutral)

  Figure 2: Per-residue hydrogen-bond counts to RNA (Δ vs WT)
            - For every ZFP36L2 residue, the average number of H-bonds with any
              RNA atom in the common 0–50 ns window
            - Heatmap of (residues × variants) and a bar chart for the
              key ARE-binding residues that Ramos et al. (2024) highlighted:
              Y154, R160, R184, S185, T187, G204, C206, C208, D210, H212

Compares qualitatively to Ramos et al. 2024 (RNA Biology 22:1):
  - R184C "abolished H-bonds with surrounding residues, 6 → 3"  (their Fig S13C)
  - Y154H / R160W "higher probability of forming 1–2 extra H-bonds"
  - G204D / C206F "reduced probability of single H-bonds"
We compute the protein-RNA H-bonds (instead of protein-protein) to directly
test the RNA-binding hypothesis on our AF3 starting structures.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import mdtraj as md
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["font.size"] = 9
mpl.rcParams["axes.titlesize"] = 10
mpl.rcParams["axes.labelsize"] = 9
mpl.rcParams["legend.fontsize"] = 8
mpl.rcParams["savefig.dpi"] = 200
mpl.rcParams["figure.dpi"] = 110

ROOT = Path(__file__).resolve().parent.parent
EXPANSE_DIR = ROOT / "expanse_md_results"
OUT_DIR = ROOT / "expanse_md_analysis"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PLOTS = OUT_DIR / "plots"
PLOTS.mkdir(exist_ok=True)

VARIANTS = ["WT", "C206F", "C212A", "D210E", "G204D",
            "R160W", "R184C", "S185N", "T187S", "Y154H"]

# Class assignment (Y154H labeled deleterious in *this* dataset per Block_Design.txt;
# the Ramos paper labels it deleterious for binding too).
# C212A was previously mislabeled C208A; mutation site is C212 (a Zn-coord Cys of ZF2)
# D210E run actually contains D219E (folder name kept for continuity)
CLASS = {
    "WT": "WT",
    "C206F": "deleterious",
    "C212A": "zn_knockout_control",
    "G204D": "deleterious",
    "R160W": "deleterious",
    "R184C": "deleterious",
    "Y154H": "deleterious",
    "D210E": "neutral",
    "S185N": "neutral",
    "T187S": "neutral",
}
COLOR = {"WT": "#0c4d8a", "deleterious": "#c0392b",
         "neutral": "#27ae60", "zn_knockout_control": "#7f8c8d"}

# Ramos-paper key ARE-binding residues (ZFP36L2 numbering) → PDB resSeq (–150)
KEY_RES_ZFP = [154, 156, 160, 170, 171, 184, 185, 187, 200, 204, 206, 208, 210, 212]
ZFP_TO_PDB = {z: z - 150 for z in KEY_RES_ZFP}

COMMON_WINDOW_NS = 269.0  # all variants completed 400ns; DCD common floor is 269.7ns


def is_rna(res) -> bool:
    n = res.name.upper().strip()
    return n in {"A", "C", "G", "U", "A5", "A3", "C5", "C3",
                 "G5", "G3", "U5", "U3"}


def load_solute_traj(variant: str, max_ns: float = COMMON_WINDOW_NS):
    """Load only protein+RNA+Zn atoms, only the frames covering 0..max_ns."""
    var_dir = EXPANSE_DIR / variant
    top = var_dir / "equilibrated.pdb"
    dcd = var_dir / "production.dcd"
    if not (top.exists() and dcd.exists()):
        return None
    full_top = md.load_topology(str(top))
    keep = np.unique(np.concatenate([
        full_top.select("protein"),
        full_top.select(" or ".join(f"resname {r}" for r in
                        ("A","C","G","U","A5","A3","C5","C3","G5","G3","U5","U3"))),
        full_top.select("name ZN or element Zn"),
    ]))
    # Production.dcd was written every 100 ps -> 10 frames per ns.
    # Some DCDs are partial 4-ns segments; we use whatever we have but cap at max_ns.
    traj = md.load(str(dcd), top=str(top), atom_indices=keep)
    # 100 ps per frame:
    n_max = int(max_ns * 10) + 1
    if traj.n_frames > n_max:
        traj = traj[:n_max]
    return traj


def compute_rmsd_to_first(traj, ref_idx, target_idx) -> np.ndarray:
    traj = traj.superpose(traj, frame=0, atom_indices=ref_idx, ref_atom_indices=ref_idx)
    diff = traj.xyz[:, target_idx, :] - traj.xyz[0:1, target_idx, :]
    return np.sqrt((diff ** 2).sum(axis=2).mean(axis=1)) * 10.0  # Å


def compute_protein_rna_hbonds_per_residue(traj):
    """Use Baker-Hubbard H-bond detection on the union of protein + RNA atoms.
    For each ZFP36L2 protein residue, count the fraction of frames in which
    *any* H-bond connects that residue to *any* RNA atom.
    Returns: {pdb_resSeq: mean_count_per_frame}
    """
    top = traj.topology
    rna_atoms = {a.index for a in top.atoms if is_rna(a.residue)}
    prot_residues = [r for r in top.residues if r.is_protein]
    pdb_resids = [r.resSeq for r in prot_residues]
    n_res = len(prot_residues)
    res_idx_of_atom = {a.index: a.residue.resSeq
                       for r in prot_residues for a in r.atoms}

    counts = np.zeros((traj.n_frames, n_res), dtype=int)
    res_pos = {r.resSeq: i for i, r in enumerate(prot_residues)}

    # Run BH per frame (cheap because traj is small after stripping water)
    for fi in range(traj.n_frames):
        hb = md.baker_hubbard(traj[fi], freq=0.0, exclude_water=True)
        # hb rows: (donor_h_idx? actually donor-h-acceptor)
        # mdtraj.baker_hubbard returns array of (donor, hydrogen, acceptor) atom indices
        for donor, hydrogen, acceptor in hb:
            d_res = res_idx_of_atom.get(donor)
            a_res = res_idx_of_atom.get(acceptor)
            d_in_rna = donor in rna_atoms
            a_in_rna = acceptor in rna_atoms
            # We only want protein↔RNA H-bonds (any direction)
            if d_in_rna and a_res is not None:
                counts[fi, res_pos[a_res]] += 1
            elif a_in_rna and d_res is not None:
                counts[fi, res_pos[d_res]] += 1

    return pd.Series(counts.mean(axis=0), index=pdb_resids,
                     name="mean_hbonds_to_rna")


# ============================================================================
# 1) Compute per-variant metrics on the common 50 ns window
# ============================================================================

def compute_all():
    rows = []
    hbonds_by_variant = {}
    rmsd_series_by_variant = {}
    contacts_by_variant = {}

    for variant in VARIANTS:
        print(f"[{variant}] loading <={COMMON_WINDOW_NS} ns ...")
        traj = load_solute_traj(variant, COMMON_WINDOW_NS)
        if traj is None:
            print(f"  skip (missing files)")
            continue
        top = traj.topology
        ca = top.select("name CA and protein")
        rna_p = top.select(" or ".join(
            f"(resname {r} and name P)" for r in
            ("A","C","G","U","A5","A3","C5","C3","G5","G3","U5","U3")))
        zf1 = top.select("name CA and protein and resSeq >= 2 and resSeq <= 29")
        zf2 = top.select("name CA and protein and resSeq >= 40 and resSeq <= 67")
        prot_heavy = top.select("protein and not name H and not type H")
        rna_heavy = top.select(" or ".join(
            f"(resname {r} and not name H and not type H)" for r in
            ("A","C","G","U","A5","A3","C5","C3","G5","G3","U5","U3")))

        rmsd_ca  = compute_rmsd_to_first(traj, ca, ca)
        rmsd_zf1 = compute_rmsd_to_first(traj, ca, zf1)
        rmsd_zf2 = compute_rmsd_to_first(traj, ca, zf2)
        rmsd_rna = compute_rmsd_to_first(traj, ca, rna_p) if len(rna_p) else np.full(traj.n_frames, np.nan)

        # contacts (protein heavy ↔ RNA heavy within 4 Å) every ~100 ps
        pairs = np.array(np.meshgrid(prot_heavy, rna_heavy)).T.reshape(-1, 2)
        # Compute in chunks to limit memory
        n_pr = np.zeros(traj.n_frames, dtype=int)
        chunk = 100
        for start in range(0, traj.n_frames, chunk):
            sub = traj[start:start+chunk]
            d = md.compute_distances(sub, pairs, periodic=True)
            n_pr[start:start+sub.n_frames] = (d < 0.4).sum(axis=1)

        # Time axis (100 ps per frame)
        time_ns = np.arange(traj.n_frames) * 0.1
        rmsd_series_by_variant[variant] = {
            "time_ns": time_ns, "rmsd_ca": rmsd_ca,
            "rmsd_zf1": rmsd_zf1, "rmsd_zf2": rmsd_zf2,
            "rmsd_rnaP": rmsd_rna,
        }
        contacts_by_variant[variant] = (time_ns, n_pr)

        rows.append({
            "variant": variant,
            "class": CLASS[variant],
            "n_frames": traj.n_frames,
            "t_final_ns": float(time_ns[-1]),
            "rmsd_ca_mean": float(rmsd_ca.mean()),
            "rmsd_ca_last20pct_mean": float(rmsd_ca[int(0.8*len(rmsd_ca)):].mean()),
            "rmsd_zf1_last20pct_mean": float(rmsd_zf1[int(0.8*len(rmsd_zf1)):].mean()),
            "rmsd_zf2_last20pct_mean": float(rmsd_zf2[int(0.8*len(rmsd_zf2)):].mean()),
            "rmsd_rnaP_last20pct_mean": float(np.nanmean(rmsd_rna[int(0.8*len(rmsd_rna)):])),
            "contacts_last20pct_mean": float(n_pr[int(0.8*len(n_pr)):].mean()),
        })

        # H-bond per residue
        print(f"  computing protein↔RNA H-bonds per residue ...")
        hb = compute_protein_rna_hbonds_per_residue(traj)
        hbonds_by_variant[variant] = hb

    summary = pd.DataFrame(rows)
    win_tag = f"common{int(COMMON_WINDOW_NS)}"
    summary.to_csv(OUT_DIR / f"{win_tag}_summary.csv", index=False)

    # Build delta-vs-WT table
    if "WT" in summary["variant"].values:
        wt = summary[summary["variant"] == "WT"].iloc[0]
        for col in ["rmsd_ca_last20pct_mean", "rmsd_zf1_last20pct_mean",
                    "rmsd_zf2_last20pct_mean", "rmsd_rnaP_last20pct_mean",
                    "contacts_last20pct_mean"]:
            summary["d_" + col + "_vs_WT"] = summary[col] - wt[col]
        summary.to_csv(OUT_DIR / f"{win_tag}_summary.csv", index=False)
        # also write to legacy filename so downstream scripts keep working
        summary.to_csv(OUT_DIR / "common50_summary.csv", index=False)

    # Build per-residue H-bond table
    if "WT" in hbonds_by_variant:
        all_res = sorted(set().union(*[set(s.index) for s in hbonds_by_variant.values()]))
        hb_df = pd.DataFrame(index=all_res)
        for v, s in hbonds_by_variant.items():
            hb_df[v] = s.reindex(all_res)
        hb_df.to_csv(OUT_DIR / f"hbonds_per_residue_{win_tag}.csv")
        hb_df.to_csv(OUT_DIR / "hbonds_per_residue_common50.csv")  # legacy
        # delta vs WT
        hb_delta = hb_df.sub(hb_df["WT"], axis=0).drop(columns="WT")
        hb_delta.to_csv(OUT_DIR / "hbonds_per_residue_delta_vs_WT.csv")
    else:
        hb_df = pd.DataFrame()
        hb_delta = pd.DataFrame()

    return summary, hb_df, hb_delta, rmsd_series_by_variant, contacts_by_variant


# ============================================================================
# 2) Figure 1: drift deltas vs WT (5 panels)
# ============================================================================

def fig_drift(summary: pd.DataFrame, rmsd_series, contacts):
    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2])

    # --- Panel 1+2+3+4: time series (Cα, ZF1, ZF2, RNA-P RMSD) ---
    series_axes = [(0,0,"rmsd_ca","Cα RMSD"),
                   (0,1,"rmsd_zf1","ZF1 (152–179) Cα RMSD"),
                   (1,0,"rmsd_zf2","ZF2 (190–217) Cα RMSD"),
                   (1,1,"rmsd_rnaP","RNA P-atom RMSD")]
    for r, c, key, title in series_axes:
        ax = fig.add_subplot(gs[r, c])
        for v in VARIANTS:
            if v not in rmsd_series: continue
            d = rmsd_series[v]
            ax.plot(d["time_ns"], d[key], lw=1.0,
                    color=COLOR[CLASS[v]], alpha=0.6 if v != "WT" else 1.0,
                    label=v)
        ax.set_xlabel("time (ns)")
        ax.set_ylabel("RMSD (Å)")
        ax.set_title(title + f"   (0–{COMMON_WINDOW_NS:.0f} ns common window)")
        ax.grid(alpha=0.3)
        if r == 0 and c == 0:
            ax.legend(ncol=2, loc="upper left", fontsize=7, frameon=False)

    # --- Bottom panel: delta-vs-WT bar chart (5 metrics) ---
    ax = fig.add_subplot(gs[2, :])
    metric_cols = [
        ("d_rmsd_ca_last20pct_mean_vs_WT",   "ΔCα RMSD"),
        ("d_rmsd_zf1_last20pct_mean_vs_WT",  "ΔZF1 RMSD"),
        ("d_rmsd_zf2_last20pct_mean_vs_WT",  "ΔZF2 RMSD"),
        ("d_rmsd_rnaP_last20pct_mean_vs_WT", "ΔRNA-P RMSD"),
        ("d_contacts_last20pct_mean_vs_WT",  "Δprotein–RNA contacts (4 Å)"),
    ]
    # Normalize contacts axis by 50 for visibility
    df = summary[summary["variant"] != "WT"].sort_values(
        ["class", "variant"],
        key=lambda c: c if c.name != "class" else c.map({"deleterious":0,"neutral":1}))
    variants_plot = list(df["variant"])
    x = np.arange(len(variants_plot))
    width = 0.16
    for i, (col, label) in enumerate(metric_cols):
        vals = df[col].values
        if "contacts" in col:
            vals = vals / 50.0  # rescale so it shares a y-axis
            label = label + " /50"
        bars = ax.bar(x + (i-2)*width, vals, width=width, label=label,
                      edgecolor="black", linewidth=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{v}\n({CLASS[v]})" for v in variants_plot], rotation=0)
    ax.axhline(0, color="black", lw=0.7)
    ax.set_ylabel("Δ (variant – WT)\n(Å for RMSDs; contacts /50)")
    ax.set_title("MD drift deltas vs WT (mean of last 20 % of 50-ns common window)")
    ax.legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
    ax.grid(axis="y", alpha=0.3)

    fig.suptitle("ZFP36L2 + nonamer MD: per-variant drift deltas vs WT  "
                 f"(common {COMMON_WINDOW_NS:.0f}-ns window)", fontsize=12)
    out = PLOTS / f"fig_drift_deltas_vs_WT_common{int(COMMON_WINDOW_NS)}ns.png"
    fig.savefig(out)
    fig.savefig(PLOTS / "fig_drift_deltas_vs_WT_common50ns.png")  # legacy filename
    plt.close(fig)
    print(f"saved {out}")


# ============================================================================
# 3) Figure 2: per-residue H-bond counts (heatmap + key residue bar chart)
# ============================================================================

def fig_hbonds(hb_df: pd.DataFrame, hb_delta: pd.DataFrame):
    if hb_df.empty:
        print("no H-bond data; skipping fig_hbonds")
        return

    # Order variants: WT first, then deleterious, then neutral
    ordered = ["WT"] + [v for v in VARIANTS if v != "WT" and CLASS[v] == "deleterious"] + \
              [v for v in VARIANTS if v != "WT" and CLASS[v] == "neutral"]
    ordered = [v for v in ordered if v in hb_df.columns]
    hb_o = hb_df[ordered]

    # ----- (a) heatmap of mean H-bonds to RNA per residue -----
    fig = plt.figure(figsize=(14, 10), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[2, 1])

    ax_h = fig.add_subplot(gs[0, 0])
    matrix = hb_o.T.values  # rows=variant, cols=residue
    im = ax_h.imshow(matrix, aspect="auto", cmap="viridis",
                     interpolation="nearest")
    ax_h.set_yticks(np.arange(len(ordered)))
    ax_h.set_yticklabels([f"{v} ({CLASS[v]})" for v in ordered])
    # Residue labels: every 5 residues, in ZFP36L2 numbering (= pdb + 150)
    ax_h.set_xticks(np.arange(len(hb_o.index))[::5])
    ax_h.set_xticklabels([str(r + 150) for r in hb_o.index[::5]], rotation=0)
    ax_h.set_xlabel("ZFP36L2 residue (UniProt numbering)")
    ax_h.set_title("Mean protein–RNA hydrogen bonds per residue  "
                   f"(common {COMMON_WINDOW_NS:.0f}-ns window)")
    cbar = fig.colorbar(im, ax=ax_h, pad=0.01, shrink=0.8)
    cbar.set_label("mean H-bonds / frame")

    # annotate key ARE-binding residues
    for zres, pres in ZFP_TO_PDB.items():
        if pres in hb_o.index:
            x = list(hb_o.index).index(pres)
            ax_h.axvline(x, color="white", lw=0.4, alpha=0.6)
            ax_h.text(x, -0.7, str(zres), color="black", ha="center",
                      fontsize=7, rotation=90)

    # ----- (b) bar chart on key residues, ΔvsWT -----
    ax_b = fig.add_subplot(gs[1, 0])
    key_pdb = [ZFP_TO_PDB[z] for z in KEY_RES_ZFP if ZFP_TO_PDB[z] in hb_delta.index]
    key_zfp = [pres + 150 for pres in key_pdb]
    plot_variants = [v for v in ordered if v != "WT"]
    width = 0.85 / max(1, len(plot_variants))
    x = np.arange(len(key_pdb))
    for i, v in enumerate(plot_variants):
        vals = hb_delta.loc[key_pdb, v].values
        ax_b.bar(x + (i - (len(plot_variants)-1)/2) * width, vals,
                 width=width, label=v, color=COLOR[CLASS[v]],
                 alpha=0.95 if CLASS[v] == "deleterious" else 0.55,
                 edgecolor="black", linewidth=0.4)
    ax_b.axhline(0, color="black", lw=0.6)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([f"{aa3}{z}" for aa3, z in zip(
        ["Y","F","R","H","R","R","S","T","R","G","C","C","D","H"], key_zfp)])
    ax_b.set_ylabel("Δ mean H-bonds vs WT")
    ax_b.set_title("Per-residue H-bond change vs WT at canonical ARE-binding residues "
                   "(Ramos et al. 2024 highlights: Y154, R160, R184, G204, C206)")
    ax_b.legend(ncol=5, fontsize=7, frameon=False, loc="upper right")
    ax_b.grid(axis="y", alpha=0.3)

    out = PLOTS / f"fig_hbonds_per_residue_common{int(COMMON_WINDOW_NS)}ns.png"
    fig.savefig(out)
    fig.savefig(PLOTS / "fig_hbonds_per_residue_common50ns.png")  # legacy filename
    plt.close(fig)
    print(f"saved {out}")

    # Also write a compact summary table for the discussion
    key_table = hb_o.loc[[p for p in key_pdb if p in hb_o.index]].copy()
    key_table.index = [pdb + 150 for pdb in key_table.index]
    key_table.index.name = "ZFP36L2_resnum"
    key_table.to_csv(OUT_DIR / f"hbonds_key_residues_common{int(COMMON_WINDOW_NS)}.csv")
    key_table.to_csv(OUT_DIR / "hbonds_key_residues_common50.csv")  # legacy


# ============================================================================
# 4) Comparison with Ramos et al.
# ============================================================================

RAMOS_EXPECTATION = {
    # variant: "expected_sign_of_per-residue Hbond change vs WT  (paper context)"
    "Y154H": ("Y154", "small ↑ or ≈0 (paper reports ↑ intra-protein H-bonds)"),
    "R160W": ("R160", "loss of guanidinium → expected ↓ at R160 with RNA"),
    "R184C": ("R184", "from 6 → 3 H-bonds (paper Figure S13C) → strong ↓"),
    "G204D": ("G204", "small ↓ (paper says 'reduced probability')"),
    "C206F": ("C206", "small ↓ (paper says 'reduced probability')"),
    "C212A": ("C212", "Zn coordination loss in ZF2 (this run; mislabeled as C208A)"),
    "D210E": ("D210", "neutral – expected ≈0 (note: real mutation is D219E)"),
    "S185N": ("S185", "neutral – expected ≈0"),
    "T187S": ("T187", "neutral – expected ≈0"),
}


def comparison_table(hb_delta: pd.DataFrame):
    rows = []
    for v, (resname, expectation) in RAMOS_EXPECTATION.items():
        zres = int(resname[1:])
        pdb_res = zres - 150
        if v not in hb_delta.columns or pdb_res not in hb_delta.index:
            continue
        delta = float(hb_delta.loc[pdb_res, v])
        agree = "match" if (
            (delta < -0.05 and "↓" in expectation) or
            (delta > 0.05 and "↑" in expectation) or
            (-0.05 <= delta <= 0.05 and "≈0" in expectation)
        ) else "differs"
        rows.append({
            "variant": v, "key_residue": resname,
            "delta_Hbonds_vs_WT": round(delta, 3),
            "paper_expectation": expectation,
            "qualitative_agreement": agree,
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT_DIR / "ramos_comparison_table.csv", index=False)
    print("\n=== Ramos et al. (2024) qualitative comparison ===")
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print(df.to_string(index=False))


# ============================================================================
if __name__ == "__main__":
    summary, hb_df, hb_delta, rmsd_series, contacts = compute_all()
    fig_drift(summary, rmsd_series, contacts)
    fig_hbonds(hb_df, hb_delta)
    comparison_table(hb_delta)
    print("\nDone.  See: expanse_md_analysis/plots/")
