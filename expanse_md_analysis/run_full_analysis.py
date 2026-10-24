#!/usr/bin/env python3
"""
End-to-end refreshed analysis runner for the ZFP36L2 + nonamer MD project.

Runs:
  1. parse_logs.py          → md_run_summary.csv  (T, E, density sanity per variant)
  2. compute_rmsd.py        → rmsd_summary.csv    (full-trajectory RMSD/contacts)
  3. plot_drift_and_hbonds.py → plots + 50-ns common-window tables
  4. AF3 vs MD AUROC comparison → feature_aurocs.csv

We invoke (1)-(3) as standalone modules so they remain re-runnable.
Step (4) is implemented inline below — it needs the outputs of step (3).

Also writes plots/fig_rmsd_full_trajectory.png with each variant's full
trajectory (regardless of length), and a per-variant trajectory-length banner.
"""
from __future__ import annotations
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "expanse_md_analysis"
PLOTS = OUT / "plots"
PLOTS.mkdir(exist_ok=True, parents=True)

VARIANTS = ["WT", "C206F", "C212A", "D210E", "G204D",
            "R160W", "R184C", "S185N", "T187S", "Y154H"]
CLASS = {
    "WT": "WT", "C206F": "deleterious", "C212A": "zn_knockout_control",
    "G204D": "deleterious", "R160W": "deleterious", "R184C": "deleterious",
    "Y154H": "deleterious",
    "D210E": "neutral", "S185N": "neutral", "T187S": "neutral",
}
TRUTH_DEL = {"Y154H", "R160W", "R184C", "G204D", "C206F", "C212A"}
TRUTH_NEU = {"S185N", "T187S", "D210E"}
COLOR = {"WT": "#0c4d8a", "deleterious": "#c0392b",
         "neutral": "#27ae60", "zn_knockout_control": "#7f8c8d"}


def _run(cmd: str):
    print(f"\n=== RUNNING: {cmd}")
    r = subprocess.run(cmd, shell=True, cwd=str(ROOT))
    if r.returncode != 0:
        sys.exit(f"step failed: {cmd}")


def step_1_2_3():
    _run("/Users/eric/miniconda3/bin/python3 expanse_md_analysis/parse_logs.py")
    _run("/Users/eric/miniconda3/bin/python3 expanse_md_analysis/compute_rmsd.py")
    _run("/Users/eric/miniconda3/bin/python3 expanse_md_analysis/plot_drift_and_hbonds.py")


def auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    n_concord = sum(1 for p in pos for n in neg if p > n)
    n_tie = sum(1 for p in pos for n in neg if p == n)
    return (n_concord + 0.5 * n_tie) / (len(pos) * len(neg))


def step_4_aurocs():
    print("\n=== STEP 4: AF3 vs MD AUROC catalog ===")
    # Use the most recent common-window summary (highest ns); fall back to legacy alias
    summary_candidates = sorted(OUT.glob("common*_summary.csv"),
                                key=lambda p: int("".join(c for c in p.stem if c.isdigit()) or 0),
                                reverse=True)
    summary_path = summary_candidates[0] if summary_candidates else OUT / "common50_summary.csv"
    md = pd.read_csv(summary_path)
    af3 = pd.read_csv(ROOT / "results/variant_avg_all.csv")
    af3 = af3.rename(columns={"Protein": "variant"})
    # The AF3 csv still uses the legacy "C208A" label; remap to "C212A" so it
    # joins onto the renamed MD folder.
    af3.loc[af3["variant"] == "C208A", "variant"] = "C212A"
    combined = md.merge(af3, on="variant", how="left")
    combined.to_csv(OUT / "af3_vs_md_combined.csv", index=False)

    mut = combined[combined["variant"] != "WT"].copy()
    y = mut["variant"].apply(lambda v: 1 if v in TRUTH_DEL else 0).values

    # bulk features
    feat = [
        ("rmsd_ca_last20pct_mean",  +1, "MD"),
        ("rmsd_zf1_last20pct_mean", +1, "MD"),
        ("rmsd_zf2_last20pct_mean", +1, "MD"),
        ("rmsd_rnaP_last20pct_mean", +1, "MD"),
        ("contacts_last20pct_mean", -1, "MD"),
        ("pair_ipTM_protein_rna",   -1, "AF3"),
        ("pair_pae_min_cross_mean", +1, "AF3"),
        ("pLDDT",                   -1, "AF3"),
        ("pLDDT_ZF1",               -1, "AF3"),
        ("pLDDT_ZF2",               -1, "AF3"),
        ("pLDDT_min",               -1, "AF3"),
        ("RMSD_whole",              +1, "AF3"),
    ]

    # per-residue H-bonds
    hb = pd.read_csv(OUT / "hbonds_per_residue_common50.csv", index_col=0)
    ZFP_TO_PDB = {154: 4, 156: 6, 160: 10, 170: 20, 171: 21,
                  184: 34, 185: 35, 187: 37, 200: 50, 204: 54,
                  206: 56, 208: 58, 210: 60, 212: 62}
    name = {154: "Y154", 156: "F156", 160: "R160", 170: "H170", 171: "R171",
            184: "R184", 185: "S185", 187: "T187", 200: "R200", 204: "G204",
            206: "C206", 208: "C208", 210: "D210", 212: "C212"}
    for z, p in ZFP_TO_PDB.items():
        if p in hb.index:
            col = f"hbond_RNA_{name[z]}"
            mut[col] = hb.loc[p, mut["variant"].values].values
            feat.append((col, -1, "MD-Hbond"))

    rows = []
    for f, sign, origin in feat:
        if f not in mut.columns:
            continue
        v = mut[f].values.astype(float) * sign
        if not np.isfinite(v).all():
            continue
        rows.append({"feature": f, "AUROC": round(auroc(v, y), 3),
                     "origin": origin})

    df = pd.DataFrame(rows).sort_values("AUROC", ascending=False)
    df.to_csv(OUT / "feature_aurocs.csv", index=False)
    print(df.to_string(index=False))


def step_5_full_traj_plot():
    """Time-series of the FULL trajectory (whatever ns each variant has)."""
    print("\n=== STEP 5: full-trajectory time series plot ===")
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)
    metrics = [("rmsd_ca_A",   "Cα RMSD"),
               ("rmsd_zf1_A",  "ZF1 Cα RMSD"),
               ("rmsd_zf2_A",  "ZF2 Cα RMSD"),
               ("rmsd_rnaP_A", "RNA P-atom RMSD")]
    for (m_key, m_title), ax in zip(metrics, axes.flat):
        for v in VARIANTS:
            csv = OUT / f"{v}_rmsd.csv"
            if not csv.exists():
                continue
            df = pd.read_csv(csv)
            ax.plot(df["time_ns"], df[m_key], lw=0.8,
                    color=COLOR[CLASS[v]],
                    alpha=0.7 if v != "WT" else 1.0, label=v)
        ax.set_xlabel("time (ns)")
        ax.set_ylabel("RMSD (Å)")
        ax.set_title(f"{m_title} (full available trajectory)")
        ax.grid(alpha=0.3)
        ax.legend(ncol=2, fontsize=7, frameon=False, loc="best")
    fig.suptitle("ZFP36L2 + nonamer: full available trajectory RMSDs",
                 fontsize=12)
    out = PLOTS / "fig_rmsd_full_trajectory.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"saved {out}")


def step_6_length_banner():
    """Bar chart showing how many ns each variant currently has."""
    print("\n=== STEP 6: trajectory length banner ===")
    fig, ax = plt.subplots(figsize=(11, 3.5), constrained_layout=True)
    rows = []
    for v in VARIANTS:
        prog = ROOT / "expanse_md_results" / v / "progress.txt"
        ns = float("nan")
        status = "?"
        if prog.exists():
            info = dict(line.split("=", 1)
                        for line in prog.read_text().splitlines()
                        if "=" in line)
            ns = float(info.get("ns_completed", "nan"))
            status = info.get("status", "?")
        rows.append({"variant": v, "ns": ns, "status": status,
                     "class": CLASS[v]})
    df = pd.DataFrame(rows)
    df = df.sort_values(["class", "variant"],
                        key=lambda c: c if c.name != "class" else c.map(
                            {"WT": 0, "deleterious": 1,
                             "zn_knockout_control": 2, "neutral": 3}))
    colors = [COLOR[c] for c in df["class"]]
    bars = ax.bar(df["variant"], df["ns"], color=colors, edgecolor="black",
                  linewidth=0.5)
    for b, st, ns in zip(bars, df["status"], df["ns"]):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 6,
                f"{ns:.0f} ns\n({st})", ha="center", fontsize=8)
    ax.set_ylabel("ns of production")
    ax.set_title("Production-trajectory length per variant (target = 400 ns)")
    ax.axhline(400, color="black", lw=0.6, ls="--", alpha=0.6)
    ax.axhline(50, color="grey", lw=0.6, ls="--", alpha=0.6)
    ax.set_ylim(0, 480)
    out = PLOTS / "fig_trajectory_length_status.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"saved {out}")


if __name__ == "__main__":
    step_1_2_3()
    step_4_aurocs()
    step_5_full_traj_plot()
    step_6_length_banner()
    print("\nAll steps complete.")
