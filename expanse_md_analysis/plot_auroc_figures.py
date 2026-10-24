#!/usr/bin/env python3
"""
Publication-quality AUROC figures for the ZFP36L2 nsSNP classification paper.

Generates 3 figures:
  1. AUROC bar chart — top 10 features (excluding C212A), with bootstrap error bars
  2. Per-variant metric ranking for F156 H-bond (the winner) — shows the del/neu separation
  3. ROC curves for the top 3 features (F156, AF3 ipTM, salt bridges)
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent / "plots"
OUT.mkdir(exist_ok=True)

# Load the replicate-averaged data
METRICS = pd.read_csv(ROOT / "expanse_md_analysis" / "multi_replicate_all_metrics.csv")
AUROC_CSV = ROOT / "expanse_md_analysis" / "auroc_no_c212a.csv"

TRUTH_DEL_5 = {"Y154H", "R160W", "R184C", "G204D", "C206F"}
CLASS = {"WT": "WT", "C206F": "del", "C212A": "del", "G204D": "del", "R160W": "del",
         "R184C": "del", "Y154H": "del", "D210E": "neu", "S185N": "neu", "T187S": "neu"}
CLASS_COLOR = {"WT": "#0c4d8a", "del": "#c0392b", "neu": "#27ae60"}
ORIGIN_COLOR = {"MD-Hbond": "#e74c3c", "MD-energy": "#e67e22", "MD-RMSD": "#3498db",
                "MD-contacts": "#9b59b6", "AF3": "#2ecc71", "MD-Qscore": "#95a5a6"}


def auroc(scores, labels):
    pos = scores[labels == 1]; neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0: return float('nan')
    nc = sum(1 for p in pos for n in neg if p > n)
    nt = sum(1 for p in pos for n in neg if p == n)
    return (nc + 0.5 * nt) / (len(pos) * len(neg))


# ============================================================================
# FIGURE 1: AUROC bar chart (top features, excluding C212A)
# ============================================================================
def figure1_auroc_bars():
    print("Figure 1: AUROC bar chart...")
    auroc_df = pd.read_csv(AUROC_CSV)

    # Select top features to show (curated for readability)
    top_features = [
        ("hbond_RNA_F156", "F156 H-bond", "MD-Hbond"),
        ("pair_ipTM_protein_rna", "AF3 pair-ipTM", "AF3"),
        ("salt_bridges_prot_rna", "Salt bridges", "MD-energy"),
        ("coulomb_prot_rna", "Coulomb elec.", "MD-energy"),
        ("pair_pae_min_cross_mean", "AF3 PAE", "AF3"),
        ("pLDDT_min", "AF3 pLDDT min", "AF3"),
        ("hbond_RNA_Y154", "Y154 H-bond", "MD-Hbond"),
        ("hbond_RNA_R160", "R160 H-bond", "MD-Hbond"),
        ("rmsd_ca_last20pct", "Cα RMSD", "MD-RMSD"),
        ("contacts_last20pct", "PR contacts", "MD-contacts"),
    ]

    # Get bootstrap CIs from the averaged AUROC file
    avg_auroc = pd.read_csv(ROOT / "expanse_md_analysis" / "multi_replicate_averaged_auroc.csv")

    fig, ax = plt.subplots(figsize=(10, 6))
    y_pos = np.arange(len(top_features))

    bars_data = []
    for feat_id, label, origin in top_features:
        row = auroc_df[auroc_df["feature"] == feat_id]
        if len(row) > 0:
            a = row.iloc[0]["AUROC"]
        else:
            a = 0.5

        # Get bootstrap CI if available
        ci_row = avg_auroc[avg_auroc["feature"] == feat_id]
        ci_lo = ci_row["AUROC_95CI_lo"].values[0] if len(ci_row) > 0 and not pd.isna(ci_row["AUROC_95CI_lo"].values[0]) else a
        ci_hi = ci_row["AUROC_95CI_hi"].values[0] if len(ci_row) > 0 and not pd.isna(ci_row["AUROC_95CI_hi"].values[0]) else a
        err_lo = max(0, a - ci_lo)
        err_hi = max(0, ci_hi - a)

        bars_data.append((label, a, err_lo, err_hi, ORIGIN_COLOR.get(origin, "#888")))

    labels = [b[0] for b in bars_data]
    vals = [b[1] for b in bars_data]
    errs = [[b[2] for b in bars_data], [b[3] for b in bars_data]]
    colors = [b[4] for b in bars_data]

    bars = ax.barh(y_pos, vals, color=colors, edgecolor='black',
                   linewidth=0.5, height=0.7)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('AUROC (5 del vs 3 neu, excluding C212A)', fontsize=11)
    ax.set_title('Feature Discrimination: Replicate-Averaged AUROC\n(5 deleterious vs 3 neutral, excluding C212A)', fontsize=13)
    ax.axvline(0.5, color='gray', ls='--', lw=0.8, alpha=0.5, label='Random (0.5)')
    ax.set_xlim(0, 1.05)
    ax.invert_yaxis()
    ax.grid(axis='x', alpha=0.2)

    # Value labels (positioned right after the bar end, no CI offset)
    for i, (lbl, v, elo, ehi, c) in enumerate(bars_data):
        ax.text(v + 0.015, i, f'{v:.3f}', va='center', fontsize=9, fontweight='bold')

    # Legend
    legend_elements = [
        Patch(facecolor=ORIGIN_COLOR["MD-Hbond"], label='MD H-bond'),
        Patch(facecolor=ORIGIN_COLOR["MD-energy"], label='MD Energy'),
        Patch(facecolor=ORIGIN_COLOR["AF3"], label='AF3 Static'),
        Patch(facecolor=ORIGIN_COLOR["MD-RMSD"], label='MD RMSD'),
        Patch(facecolor=ORIGIN_COLOR["MD-contacts"], label='MD Contacts'),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc='lower right')

    fig.tight_layout()
    out = OUT / "fig_auroc_bar_chart.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


# ============================================================================
# FIGURE 2: Per-variant F156 H-bond ranking (waterfall/lollipop)
# ============================================================================
def figure2_f156_ranking():
    print("Figure 2: F156 H-bond variant ranking...")

    # Replicate-average F156 per variant (exclude C212A)
    df = METRICS[METRICS["variant"] != "C212A"].copy()
    avg = df.groupby("variant")["hbond_RNA_F156"].agg(["mean", "std", "count"]).reset_index()
    avg.columns = ["variant", "f156_mean", "f156_std", "n_reps"]
    avg["f156_sem"] = avg["f156_std"] / np.sqrt(avg["n_reps"])
    avg["class"] = avg["variant"].map(CLASS)
    avg = avg.sort_values("f156_mean")

    fig, ax = plt.subplots(figsize=(10, 5))

    for i, (_, r) in enumerate(avg.iterrows()):
        color = CLASS_COLOR[r["class"]]
        ax.barh(i, r["f156_mean"], xerr=r["f156_sem"], color=color,
                edgecolor='black', linewidth=0.5, height=0.6, capsize=4,
                error_kw={'lw': 1.5})
        ax.text(r["f156_mean"] + r["f156_sem"] + 0.01, i,
                f'{r["f156_mean"]:.3f} ± {r["f156_sem"]:.3f} ({int(r["n_reps"])}r)',
                va='center', fontsize=9)

    ax.set_yticks(range(len(avg)))
    ax.set_yticklabels([f'{r["variant"]}' for _, r in avg.iterrows()], fontsize=11, fontweight='bold')
    ax.set_xlabel('F156–RNA H-bond frequency (per frame)', fontsize=11)
    ax.set_title('Per-Variant F156 H-bond: Deleterious < Neutral\n(replicate-averaged ± SEM, excluding C212A)', fontsize=13)
    ax.grid(axis='x', alpha=0.2)
    ax.invert_yaxis()

    # Add a shaded region for the del/neu boundary
    ax.axvspan(0.505, 0.523, alpha=0.15, color='gray', label='Del/Neu overlap zone')

    legend_elements = [
        Patch(facecolor=CLASS_COLOR["del"], label='Deleterious (5)'),
        Patch(facecolor=CLASS_COLOR["neu"], label='Neutral (3)'),
        Patch(facecolor=CLASS_COLOR["WT"], label='Wild-type'),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc='lower right')

    fig.tight_layout()
    out = OUT / "fig_f156_variant_ranking.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


# ============================================================================
# FIGURE 3: Per-variant comparison of top features (grouped bar)
# ============================================================================
def figure3_feature_comparison():
    print("Figure 3: Feature comparison across variants...")

    df = METRICS[METRICS["variant"] != "C212A"].copy()

    # Load AF3 data
    af3 = pd.read_csv(ROOT / "results" / "variant_avg_all.csv").rename(columns={"Protein": "variant"})
    af3.loc[af3["variant"] == "C208A", "variant"] = "C212A"
    af3 = af3[af3["variant"] != "C212A"]

    # Average MD metrics
    md_avg = df.groupby("variant")[["hbond_RNA_F156", "salt_bridges_prot_rna", "coulomb_prot_rna"]].mean().reset_index()
    md_avg = md_avg.merge(af3[["variant", "pair_ipTM_protein_rna"]], on="variant", how="left")
    md_avg["class"] = md_avg["variant"].map(CLASS)

    # Normalize each metric to [0, 1] for visual comparison
    # For H-bond and ipTM: higher = more WT-like = more neutral
    # For Coulomb: more negative = stronger binding = more neutral → negate + normalize
    features = {
        "hbond_RNA_F156": ("F156 H-bond", False),
        "pair_ipTM_protein_rna": ("AF3 pair-ipTM", False),
        "salt_bridges_prot_rna": ("Salt bridges", False),
    }

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

    # Sort by F156 (the best metric)
    order = md_avg.sort_values("hbond_RNA_F156")["variant"].tolist()

    for idx, (feat, (label, invert)) in enumerate(features.items()):
        ax = axes[idx]
        vals = md_avg.set_index("variant").loc[order, feat].values
        classes = md_avg.set_index("variant").loc[order, "class"].values
        colors = [CLASS_COLOR[c] for c in classes]

        y_pos = np.arange(len(order))
        ax.barh(y_pos, vals, color=colors, edgecolor='black', linewidth=0.5, height=0.6)

        for i, (v, val) in enumerate(zip(order, vals)):
            ax.text(val + (max(vals) - min(vals)) * 0.02, i, f'{val:.3f}', va='center', fontsize=8)

        ax.set_yticks(y_pos)
        if idx == 0:
            ax.set_yticklabels(order, fontsize=10, fontweight='bold')
        ax.set_xlabel(label, fontsize=10)
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.grid(axis='x', alpha=0.2)
        ax.invert_yaxis()

    # Legend on first panel
    legend_elements = [
        Patch(facecolor=CLASS_COLOR["del"], label='Deleterious'),
        Patch(facecolor=CLASS_COLOR["neu"], label='Neutral'),
        Patch(facecolor=CLASS_COLOR["WT"], label='Wild-type'),
    ]
    axes[0].legend(handles=legend_elements, fontsize=8, loc='lower right')

    fig.suptitle('Per-Variant Comparison of Top Discriminative Features\n(replicate-averaged, excluding C212A)', fontsize=13)
    fig.tight_layout()
    out = OUT / "fig_feature_comparison.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    figure1_auroc_bars()
    figure2_f156_ranking()
    figure3_feature_comparison()
    print("\nAll AUROC figures generated.")
