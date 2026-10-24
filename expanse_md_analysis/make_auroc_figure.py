#!/usr/bin/env python3
"""
Build a single AUROC summary figure that compiles every feature side-by-side.

  - Horizontal bar chart, sorted by AUROC, color-coded by feature origin
    (MD-Hbond, MD-bulk, AF3).
  - Reference lines at AUROC = 0.5 (random) and AUROC = the best AF3 feature.
  - Saves: plots/fig_feature_aurocs.png
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["savefig.dpi"] = 200
mpl.rcParams["figure.dpi"] = 110
mpl.rcParams["font.size"] = 10

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "expanse_md_analysis"
PLOTS = OUT / "plots"

df = pd.read_csv(OUT / "feature_aurocs.csv")
# Normalise origin and add a colour
origin_colour = {
    "MD-Hbond": "#c0392b",
    "MD": "#e67e22",
    "AF3": "#2980b9",
}
df["colour"] = df["origin"].map(origin_colour).fillna("#7f8c8d")

# Pretty feature names
PRETTY = {
    "hbond_RNA_F156": "MD H-bond rate at F156 ↔ RNA",
    "hbond_RNA_R160": "MD H-bond rate at R160 ↔ RNA",
    "hbond_RNA_H170": "MD H-bond rate at H170 ↔ RNA",
    "hbond_RNA_R184": "MD H-bond rate at R184 ↔ RNA",
    "hbond_RNA_Y154": "MD H-bond rate at Y154 ↔ RNA",
    "hbond_RNA_S185": "MD H-bond rate at S185 ↔ RNA",
    "hbond_RNA_T187": "MD H-bond rate at T187 ↔ RNA",
    "hbond_RNA_G204": "MD H-bond rate at G204 ↔ RNA",
    "hbond_RNA_C206": "MD H-bond rate at C206 ↔ RNA",
    "hbond_RNA_C208": "MD H-bond rate at C208 ↔ RNA",
    "hbond_RNA_C212": "MD H-bond rate at C212 ↔ RNA",
    "hbond_RNA_D210": "MD H-bond rate at D210 ↔ RNA",
    "hbond_RNA_R171": "MD H-bond rate at R171 ↔ RNA",
    "hbond_RNA_R200": "MD H-bond rate at R200 ↔ RNA",
    "rmsd_rnaP_last20pct_mean": "MD RNA-P RMSD (whole RNA)",
    "rmsd_ca_last20pct_mean":   "MD Cα RMSD (whole protein)",
    "rmsd_zf1_last20pct_mean":  "MD Cα RMSD (ZF1, res 152–179)",
    "rmsd_zf2_last20pct_mean":  "MD Cα RMSD (ZF2, res 190–217)",
    "contacts_last20pct_mean":  "MD protein–RNA heavy-atom contacts (4 Å)",
    "pair_ipTM_protein_rna":    "AF3 pair-ipTM (protein↔RNA)",
    "pair_pae_min_cross_mean":  "AF3 PAE_min (cross-chain)",
    "pLDDT":                    "AF3 pLDDT (whole)",
    "pLDDT_ZF1":                "AF3 pLDDT (ZF1)",
    "pLDDT_ZF2":                "AF3 pLDDT (ZF2)",
    "pLDDT_min":                "AF3 pLDDT_min",
    "RMSD_whole":               "AF3 whole-complex RMSD",
}
df["label"] = df["feature"].map(PRETTY).fillna(df["feature"])

# Best AF3 feature for the dashed reference line
af3 = df[df["origin"] == "AF3"]
af3_best = af3["AUROC"].max() if len(af3) else 0.5

# Sort ascending so the largest bar is at the top of a horizontal bar chart
df = df.sort_values("AUROC", ascending=True).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(11, 0.32 * len(df) + 1.5),
                       constrained_layout=True)
ax.barh(df["label"], df["AUROC"], color=df["colour"], edgecolor="black", lw=0.4)

# Annotate
for i, (lbl, auc, col) in enumerate(zip(df["label"], df["AUROC"], df["colour"])):
    ax.text(auc + 0.005, i, f"{auc:.3f}", va="center", fontsize=8)

ax.axvline(0.5, color="black", lw=0.6, ls="--", alpha=0.7)
ax.text(0.5, len(df) + 0.4, " random (0.5)", fontsize=8, va="center")
ax.axvline(af3_best, color=origin_colour["AF3"], lw=0.8, ls=":")
ax.text(af3_best, -0.6, f" best AF3 feature = {af3_best:.3f}",
        color=origin_colour["AF3"], fontsize=8, va="center")

ax.set_xlim(0.0, 1.02)
ax.set_xlabel("AUROC for deleterious vs neutral (6 vs 3 mutants, common 139-ns window)")
ax.set_title("ZFP36L2 + nonamer: per-feature ability to discriminate "
             "deleterious from neutral nsSNPs", fontsize=11)
# Legend
import matplotlib.patches as mpatches
handles = [mpatches.Patch(color=v, label=k) for k, v in origin_colour.items()]
ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)

out = PLOTS / "fig_feature_aurocs.png"
fig.savefig(out)
plt.close(fig)
print(f"saved {out}")
