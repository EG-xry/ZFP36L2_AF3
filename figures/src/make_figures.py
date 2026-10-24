#!/usr/bin/env python3
"""
Build every data figure in the write-up from files committed to this repo.

    python figures/src/make_figures.py        (repo .venv: numpy, pandas, matplotlib)

Inputs
    figures/data/af3_per_model.csv               af3_metrics.py (from Data/)
    figures/data/t156_rna_hbond_series_1ns.csv   extract_t156_series.py (from DCDs)
    expanse_md_analysis/*.csv                    the project's MD analysis outputs
    analysis/interface_metrics_models.csv        non-binding RNA control (AF3)

Outputs: figures/<name>.svg (editable text) + figures/<name>.png, and
figures/data/tables.md with the numbers behind each figure.

Naming: the MD code calls resSeq 6 "F156", but that residue is Thr156 of the
RYKTEL motif. Five other per-residue columns carry the same kind of shifted
name; FEATURE_LABELS below uses the residue actually present at each position.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

import style as S

ROOT = Path(__file__).resolve().parents[2]
MDA = ROOT / "expanse_md_analysis"
FIGDATA = ROOT / "figures" / "data"

# MD analysis code names C208A "C212A" (the residue it mutates in construct numbering).
MD_NAME = {"C212A": "C208A"}
N_FEATURES_SCREENED = 32  # 30 in auroc_no_c212a.csv + Q-score + MM-GBSA (figure M2)

FEATURE_LABELS = {
    "hbond_RNA_F156": "T156 H-bonds to RNA",
    "hbond_RNA_Y154": "Y154 H-bonds", "hbond_RNA_R160": "R160 H-bonds",
    "hbond_RNA_H170": "Y170 H-bonds", "hbond_RNA_R171": "G171 H-bonds",
    "hbond_RNA_R184": "R184 H-bonds", "hbond_RNA_S185": "S185 H-bonds",
    "hbond_RNA_T187": "T187 H-bonds", "hbond_RNA_R200": "F200 H-bonds",
    "hbond_RNA_G204": "G204 H-bonds", "hbond_RNA_C206": "C206 H-bonds",
    "hbond_RNA_C208": "Y208 H-bonds", "hbond_RNA_D210": "P210 H-bonds",
    "hbond_RNA_C212": "C212 H-bonds",
    "pair_ipTM_protein_rna": "pair-ipTM (protein–RNA)",
    "pair_pae_min_cross_mean": "Cross-chain PAE",
    "pLDDT": "Mean pLDDT", "pLDDT_ZF1": "pLDDT, finger 1",
    "pLDDT_ZF2": "pLDDT, finger 2", "pLDDT_min": "Min pLDDT (residue E220)",
    "RMSD_whole": "RMSD to WT model",
    "rmsd_ca_last20pct": "Cα RMSD", "rmsd_zf1_last20pct": "Finger-1 RMSD",
    "rmsd_zf2_last20pct": "Finger-2 RMSD", "rmsd_rnaP_last20pct": "RNA phosphate RMSD",
    "contacts_last20pct": "Protein–RNA contacts",
    "coulomb_prot_rna": "Protein–RNA Coulomb energy",
    "salt_bridges_prot_rna": "Protein–RNA salt bridges",
    "interface_atom_count_5A": "Interface atoms (5 Å)",
    "binding_residue_rmsf_A": "Binding-residue RMSF",
    "qscore": "Native-contact Q-score", "mmgbsa": "MM-GBSA (Zn stripped)",
}

TABLES: list[str] = []


# ── Statistics ────────────────────────────────────────────────────────────
def auroc(values: dict, sign: int) -> float:
    """P(deleterious scores above neutral); score = sign * value, sign fixed a priori."""
    hits = 0.0
    for d in S.DELETERIOUS:
        for n in S.NEUTRAL:
            a, b = sign * values[d], sign * values[n]
            hits += 1.0 if a > b else 0.5 if a == b else 0.0
    return hits / (len(S.DELETERIOUS) * len(S.NEUTRAL))


def mann_whitney_null(m: int, n: int) -> np.ndarray:
    """Exact null distribution of U for groups of size m, n (probabilities, U = 0..mn)."""
    f = [[None] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 or j == 0:
                f[i][j] = np.ones(1)
                continue
            # The largest value belongs to group 1 (beats all j) or to group 2.
            a = np.zeros(i * j + 1)
            prev = f[i - 1][j]
            a[j:j + len(prev)] += prev
            left = f[i][j - 1]
            a[:len(left)] += left
            f[i][j] = a
    counts = f[m][n]
    return counts / counts.sum()


def tail_p(m: int, n: int) -> np.ndarray:
    """One-sided p = P(U >= u) for u = 0..mn."""
    return np.cumsum(mann_whitney_null(m, n)[::-1])[::-1]


def critical_auroc(m: int, n: int, alpha: float) -> float:
    p = tail_p(m, n)
    ok = np.nonzero(p <= alpha)[0]
    return ok[0] / (m * n) if len(ok) else np.nan


# ── Small drawing helpers ─────────────────────────────────────────────────
def dot(ax, x, y, cls, s=46, ring=True, z=3, alpha=1.0, **kw):
    ax.scatter(x, y, s=s, color=S.CLASS_COLOR[cls], marker=S.CLASS_MARKER[cls],
               edgecolors=S.SURFACE if ring else "none", linewidths=1.3 if ring else 0,
               zorder=z, alpha=alpha, **kw)


def legend(fig, items, top_in, x=0.012):
    """items: (color, marker, label). One row, left-aligned under the header."""
    handles = [Line2D([], [], ls="", marker=m, markersize=7, markerfacecolor=c,
                      markeredgecolor=S.SURFACE, label=l) for c, m, l in items]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(x, S.top(fig, top_in)),
               ncol=len(items), handletextpad=0.3, columnspacing=1.6, borderaxespad=0)


def class_items(classes):
    return [(S.CLASS_COLOR[c], S.CLASS_MARKER[c], S.CLASS_LABEL[c]) for c in classes]


METHOD_ITEMS = [(S.MD, "o", "MD feature"), (S.AF3, "o", "AlphaFold 3 confidence")]


def significance_band(ax, crit: float, label_y: float = 1.01) -> None:
    """Wash over AUROC values that would reach uncorrected p < 0.05; label above the axes."""
    lo = crit - 1 / 30  # halfway to the next attainable value below (step = 1/15)
    ax.axvspan(lo, 1.02, color=S.WASH, lw=0, zorder=0)
    ax.text(lo + 0.004, label_y, "p < 0.05, uncorrected", transform=ax.get_xaxis_transform(),
            fontsize=7.5, color=S.MUTED, va="bottom", ha="left")


def strip_y(ax) -> None:
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)


def table(title: str, df: pd.DataFrame, floatfmt: str = ".3f") -> None:
    TABLES.append(f"### {title}\n\n{df.to_markdown(index=False, floatfmt=floatfmt)}\n")


# ── Data loading ──────────────────────────────────────────────────────────
def load_af3() -> pd.DataFrame:
    return pd.read_csv(FIGDATA / "af3_per_model.csv")


def load_t156_reps() -> pd.DataFrame:
    """Per-trajectory T156 means (standardized re-analysis, 1 ns sampling)."""
    d = pd.read_csv(MDA / "diagnostic_f156_convergence.csv")
    d["variant"] = d["variant"].replace(MD_NAME)
    d["cls"] = d["variant"].map(S.VARIANT_CLASS)
    d["used"] = d["ns"] >= 180
    return d


# ══════════════════════════════════════════════════════════════════════════
#  Structure (README)
# ══════════════════════════════════════════════════════════════════════════
# Label position (render pixels) for each anchor written by render_structure.py.
STRUCTURE_LABELS = {
    "D210E": (735, 404, "right", "D210E†"), "T187S": (1085, 318, "center", None),
    "S185N": (1255, 318, "center", None), "R184C": (1425, 318, "center", None),
    "T156": (1290, 690, "right", None),
    "Y154H": (1905, 655, "left", None), "R160W": (1905, 735, "left", None),
    "G204D": (500, 735, "right", None), "C206F": (900, 1100, "center", None),
    "C208A": (1150, 1100, "center", "C208A†"),
}


def fig_structure() -> None:
    img = plt.imread(FIGDATA / "structure_render.png")
    anchors = pd.read_csv(FIGDATA / "structure_anchors.csv").set_index("label")
    fig = plt.figure(figsize=(8.6, 5.6))
    ax = fig.add_axes([0.0, 0.035, 1.0, S.top(fig, 1.0) - 0.035])
    S.header(fig, "The AlphaFold 3 starting structure",
             "ZFP36L2 tandem zinc fingers (151–220) on the ARE 5′-UUAUUUAUU-3′ with two "
             "Zn$^{2{+}}$. Spheres mark variant sites (Cα).")
    ax.imshow(img, interpolation="lanczos")
    ax.set_xlim(330, 2080)
    ax.set_ylim(1150, 250)
    ax.axis("off")
    for key, (lx, ly, ha, text) in STRUCTURE_LABELS.items():
        r = anchors.loc[key]
        label = text or key
        ax.annotate(label, xy=(r["x_px"], r["y_px"]), xytext=(lx, ly), ha=ha, va="center",
                    fontsize=8.6, color=S.INK, fontweight="bold" if key == "T156" else None,
                    arrowprops=dict(arrowstyle="-", color=S.MUTED, lw=0.7, shrinkA=3,
                                    shrinkB=17))
    for key, (dx, dy) in {"5'": (-55, 0), "3'": (80, 0)}.items():
        r = anchors.loc[key]
        ax.text(r["x_px"] + dx, r["y_px"] + dy, key.replace("'", "′"), ha="center",
                va="center", fontsize=9.5, color=S.INK, fontweight="bold")
    ax.text(1735, 405, "finger 1\n159–178", ha="center", va="center", fontsize=8.6,
            color=S.INK2, linespacing=1.25)
    ax.text(620, 560, "finger 2\n197–216", ha="center", va="center", fontsize=8.6,
            color=S.INK2, linespacing=1.25)
    ax.text(1330, 1100, "RNA", ha="left", va="center", fontsize=8.6, color=S.INK2)
    legend(fig, class_items(["deleterious", "neutral", "zn_ctrl"])
           + [(S.MD, "o", "T156 (best MD feature)"), ("#dfe2e6", "o", "Zn$^{2{+}}$")],
           top_in=0.82)
    fig.text(0.012, 0.012, "† C208A and D210E follow Akçeşme et al. 2025; in the construct "
             "numbering used here they mutate C212 and D219.", fontsize=7.6, color=S.MUTED)
    S.save(fig, "fig1_structure")


# ══════════════════════════════════════════════════════════════════════════
#  AF3
# ══════════════════════════════════════════════════════════════════════════
AF3_METRICS = [
    ("pair_ipTM_protein_rna", "pair-ipTM", "higher = more confident", -1),
    ("pae_cross_mean", "Cross-chain PAE (Å)", "lower = more confident", +1),
    ("pLDDT", "Mean pLDDT", "higher = more confident", -1),
    ("pLDDT_min", "Min pLDDT", "lowest residue: E220, C-terminus", -1),
]


def fig_af3_confidence(af3: pd.DataFrame) -> None:
    sub = af3[af3["variant"].isin(S.DISPLAY_ORDER)]
    seeds = sub.groupby(["variant", "run"])[[m for m, *_ in AF3_METRICS]].mean().reset_index()
    means = seeds.drop(columns="run").groupby("variant").mean()
    rows = S.DISPLAY_ORDER
    ypos = {v: len(rows) - 1 - i for i, v in enumerate(rows)}
    rng = np.random.default_rng(7)

    fig, axes = plt.subplots(1, 4, figsize=(8.6, 5.6), sharey=True)
    fig.subplots_adjust(left=0.105, right=0.985, top=S.top(fig, 1.95),
                        bottom=0.35 / 5.6, wspace=0.16)
    S.header(fig, "AlphaFold 3 is almost equally confident about every variant",
             "Eight AF3 jobs per construct: small dots = one seed (mean of its 5 models), "
             "large dot = construct mean, line = WT.\nOnly C206F, which removes a zinc ligand "
             "of finger 2, separates from wild type on every metric.")
    legend(fig, class_items(["WT", "deleterious", "neutral", "zn_ctrl"]), top_in=1.0)

    out = []
    for ax, (col, name, hint, sign) in zip(axes, AF3_METRICS):
        ax.axhspan(ypos["C206F"] - 0.5, ypos["Y154H"] + 0.5, color=S.WASH, zorder=0, lw=0)
        ax.axvline(means.loc["WT", col], color=S.INK2, lw=0.8, alpha=0.55, zorder=1)
        for v in rows:
            cls = S.VARIANT_CLASS[v]
            vals = seeds.loc[seeds["variant"] == v, col].to_numpy()
            jitter = rng.uniform(-0.17, 0.17, len(vals))
            ax.scatter(vals, ypos[v] + jitter, s=13, color=S.CLASS_COLOR[cls],
                       alpha=0.45, lw=0, zorder=2)
            dot(ax, means.loc[v, col], ypos[v], cls, s=62)
        a = auroc(means[col].to_dict(), sign)
        ax.set_title(f"{name}\n\n", fontsize=9.8)
        ax.text(0, 1.012, f"{hint}\nAUROC {a:.2f}", transform=ax.transAxes,
                fontsize=7.8, color=S.MUTED, va="bottom", linespacing=1.35)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        S.vgrid(ax)
        strip_y(ax)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        out.append((name, a))
    axes[0].set_yticks([ypos[v] for v in rows], rows)
    S.save(fig, "figA2_af3_confidence")

    t = means.loc[rows].reset_index()
    t.insert(1, "class", t["variant"].map(S.VARIANT_CLASS))
    table("AF3 confidence per construct (mean of 40 models = 8 seeds × 5)", t.round(3))
    table("AF3 AUROC, 5 deleterious vs 3 neutral (all-model means)",
          pd.DataFrame(out, columns=["metric", "AUROC"]))


def fig_af3_rna_control(af3: pd.DataFrame) -> None:
    ctrl = pd.read_csv(ROOT / "analysis" / "interface_metrics_models.csv")
    nonbind = ctrl.loc[ctrl["ComplexBase"].str.contains("nonbind"), "pair_ipTM_protein_rna"]
    wt = af3.loc[af3["variant"] == "WT", "pair_ipTM_protein_rna"]
    labeled = af3[af3["variant"].isin(S.DELETERIOUS + S.NEUTRAL)]
    var_means = labeled.groupby("variant")["pair_ipTM_protein_rna"].mean()
    rng = np.random.default_rng(3)

    fig, ax = plt.subplots(figsize=(8.6, 3.7))
    fig.subplots_adjust(left=0.3, right=0.975, top=S.top(fig, 1.45), bottom=0.6 / 3.7)
    S.header(fig, "AF3 notices a bad RNA, not a bad protein",
             "Swapping the AU-rich element for a non-binding nonamer halves pair-ipTM. "
             "Point mutations move it a few hundredths,\nand deleterious and neutral variants "
             "overlap completely.")
    legend(fig, class_items(["WT", "deleterious", "neutral"]), top_in=1.0)

    rows = [
        ("WT protein + non-binding RNA", "UACAAACAU · 1 seed, 5 models", nonbind, "WT", None),
        ("WT protein + ARE", "UUAUUUAUU · 8 seeds, 40 models", wt, "WT", None),
        ("5 deleterious variants + ARE", "200 models · large dots = variant means",
         labeled.loc[labeled["variant"].isin(S.DELETERIOUS), "pair_ipTM_protein_rna"],
         "deleterious", S.DELETERIOUS),
        ("3 neutral variants + ARE", "120 models · large dots = variant means",
         labeled.loc[labeled["variant"].isin(S.NEUTRAL), "pair_ipTM_protein_rna"],
         "neutral", S.NEUTRAL),
    ]
    for i, (lab, sub, vals, cls, members) in enumerate(rows):
        y = len(rows) - 1 - i
        ax.scatter(vals, y + rng.uniform(-0.16, 0.16, len(vals)), s=14 if members else 22,
                   color=S.CLASS_COLOR[cls], alpha=0.35 if members else 0.6, lw=0, zorder=2)
        if members:
            for v in members:
                dot(ax, var_means[v], y, cls, s=52, z=4)
        else:
            dot(ax, vals.mean(), y, cls, s=74, z=4)
        ax.text(-0.015, y + 0.06, lab, transform=ax.get_yaxis_transform(), ha="right",
                va="bottom", fontsize=9, color=S.INK)
        ax.text(-0.015, y - 0.04, sub, transform=ax.get_yaxis_transform(), ha="right",
                va="top", fontsize=7.6, color=S.MUTED)

    drop = wt.mean() - nonbind.mean()
    spread = var_means.max() - var_means.min()
    ax.annotate("", xy=(nonbind.mean() + 0.01, 2.55), xytext=(wt.mean() - 0.01, 2.55),
                arrowprops=dict(arrowstyle="-|>", color=S.INK2, lw=1, mutation_scale=9))
    ax.text((wt.mean() + nonbind.mean()) / 2, 2.46, f"pair-ipTM drops by {drop:.2f}",
            ha="center", va="top", fontsize=8.5, color=S.INK2)
    ax.plot([var_means.min(), var_means.max()], [-0.5, -0.5], color=S.INK2, lw=1)
    for x in (var_means.min(), var_means.max()):
        ax.plot([x, x], [-0.56, -0.44], color=S.INK2, lw=1)
    ax.text(var_means.min() - 0.01, -0.5, f"all 8 variant means within {spread:.2f}",
            ha="right", va="center", fontsize=8, color=S.INK2)
    ax.set_xlim(0.2, 0.86)
    ax.set_ylim(-0.8, len(rows) - 0.35)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("AF3 pair-ipTM, protein–RNA interface")
    S.vgrid(ax)
    S.save(fig, "figA3_af3_rna_control")
    table("AF3 RNA control (pair-ipTM)", pd.DataFrame({
        "complex": ["WT + ARE (UUAUUUAUU)", "WT + non-binding (UACAAACAU)",
                    "lowest variant mean", "highest variant mean"],
        "pair_ipTM": [wt.mean(), nonbind.mean(), var_means.min(), var_means.max()],
        "n_models": [len(wt), len(nonbind), 40, 40]}))


# ══════════════════════════════════════════════════════════════════════════
#  MD
# ══════════════════════════════════════════════════════════════════════════
SOURCES = [("expanse_orig", "Expanse V100", "Replicate 1"),
           ("anvil_rep1", "Anvil A100", "Replicate 2"),
           ("anvil_rep2", "Anvil A100", "Replicate 3"),
           ("expanse_rep2", "Expanse V100", "Replicate 4")]


def fig_md_sampling(reps: pd.DataFrame) -> None:
    rows = S.DISPLAY_ORDER
    used = reps[reps["used"]]
    fig, axes = plt.subplots(1, 4, figsize=(8.6, 4.7), sharey=True)
    fig.subplots_adjust(left=0.085, right=0.92, top=S.top(fig, 1.62), bottom=0.62 / 4.7,
                        wspace=0.1)
    S.header(fig, f"{len(reps)} production runs, {len(used)} used in the analysis",
             f"OpenMM · AMBER14SB + RNA.OL3 · TIP3P, 0.15 M NaCl · 310 K, NPT. Runs of at least "
             f"180 ns count, capped at 200 ns:\n{used['ns'].clip(upper=200).sum() / 1000:.1f} µs "
             f"analyzed. Runs that stopped early (light bars) are left out of every average.")
    handles = [plt.Rectangle((0, 0), 1, 1, color=S.INK2, label="used (≥ 180 ns)"),
               plt.Rectangle((0, 0), 1, 1, color=S.AXIS, label="partial, excluded")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.012, S.top(fig, 1.0)),
               ncol=2, handlelength=1.4, handleheight=0.7, columnspacing=1.6, borderaxespad=0)
    ypos = {v: len(rows) - 1 - i for i, v in enumerate(rows)}
    for ax, (src, platform, label) in zip(axes, SOURCES):
        s = reps[reps["source"] == src].set_index("variant")
        for v in rows:
            if v not in s.index:
                continue
            ns, ok = min(s.loc[v, "ns"], 200), s.loc[v, "used"]
            ax.barh(ypos[v], ns, height=0.56, color=S.INK2 if ok else S.AXIS, lw=0)
            if not ok:
                ax.text(ns + 6, ypos[v], f"{ns:.0f} ns", va="center", fontsize=7.3, color=S.MUTED)
        ax.set_xlim(0, 215)
        ax.set_xticks([0, 100, 200])
        ax.set_title(f"{label}\n", fontsize=9.8)
        ax.text(0, 1.012, platform, transform=ax.transAxes, fontsize=7.8, color=S.MUTED,
                va="bottom")
        strip_y(ax)
        S.vgrid(ax)
    axes[0].set_yticks([ypos[v] for v in rows], rows)
    axes[0].set_ylim(-0.6, len(rows) - 0.4)
    fig.text(0.5, 0.025, "analyzed production time per run (ns)", ha="center", fontsize=9,
             color=S.INK2)
    n_used = used.groupby("variant").size()
    for v in rows:
        axes[-1].text(1.06, ypos[v], f"{n_used.get(v, 0)}",
                      transform=axes[-1].get_yaxis_transform(), va="center", ha="left",
                      fontsize=9, color=S.INK, fontweight="bold")
    axes[-1].text(1.06, len(rows) - 0.35, "used", transform=axes[-1].get_yaxis_transform(),
                  va="bottom", ha="left", fontsize=7.8, color=S.MUTED)
    S.save(fig, "figM1_md_sampling")
    t = reps.pivot(index="variant", columns="source", values="ns").loc[
        rows, [s for s, *_ in SOURCES]]
    table("Analyzed ns per run (≥ 180 ns = used)", t.clip(upper=200).reset_index(),
          floatfmt=".0f")


def load_auroc_table() -> pd.DataFrame:
    a = pd.read_csv(MDA / "auroc_no_c212a.csv")
    q = pd.read_csv(MDA / "mmgbsa_qscore_averaged.csv").set_index("variant")
    extra = pd.DataFrame([
        {"feature": "qscore", "AUROC": auroc(q["qscore_last20pct"].to_dict(), -1),
         "origin": "MD-Qscore"},
        {"feature": "mmgbsa", "AUROC": auroc(q["mmgbsa_mean_kJmol"].to_dict(), +1),
         "origin": "MD-energy"},
    ])
    a = pd.concat([a, extra], ignore_index=True)
    a["method"] = np.where(a["origin"] == "AF3", "AF3", "MD")
    a["label"] = a["feature"].map(FEATURE_LABELS)
    return a.sort_values("AUROC", ascending=True).reset_index(drop=True)


def fig_auroc_ranking(t156_std: float) -> None:
    a = load_auroc_table()
    n = len(a)
    crit = critical_auroc(5, 3, 0.05)
    fig, ax = plt.subplots(figsize=(8.6, 8.9))
    fig.subplots_adjust(left=0.34, right=0.955, top=S.top(fig, 1.6), bottom=0.6 / 8.9)
    S.header(fig, "One MD feature edges out AF3, and nothing reaches significance",
             f"Replicate-averaged AUROC for all {n} features screened (5 deleterious vs 3 neutral), "
             f"direction fixed in advance.\nBelow 0.5 = leaned the wrong way. Open circle = T156 "
             f"after the standardized re-analysis ({t156_std:.3f}).")
    legend(fig, METHOD_ITEMS, top_in=1.0)
    significance_band(ax, crit)
    ax.axvline(0.5, color=S.INK2, lw=0.9, zorder=1)
    ax.text(0.5, 1.01, "chance", transform=ax.get_xaxis_transform(), fontsize=7.5,
            color=S.MUTED, ha="center", va="bottom")
    for i, r in a.iterrows():
        c = S.AF3 if r["method"] == "AF3" else S.MD
        ax.plot([0.5, r["AUROC"]], [i, i], color=c, lw=1.4, alpha=0.55, zorder=2,
                solid_capstyle="butt")
        ax.scatter(r["AUROC"], i, s=40, color=c, edgecolors=S.SURFACE, lw=1.2, zorder=3)
        if r["AUROC"] >= 0.7 or r["AUROC"] <= 0.25:
            right = r["AUROC"] >= 0.5
            ax.text(r["AUROC"] + (0.018 if right else -0.018), i, f"{r['AUROC']:.3f}",
                    va="center", ha="left" if right else "right", fontsize=8, color=S.INK,
                    fontweight="bold" if right else None)
    top = a.index[a["feature"] == "hbond_RNA_F156"][0]
    ax.scatter(t156_std, top, s=44, facecolors=S.SURFACE, edgecolors=S.MD, lw=1.5, zorder=4)
    ax.set_yticks(range(n), a["label"])
    for tick, (_, r) in zip(ax.get_yticklabels(), a.iterrows()):
        if r["feature"] in ("hbond_RNA_F156", "pair_ipTM_protein_rna"):
            tick.set_color(S.INK)
            tick.set_fontweight("bold")
    ax.set_xlim(0.1, 1.02)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel("AUROC (replicate-averaged, C208A control excluded)")
    strip_y(ax)
    S.vgrid(ax)
    S.save(fig, "figM2_auroc_ranking")
    t = a.sort_values("AUROC", ascending=False)[["label", "method", "AUROC", "feature"]]
    t = t.rename(columns={"feature": "column in analysis CSVs"})
    table("Feature AUROCs (5 deleterious vs 3 neutral)", t)


def fig_t156_replicates(reps: pd.DataFrame) -> None:
    used = reps[reps["used"]]
    means = used.groupby("variant")["f156_mean"].mean().sort_values()
    rows = list(means.index)
    ypos = {v: i for i, v in enumerate(rows)}
    a = auroc(means.to_dict(), -1)
    p = tail_p(5, 3)[int(round(a * 15))]
    del_mean = means[S.DELETERIOUS].mean()
    neu_mean = means[S.NEUTRAL].mean()

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    fig.subplots_adjust(left=0.105, right=0.86, top=S.top(fig, 1.55), bottom=0.62 / 5.0)
    S.header(fig, "The best MD signal is smaller than the noise between replicates",
             f"T156–RNA hydrogen bonds per frame, one small dot per 200 ns replicate. AUROC {a:.3f} "
             f"({int(round(a * 15))} of 15 pairs), one-sided p = {p:.3f}.\nClass means differ by "
             f"{neu_mean - del_mean:.3f}; replicates of a single variant spread 0.10–0.26.")
    legend(fig, class_items(["WT", "deleterious", "neutral", "zn_ctrl"]), top_in=1.0)
    ax.axvspan(del_mean, neu_mean, color=S.WASH, zorder=0, lw=0)
    ax.axvline(del_mean, color=S.DEL, lw=1.0, alpha=0.85, zorder=1)
    ax.axvline(neu_mean, color=S.NEU, lw=1.0, alpha=0.85, zorder=1)
    ax.text(del_mean - 0.003, 1.01, "deleterious mean", transform=ax.get_xaxis_transform(),
            ha="right", va="bottom", fontsize=7.5, color=S.INK2)
    ax.text(neu_mean + 0.003, 1.01, "neutral mean", transform=ax.get_xaxis_transform(),
            ha="left", va="bottom", fontsize=7.5, color=S.INK2)
    for v in rows:
        cls = S.VARIANT_CLASS[v]
        vals = used.loc[used["variant"] == v, "f156_mean"].to_numpy()
        ax.plot([vals.min(), vals.max()], [ypos[v]] * 2, color=S.CLASS_COLOR[cls], lw=1.3,
                alpha=0.35, zorder=2, solid_capstyle="round")
        ax.scatter(vals, [ypos[v]] * len(vals), s=30, color=S.CLASS_COLOR[cls], alpha=0.55,
                   lw=0, zorder=3, marker=S.CLASS_MARKER[cls])
        dot(ax, means[v], ypos[v], cls, s=74, z=4)
        ax.text(1.015, ypos[v], f"{means[v]:.3f}", transform=ax.get_yaxis_transform(),
                va="center", fontsize=8.5, color=S.INK)
        ax.text(1.105, ypos[v], f"n={len(vals)}", transform=ax.get_yaxis_transform(),
                va="center", fontsize=7.8, color=S.MUTED)
    ax.text(1.015, len(rows) - 0.4, "mean", transform=ax.get_yaxis_transform(), va="bottom",
            fontsize=7.8, color=S.MUTED)
    ax.set_yticks(range(len(rows)), rows)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlim(0.33, 0.75)
    ax.set_xlabel("T156–RNA hydrogen bonds per frame (1 ns sampling, ≤ 200 ns per run)")
    strip_y(ax)
    S.vgrid(ax)
    S.save(fig, "figM3_t156_replicates")
    t = used.groupby("variant")["f156_mean"].agg(["mean", "min", "max", "count"]).loc[rows[::-1]]
    t.insert(0, "class", [S.VARIANT_CLASS[v] for v in t.index])
    table("T156–RNA H-bonds per frame (runs ≥ 180 ns)", t.reset_index())


def fig_t156_convergence(reps: pd.DataFrame) -> None:
    s = pd.read_csv(FIGDATA / "t156_rna_hbond_series_1ns.csv")
    s["variant"] = s["variant"].replace(MD_NAME)
    used_keys = set(map(tuple, reps.loc[reps["used"], ["variant", "source"]].to_numpy()))
    order = ["WT", "Y154H", "R160W", "R184C", "G204D", "C206F", "S185N", "T187S", "D210E", "C208A"]
    means = reps[reps["used"]].groupby("variant")["f156_mean"].mean()
    del_mean, neu_mean = means[S.DELETERIOUS].mean(), means[S.NEUTRAL].mean()

    fig, axes = plt.subplots(2, 5, figsize=(8.6, 5.3), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, top=S.top(fig, 1.55), bottom=0.55 / 5.3,
                        wspace=0.2, hspace=0.34)
    legend(fig, class_items(["WT", "deleterious", "neutral", "zn_ctrl"]), top_in=1.0)
    S.header(fig, "Each replicate settles, but on its own plateau",
             "Running mean of T156–RNA H-bonds for every replicate. Gray band = gap between the "
             "deleterious and neutral class means,\nthe entire signal behind the 0.800 AUROC. "
             "Faint lines are partial runs left out of the averages.")
    for ax, v in zip(axes.flat, order):
        cls = S.VARIANT_CLASS[v]
        ax.axhspan(del_mean, neu_mean, color=S.AXIS, alpha=0.6, lw=0, zorder=0)
        for src, *_ in SOURCES:
            g = s[(s["variant"] == v) & (s["source"] == src)]
            if g.empty:
                continue
            run = np.cumsum(g["hbonds"].to_numpy()) / np.arange(1, len(g) + 1)
            t = g["t_ns"].to_numpy()
            keep = t >= 3
            full = (v, src) in used_keys
            ax.plot(t[keep], run[keep], color=S.CLASS_COLOR[cls], lw=1.5 if full else 1.0,
                    alpha=0.95 if full else 0.3, zorder=2 if full else 1)
        ax.set_title(v, fontsize=9.2, pad=4)
        ax.set_ylim(0.15, 1.0)
        ax.set_xlim(0, 200)
        ax.set_xticks([0, 100, 200])
        S.hgrid(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel("running mean")
    fig.text(0.53, 0.025, "simulation time (ns)", ha="center", fontsize=9, color=S.INK2)
    S.save(fig, "figM4_t156_convergence")


def fig_mechanism_checks() -> None:
    base = pd.read_csv(MDA / "t156_base_reading.csv", index_col=0)
    base.index = [MD_NAME.get(v, v) for v in base.index]
    stack = pd.read_csv(MDA / "stacking_metrics.csv")
    stack = stack[stack["ns"] >= 180]
    stack["variant"] = stack["variant"].replace(MD_NAME)
    aro = ["Y154", "F162", "Y170", "F176", "F180", "F200", "F205", "Y208", "F214"]
    combos = {"aromatic_any": "any aromatic", "aromatic_sum": "all aromatics, summed",
              "finger1_any": "any, finger 1", "finger2_any": "any, finger 2"}
    sm = stack.groupby("variant")[aro + list(combos)].mean()
    st = pd.DataFrame([(c, combos.get(c, c), auroc(sm[c].to_dict(), -1))
                       for c in aro + list(combos)], columns=["feature", "label", "AUROC"])
    ryk = pd.read_csv(MDA / "t156_ryktel_propagation.csv")
    ryk["variant"] = ryk["variant"].replace(MD_NAME)
    ryk = ryk.sort_values("RYKTEL_sum").reset_index(drop=True)
    u6 = base["6"]
    u6_del, u6_neu = u6[S.DELETERIOUS].mean(), u6[S.NEUTRAL].mean()

    fig = plt.figure(figsize=(8.6, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.25, 1, 1], left=0.075, right=0.985,
                          top=S.top(fig, 1.45), bottom=0.85 / 5.0, wspace=0.42)
    S.header(fig, "Three ways the signal could have localized, and didn't",
             "If deleterious variants broke one specific contact, at least one of these "
             "readouts should separate the classes. None does.")
    legend(fig, class_items(["WT", "deleterious", "neutral", "zn_ctrl"]), top_in=0.8)

    ax = fig.add_subplot(gs[0])
    rows = S.DISPLAY_ORDER
    mat = base.loc[rows].to_numpy()
    cmap = LinearSegmentedColormap.from_list("blue", ["#ffffff"] + S.BLUE_RAMP[2:11])
    ax.imshow(mat, cmap=cmap, vmin=0, vmax=0.65, aspect="auto")
    for i, v in enumerate(rows):
        ax.text(5, i, f"{mat[i, 5]:.2f}", ha="center", va="center", fontsize=7.4,
                color="#ffffff" if mat[i, 5] > 0.4 else S.INK)
    ax.set_xticks(range(9), [f"{b}{i + 1}" for i, b in enumerate("UUAUUUAUU")], fontsize=7.6)
    ax.set_yticks(range(len(rows)), rows, fontsize=8)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("a  T156 reads one base, U6", fontsize=9.6)
    ax.text(0, -0.09, f"contact occupancy per base\ndeleterious {u6_del:.2f} vs "
            f"neutral {u6_neu:.2f} at U6", transform=ax.transAxes, fontsize=7.4,
            color=S.MUTED, va="top", linespacing=1.35)

    ax = fig.add_subplot(gs[1])
    st = st.sort_values("AUROC").reset_index(drop=True)
    ax.axvline(0.5, color=S.INK2, lw=0.9, zorder=1)
    for i, r in st.iterrows():
        ax.plot([0.5, r["AUROC"]], [i, i], color=S.INK2, lw=1.3, alpha=0.4, zorder=2)
        ax.scatter(r["AUROC"], i, s=34, color=S.INK2, edgecolors=S.SURFACE, lw=1.1, zorder=3)
    ax.set_yticks(range(len(st)), st["label"], fontsize=7.8)
    ax.set_xlim(0.25, 0.95)
    ax.set_xticks([0.3, 0.5, 0.7, 0.9])
    strip_y(ax)
    S.vgrid(ax)
    ax.set_title("b  Aromatic stacking", fontsize=9.6)
    ax.set_xlabel("AUROC", fontsize=8.5)
    ax.text(0, -0.155, f"best {st['AUROC'].max():.2f}; 13 readouts tested",
            transform=ax.transAxes, fontsize=7.4, color=S.MUTED, va="top")

    ax = fig.add_subplot(gs[2])
    for i, r in ryk.iterrows():
        dot(ax, r["RYKTEL_sum"], i, S.VARIANT_CLASS[r["variant"]], s=46)
    d = ryk[ryk["class"] == "del"]["RYKTEL_sum"].mean()
    n = ryk[ryk["class"] == "neu"]["RYKTEL_sum"].mean()
    ax.axvline(d, color=S.DEL, lw=1, alpha=0.85, zorder=1)
    ax.axvline(n, color=S.NEU, lw=1, alpha=0.85, zorder=1)
    ax.set_yticks(range(len(ryk)), ryk["variant"], fontsize=8)
    strip_y(ax)
    S.vgrid(ax)
    ax.set_title("c  RYKTEL-region contacts", fontsize=9.6)
    ax.set_xlabel("RNA contacts, residues 153–158", fontsize=8.5)
    ax.text(0, -0.155, f"deleterious {d:.2f} > neutral {n:.2f}:\nthe wrong direction",
            transform=ax.transAxes, fontsize=7.4, color=S.MUTED, va="top", linespacing=1.35)
    S.save(fig, "figM5_mechanism_checks")
    table("Aromatic stacking AUROCs", st[["label", "AUROC"]])
    table("RYKTEL-region RNA contacts", ryk[["variant", "class", "RYKTEL_sum"]])


# ══════════════════════════════════════════════════════════════════════════
#  Conclusions
# ══════════════════════════════════════════════════════════════════════════
def fig_n8_ceiling(t156: float, af3_iptm: float) -> int:
    p53 = tail_p(5, 3)
    ks = np.arange(8, 16)
    bonf = 0.05 / N_FEATURES_SCREENED

    fig = plt.figure(figsize=(8.6, 4.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], left=0.085, right=0.975,
                          top=S.top(fig, 1.25), bottom=0.62 / 4.9, wspace=0.3)
    S.header(fig, "Eight labeled variants cannot carry a significant result",
             f"Exact Mann–Whitney null for 5 deleterious vs 3 neutral: 56 equally likely "
             f"rankings. Even a perfect ranking gives p = {p53[15]:.3f},\nwell above the "
             f"Bonferroni bar for {N_FEATURES_SCREENED} screened features "
             f"(0.05 / {N_FEATURES_SCREENED} = {bonf:.4f}).")

    ax = fig.add_subplot(gs[0])
    ax.axhspan(1e-4, bonf, color=S.WASH, lw=0, zorder=0)
    ax.axhline(0.05, color=S.INK2, lw=0.9)
    ax.axhline(bonf, color=S.INK, lw=0.9)
    ax.text(0.508, 0.05 * 1.12, "p = 0.05", fontsize=7.6, color=S.INK2, va="bottom")
    ax.text(0.508, bonf * 1.12, f"Bonferroni, {N_FEATURES_SCREENED} features", fontsize=7.6,
            color=S.INK, va="bottom")
    ax.text(0.508, bonf / 1.3, "unreachable at n = 8", fontsize=7.6, color=S.MUTED, va="top")
    ax.plot(ks / 15, p53[ks], color=S.INK2, lw=1.4, zorder=2, drawstyle="steps-mid", alpha=0.35)
    ax.scatter(ks / 15, p53[ks], s=22, color=S.INK2, zorder=3, edgecolors=S.SURFACE, lw=1)
    for val, c in [(af3_iptm, S.AF3), (t156, S.MD), (1.0, S.INK)]:
        ax.scatter(val, p53[int(round(val * 15))], s=60, color=c, zorder=4,
                   edgecolors=S.SURFACE, lw=1.3)
    ax.text(1.0 - 0.025, p53[15], f"perfect\np = {p53[15]:.3f}", ha="right", va="center",
            fontsize=7.4, color=S.INK)
    keys = [Line2D([], [], ls="", marker="o", markersize=7, markerfacecolor=c,
                   markeredgecolor=S.SURFACE,
                   label=f"{lab}: {val:.3f}, p = {p53[int(round(val * 15))]:.3f}")
            for val, lab, c in [(af3_iptm, "AF3 pair-ipTM", S.AF3), (t156, "T156 (MD)", S.MD)]]
    ax.legend(handles=keys, loc="center left", bbox_to_anchor=(0.0, 0.4), fontsize=7.6,
              handletextpad=0.2, borderaxespad=0.3)
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
    ax.set_ylim(4e-4, 1.5)
    ax.set_xlim(0.5, 1.035)
    ax.set_xlabel("AUROC (5 vs 3: only multiples of 1/15 exist)")
    ax.set_ylabel("one-sided exact p")
    S.hgrid(ax)
    ax.set_title("a  What eight variants allow", fontsize=9.8)

    ax = fig.add_subplot(gs[1])
    per_class = np.arange(3, 31)
    crit_u = np.array([critical_auroc(k, k, 0.05) for k in per_class])
    crit_b = np.array([critical_auroc(k, k, bonf) for k in per_class])
    ax.plot(2 * per_class, crit_u, color=S.INK2, lw=1.6)
    ax.plot(2 * per_class, crit_b, color=S.INK, lw=1.6)
    ax.text(60.5, crit_b[-1] - 0.012, "needed after Bonferroni", ha="right", va="top",
            fontsize=7.6, color=S.INK)
    ax.text(60.5, crit_u[-1] - 0.012, "needed for p < 0.05", ha="right", va="top",
            fontsize=7.6, color=S.INK2)
    ax.axhline(t156, color=S.MD, lw=1.4)
    ax.text(60.5, t156 + 0.008, f"observed T156, {t156:.3f}", ha="right", va="bottom",
            fontsize=7.8, color=S.INK)
    need = int(next(2 * k for k, c in zip(per_class, crit_b) if c <= t156))
    ax.scatter(need, t156, s=60, color=S.MD, edgecolors=S.SURFACE, lw=1.3, zorder=4)
    ax.annotate(f"≈ {need} variants ({need // 2} + {need // 2})\nto clear Bonferroni\n"
                f"at the observed effect", xy=(need, t156), xytext=(need + 3, t156 + 0.075),
                ha="left", va="bottom", fontsize=7.8, color=S.INK, linespacing=1.3,
                arrowprops=dict(arrowstyle="-", color=S.MUTED, lw=0.7, shrinkA=2, shrinkB=4))
    ax.axvline(8, color=S.MUTED, lw=0.8)
    ax.text(8.8, 0.56, "this study:\n8 variants", fontsize=7.6, color=S.MUTED, va="bottom")
    ax.set_xlim(5, 62)
    ax.set_ylim(0.55, 1.02)
    ax.set_xlabel("validated variants (half deleterious, half neutral)")
    ax.set_ylabel("minimum AUROC")
    S.hgrid(ax)
    ax.set_title("b  What it would take", fontsize=9.8)
    S.save(fig, "figC1_n8_ceiling")
    table("Exact one-sided p for 5 vs 3", pd.DataFrame(
        {"pairs_correct": ks, "AUROC": ks / 15, "p": p53[ks]}), floatfmt=(".0f", ".3f", ".4f"))
    table("Minimum AUROC by design size (balanced)", pd.DataFrame(
        {"variants_total": 2 * per_class, "AUROC_p<0.05": crit_u,
         f"AUROC_p<0.05/{N_FEATURES_SCREENED}": crit_b}), floatfmt=(".0f", ".3f", ".3f"))
    return need


def fig_forking_paths(af3: pd.DataFrame, t156_std: float) -> None:
    def window(n_ns: str) -> float:
        h = pd.read_csv(MDA / f"hbonds_per_residue_common{n_ns}.csv", index_col=0)
        return auroc(h.loc[6].to_dict(), -1)  # resSeq 6 = T156

    multi = pd.read_csv(MDA / "auroc_no_c212a.csv").set_index("feature")["AUROC"]
    top = af3[af3["model"] == 0].groupby("variant").mean(numeric_only=True)
    allm = af3.groupby("variant").mean(numeric_only=True)
    paths = [
        ("MD", "T156 · 1 run per variant, first 139 ns", window("139")),
        ("MD", "T156 · 1 run per variant, first 269 ns", window("269")),
        ("MD", "T156 · 2–4 runs, 200 ns, 0.1 ns sampling", multi["hbond_RNA_F156"]),
        ("MD", "T156 · 2–4 runs, 200 ns, 1 ns sampling, ≥ 180 ns", t156_std),
        ("AF3", "pair-ipTM · mean of all 40 models",
         auroc(allm["pair_ipTM_protein_rna"].to_dict(), -1)),
        ("AF3", "pair-ipTM · top-ranked model per seed",
         auroc(top["pair_ipTM_protein_rna"].to_dict(), -1)),
        ("AF3", "min pLDDT · mean of all 40 models", auroc(allm["pLDDT_min"].to_dict(), -1)),
        ("AF3", "min pLDDT · top-ranked model per seed", auroc(top["pLDDT_min"].to_dict(), -1)),
    ]
    crit = critical_auroc(5, 3, 0.05)
    fig, ax = plt.subplots(figsize=(8.6, 4.6))
    fig.subplots_adjust(left=0.43, right=0.96, top=S.top(fig, 1.5), bottom=0.62 / 4.6)
    S.header(fig, "Reasonable analysis choices move the answer more than biology does",
             "Same eight labeled variants, AUROC under defensible alternatives. Both values "
             "that look significant rest on a fragile\nchoice: a 139 ns window, or the "
             "confidence of the floppy C-terminal residue.")
    legend(fig, [(S.MD, "o", "MD, T156–RNA H-bonds"), (S.AF3, "o", "AlphaFold 3 confidence")],
           top_in=1.0)
    significance_band(ax, crit)
    ax.axvline(0.5, color=S.INK2, lw=0.9)
    ax.text(0.5, 1.01, "chance", transform=ax.get_xaxis_transform(), fontsize=7.5,
            color=S.MUTED, ha="center", va="bottom")
    ys = list(range(len(paths)))[::-1]
    for y, (meth, lab, val) in zip(ys, paths):
        c = S.AF3 if meth == "AF3" else S.MD
        ax.plot([0.5, val], [y, y], color=c, lw=1.4, alpha=0.55, solid_capstyle="butt")
        ax.scatter(val, y, s=46, color=c, edgecolors=S.SURFACE, lw=1.2, zorder=3)
        ax.text(val + 0.015, y, f"{val:.3f}", va="center", fontsize=8.2, color=S.INK)
    ax.axhline(3.5, color=S.GRID, lw=0.8)
    ax.set_yticks(ys, [lab for _, lab, _ in paths])
    strip_y(ax)
    ax.set_xlim(0.45, 1.02)
    ax.set_ylim(-0.6, len(paths) - 0.4)
    ax.set_xlabel("AUROC, 5 deleterious vs 3 neutral")
    S.vgrid(ax)
    S.save(fig, "figC2_forking_paths")
    table("AUROC under alternative analysis choices",
          pd.DataFrame(paths, columns=["method", "analysis", "AUROC"]))


def main() -> None:
    S.apply()
    af3 = load_af3()
    reps = load_t156_reps()
    used = reps[reps["used"]]
    t156_std = auroc(used.groupby("variant")["f156_mean"].mean().to_dict(), -1)
    af3_iptm = auroc(af3.groupby("variant")["pair_ipTM_protein_rna"].mean().to_dict(), -1)
    assert abs(t156_std - 0.8) < 1e-9 and abs(af3_iptm - 11 / 15) < 1e-9

    fig_structure()
    fig_af3_confidence(af3)
    fig_af3_rna_control(af3)
    fig_md_sampling(reps)
    fig_auroc_ranking(t156_std)
    fig_t156_replicates(reps)
    fig_t156_convergence(reps)
    fig_mechanism_checks()
    need = fig_n8_ceiling(t156_std, af3_iptm)
    fig_forking_paths(af3, t156_std)
    (FIGDATA / "tables.md").write_text(
        "# Numbers behind the figures\n\nGenerated by figures/src/make_figures.py.\n\n"
        + "\n".join(TABLES))
    print(f"T156 standardized AUROC {t156_std:.3f}; AF3 pair-ipTM {af3_iptm:.3f}; "
          f"variants needed for Bonferroni at observed effect: {need}")


if __name__ == "__main__":
    main()
