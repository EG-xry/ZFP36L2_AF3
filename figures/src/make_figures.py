#!/usr/bin/env python3
"""
Build every data figure in the write-up from files committed to this repo.

    python figures/src/make_figures.py        (repo .venv: numpy, pandas, matplotlib)

Inputs
    figures/data/af3_per_model.csv        af3_metrics.py (from Data/)
    figures/data/md_run_inventory.csv     audit_md_runs.py (protein actually simulated per run)
    expanse_md_analysis/*.csv             the project's MD analysis outputs
    analysis/interface_metrics_models.csv non-binding RNA control (AF3)

Outputs: figures/<name>.svg (editable text) + figures/<name>.png, and
figures/data/tables.md with the numbers behind each figure.

Which MD runs count: every Anvil replicate simulated the wild-type protein,
whatever its folder says (see audit_md_runs.py). Variant readouts therefore
use only runs whose topology carries the variant's mutation ("genuine"
runs); the wild-type runs, misfiled ones included, form a wild-type noise set.

Naming: the MD code calls resSeq 6 "F156", but that residue is Thr156 of the
RYKTEL motif. Five other per-residue columns carry the same kind of shifted
name; FEATURE_LABELS below uses the residue actually present at each position.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator

import style as S

ROOT = Path(__file__).resolve().parents[2]
MDA = ROOT / "expanse_md_analysis"
FIGDATA = ROOT / "figures" / "data"

# MD analysis code names C208A "C212A" (the residue it mutates in construct numbering).
MD_NAME = {"C212A": "C208A"}
STACKING = ["Y154", "F162", "Y170", "F176", "F180", "F200", "F205", "Y208", "F214",
            "aromatic_any", "aromatic_sum", "finger1_any", "finger2_any"]
# 32 readouts in the screen (figure M2) + 13 aromatic-stacking readouts.
N_FEATURES_SCREENED = 32 + len(STACKING)

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
# Direction scored as deleterious (as in multi_replicate_auroc_v2.py).
MD_SIGN = {"rmsd_ca_last20pct": 1, "rmsd_zf1_last20pct": 1, "rmsd_zf2_last20pct": 1,
           "rmsd_rnaP_last20pct": 1, "contacts_last20pct": -1, "coulomb_prot_rna": 1,
           "salt_bridges_prot_rna": -1, "interface_atom_count_5A": -1,
           "binding_residue_rmsf_A": 1, "qscore": -1, "mmgbsa": 1}

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


def exact_p(a: float) -> float:
    """One-sided exact p for an AUROC on 5 deleterious vs 3 neutral."""
    return float(tail_p(5, 3)[int(round(a * 15))])


def wt_noise_p(runs: pd.DataFrame, col: str, sign: int, b: int = 20000) -> float:
    """How often wild-type runs, dealt out as fake variants with the genuine run counts,
    rank at least as well as the real variants do."""
    labeled = runs[runs["genuine"] & runs["label"].isin(S.DELETERIOUS + S.NEUTRAL)]
    counts = labeled.groupby("label").size()
    observed = auroc(labeled.groupby("label")[col].mean().to_dict(), sign)
    pool = runs.loc[runs["wt_protein"], col].to_numpy()
    rng = np.random.default_rng(0)
    hits = 0
    for _ in range(b):
        shuffled, i, fake = rng.permutation(pool), 0, {}
        for v in S.DELETERIOUS + S.NEUTRAL:
            fake[v] = shuffled[i:i + counts[v]].mean()
            i += counts[v]
        hits += auroc(fake, sign) >= observed - 1e-12
    return hits / b


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


def significance_band(ax, crit: float, label_y: float = 1.01) -> None:
    """Wash over AUROC values that would reach uncorrected p < 0.05; label above the axes."""
    lo = crit - 1 / 30  # halfway to the next attainable value below (step = 1/15)
    ax.axvspan(lo, 1.02, color=S.WASH, lw=0, zorder=0)
    ax.text(1.02, label_y, "p < 0.05, uncorrected", transform=ax.get_xaxis_transform(),
            fontsize=7.5, color=S.MUTED, va="bottom", ha="right")


def strip_y(ax) -> None:
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)


def table(title: str, df: pd.DataFrame, floatfmt: str = ".3f") -> None:
    TABLES.append(f"### {title}\n\n{df.to_markdown(index=False, floatfmt=floatfmt)}\n")


# ── Data loading ──────────────────────────────────────────────────────────
def load_af3() -> pd.DataFrame:
    return pd.read_csv(FIGDATA / "af3_per_model.csv")


def with_audit(df: pd.DataFrame, ns_col: str | None = None) -> pd.DataFrame:
    """Attach what each run actually simulated (md_run_inventory.csv) to a per-run table.

    Adds: label (display name), wt_protein (topology carries the WT protein),
    genuine (topology carries the folder's own mutation), used (run long enough).
    """
    inv = pd.read_csv(FIGDATA / "md_run_inventory.csv")
    inv["variant"] = inv["folder"].replace({"C208A": "C212A"})
    out = df.merge(inv[["source", "variant", "protein"]], on=["source", "variant"],
                   how="left", validate="one_to_one")
    assert out["protein"].notna().all(), "run missing from md_run_inventory.csv"
    out["label"] = out["variant"].replace(MD_NAME)
    out["wt_protein"] = out["protein"] == "WT"
    out["genuine"] = ~out["wt_protein"]
    if ns_col == "ns":
        out["used"] = out["ns"] >= 180
    elif ns_col == "n_frames":
        out["used"] = out["n_frames"] >= 1800
    else:
        out["used"] = True
    return out


def load_t156_runs() -> pd.DataFrame:
    """Per-run T156 means (1 ns sampling, <= 200 ns) with the audit attached."""
    return with_audit(pd.read_csv(MDA / "diagnostic_f156_convergence.csv"), "ns")


def load_stacking_runs() -> pd.DataFrame:
    s = pd.read_csv(MDA / "stacking_metrics.csv").rename(columns={"cls": "class"})
    return with_audit(s, "ns")


# ══════════════════════════════════════════════════════════════════════════
#  Structure (README)
# ══════════════════════════════════════════════════════════════════════════
# Label position (render pixels) for each anchor written by render_structure.py.
STRUCTURE_LABELS = {
    "D210E": (735, 404, "right", "D210E†"), "T187S": (1085, 318, "center", None),
    "S185N": (1255, 318, "center", None), "R184C": (1425, 318, "center", None),
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
        ax.annotate(text or key, xy=(r["x_px"], r["y_px"]), xytext=(lx, ly), ha=ha,
                    va="center", fontsize=8.6, color=S.INK,
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
           + [("#dfe2e6", "o", "Zn$^{2{+}}$")], top_in=0.82)
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


def protein_class(protein: str) -> str:
    """Class of the protein a run actually simulated (not of its folder)."""
    if protein == "WT":
        return "WT"
    if protein == "C212A":
        return "zn_ctrl"
    if protein == "D219E":
        return "neutral"
    return S.VARIANT_CLASS[protein]


def fig_md_sampling(runs: pd.DataFrame) -> None:
    rows = S.DISPLAY_ORDER
    used = runs[runs["used"]]
    genuine = used[used["genuine"] & (used["label"] != "WT")]
    wt_set = used[used["wt_protein"]]
    misfiled = wt_set[wt_set["label"] != "WT"]
    fig, axes = plt.subplots(1, 4, figsize=(8.6, 4.9), sharey=True)
    fig.subplots_adjust(left=0.085, right=0.915, top=S.top(fig, 1.62), bottom=0.62 / 4.9,
                        wspace=0.1)
    S.header(fig, "Every Anvil run simulated the wild type",
             f"Bar color = the protein actually in each run's topology, whatever its folder "
             f"says. {len(genuine)} complete runs carry their\nvariant; {len(misfiled)} complete "
             f"runs in variant folders are wild type and join the {len(wt_set) - len(misfiled)} "
             f"real wild-type runs as a {len(wt_set)}-run wild-type set.")
    legend(fig, class_items(["WT", "deleterious", "neutral", "zn_ctrl"])
           + [(S.AXIS, "s", "stopped early (< 180 ns)")], top_in=1.0)
    ypos = {v: len(rows) - 1 - i for i, v in enumerate(rows)}
    for ax, (src, platform, label) in zip(axes, SOURCES):
        s = runs[runs["source"] == src].set_index("label")
        for v in rows:
            if v not in s.index:
                continue
            r = s.loc[v]
            ns = min(r["ns"], 200)
            color = S.CLASS_COLOR[protein_class(r["protein"])]
            ax.barh(ypos[v], ns, height=0.56, color=color if r["used"] else S.AXIS, lw=0)
            if not r["used"]:
                ax.text(ns + 6, ypos[v], f"{ns:.0f} ns", va="center", fontsize=7.3,
                        color=S.MUTED)
            elif r["wt_protein"] and v != "WT":
                ax.text(8, ypos[v], "WT", va="center", ha="left", fontsize=7.3,
                        color="#ffffff", fontweight="bold")
        ax.set_xlim(0, 215)
        ax.set_xticks([0, 100, 200])
        ax.set_title(f"{label}\n", fontsize=9.8)
        ax.text(0, 1.012, platform, transform=ax.transAxes, fontsize=7.8, color=S.MUTED,
                va="bottom")
        strip_y(ax)
        S.vgrid(ax)
    axes[0].set_yticks([ypos[v] for v in rows], [f"{v} folder" if v != "WT" else "WT folder"
                                                 for v in rows], fontsize=8)
    axes[0].set_ylim(-0.6, len(rows) - 0.4)
    fig.text(0.5, 0.025, "analyzed production time per run (ns)", ha="center", fontsize=9,
             color=S.INK2)
    n_gen = genuine.groupby("label").size()
    for v in rows:
        n = len(wt_set) if v == "WT" else n_gen.get(v, 0)
        axes[-1].text(1.07, ypos[v], f"{n}", transform=axes[-1].get_yaxis_transform(),
                      va="center", ha="left", fontsize=9, color=S.INK, fontweight="bold")
    axes[-1].text(1.07, len(rows) - 0.35, "usable", transform=axes[-1].get_yaxis_transform(),
                  va="bottom", ha="left", fontsize=7.8, color=S.MUTED)
    S.save(fig, "figM1_md_sampling")
    t = runs.pivot(index="label", columns="source", values="protein").loc[
        rows, [s for s, *_ in SOURCES]]
    table("Protein actually simulated, by folder and replicate", t.reset_index(), floatfmt="")
    n = runs.pivot(index="label", columns="source", values="ns").loc[rows, [s for s, *_ in SOURCES]]
    table("Analyzed ns per run (≥ 180 ns = complete)", n.clip(upper=200).reset_index(),
          floatfmt=".0f")


def screen_table() -> pd.DataFrame:
    """AUROC for the 32-readout screen: genuine runs only, and as first published."""
    m = with_audit(pd.read_csv(MDA / "multi_replicate_all_metrics.csv"))
    q = with_audit(pd.read_csv(MDA / "mmgbsa_qscore_all.csv"), "n_frames")
    q = q.rename(columns={"qscore_last20pct": "qscore", "mmgbsa_mean_kJmol": "mmgbsa"})
    signs = {**MD_SIGN, **{c: -1 for c in m.columns if c.startswith("hbond_RNA_")}}
    published = pd.read_csv(MDA / "auroc_no_c212a.csv").set_index("feature")["AUROC"]
    q_all = q[q["used"]].groupby("label")[["qscore", "mmgbsa"]].mean()
    for col in ("qscore", "mmgbsa"):
        published[col] = auroc(q_all[col].to_dict(), signs[col])

    rows = []
    for feature, a_pub in published.items():
        if feature not in FEATURE_LABELS:
            continue
        if feature in signs:
            src = q if feature in ("qscore", "mmgbsa") else m
            gen = src[src["used"] & src["genuine"]].groupby("label")[feature].mean()
            a_gen, method = auroc(gen.to_dict(), signs[feature]), "MD"
        else:
            a_gen, method = a_pub, "AF3"  # AF3 readouts do not depend on MD runs
        rows.append({"feature": feature, "label": FEATURE_LABELS[feature], "method": method,
                     "AUROC_genuine": a_gen, "AUROC_published": a_pub})
    return pd.DataFrame(rows).sort_values("AUROC_genuine").reset_index(drop=True)


def fig_auroc_ranking(t: pd.DataFrame) -> None:
    n = len(t)
    crit = critical_auroc(5, 3, 0.05)
    best_md = t[t["method"] == "MD"]["AUROC_genuine"].max()
    fig, ax = plt.subplots(figsize=(8.6, 9.0))
    fig.subplots_adjust(left=0.34, right=0.955, top=S.top(fig, 1.6), bottom=0.6 / 9.0)
    S.header(fig, "With only genuine runs, no MD readout beats AlphaFold 3",
             f"AUROC for all {n} screened readouts, 5 deleterious vs 3 neutral, direction fixed "
             f"in advance. Filled = genuine variant runs only\n(1–2 per variant). Open = as first "
             f"computed, with 16 wild-type runs averaged into the variants. Best MD: {best_md:.3f}.")
    legend(fig, [(S.MD, "o", "MD readout, genuine runs"), (S.AF3, "o", "AlphaFold 3 confidence"),
                 (S.SURFACE, "o", "MD, as first computed")], top_in=1.0)
    fig.legends[-1].legend_handles[2].set_markeredgecolor(S.MUTED)
    significance_band(ax, crit)
    ax.axvline(0.5, color=S.INK2, lw=0.9, zorder=1)
    ax.text(0.5, 1.01, "chance", transform=ax.get_xaxis_transform(), fontsize=7.5,
            color=S.MUTED, ha="center", va="bottom")
    for i, r in t.iterrows():
        c = S.AF3 if r["method"] == "AF3" else S.MD
        a = r["AUROC_genuine"]
        if r["method"] == "MD" and abs(r["AUROC_published"] - a) > 1e-9:
            ax.plot([r["AUROC_published"], a], [i, i], color=S.AXIS, lw=1.0, zorder=1)
            ax.scatter(r["AUROC_published"], i, s=34, facecolors=S.SURFACE,
                       edgecolors=S.MUTED, lw=1.1, zorder=2)
        ax.plot([0.5, a], [i, i], color=c, lw=1.4, alpha=0.55, zorder=2, solid_capstyle="butt")
        ax.scatter(a, i, s=40, color=c, edgecolors=S.SURFACE, lw=1.2, zorder=3)
        if a >= 0.7 or a <= 0.25:
            right = a >= 0.5
            ax.text(a + (0.018 if right else -0.018), i, f"{a:.3f}", va="center",
                    ha="left" if right else "right", fontsize=8, color=S.INK,
                    fontweight="bold" if right else None)
    ax.set_yticks(range(n), t["label"])
    for tick, (_, r) in zip(ax.get_yticklabels(), t.iterrows()):
        if r["feature"] in ("hbond_RNA_F156", "pair_ipTM_protein_rna"):
            tick.set_color(S.INK)
            tick.set_fontweight("bold")
    ax.set_xlim(0.0, 1.02)
    ax.set_ylim(-0.7, n - 0.3)
    ax.set_xlabel("AUROC (replicate-averaged, C208A control excluded)")
    strip_y(ax)
    S.vgrid(ax)
    S.save(fig, "figM2_auroc_ranking")
    out = t.sort_values("AUROC_genuine", ascending=False)[
        ["label", "method", "AUROC_genuine", "AUROC_published", "feature"]]
    table("Screen AUROCs: genuine runs vs as first computed (5 deleterious vs 3 neutral)",
          out.rename(columns={"feature": "column in analysis CSVs"}))


def noise_panel(ax, runs: pd.DataFrame, col: str, sign: int, title: str, xlabel: str,
                xlim: tuple[float, float]) -> tuple[float, float]:
    """One panel of figure M3: the wild-type run set on top, genuine variant runs below."""
    used = runs[runs["used"]]
    wt = used.loc[used["wt_protein"], col].to_numpy()
    gen = used[used["genuine"] & (used["label"] != "WT")]
    order = [v for v in S.DISPLAY_ORDER if v != "WT"]
    ypos = {v: len(order) - 1 - i for i, v in enumerate(order)}
    top = len(order) + 0.6
    rng = np.random.default_rng(5)
    ax.axvspan(wt.min(), wt.max(), color=S.WASH, lw=0, zorder=0)
    ax.axvline(np.median(wt), color=S.WT, lw=0.9, alpha=0.6, zorder=1)
    ax.scatter(wt, top + rng.uniform(-0.28, 0.28, len(wt)), s=22, color=S.WT, alpha=0.7, lw=0,
               marker="s", zorder=2)
    at_zero = int((wt < 0.02).sum())
    if at_zero > 2:
        ax.annotate(f"{at_zero} runs at 0", xy=(0, top - 0.3), xytext=(0.06, top - 0.75),
                    fontsize=7.4, color=S.INK2, va="center",
                    arrowprops=dict(arrowstyle="-", color=S.MUTED, lw=0.7))
    means = gen.groupby("label")[col].mean()
    for v in order:
        cls = S.VARIANT_CLASS[v]
        vals = gen.loc[gen["label"] == v, col].to_numpy()
        ax.scatter(vals, [ypos[v]] * len(vals), s=24, color=S.CLASS_COLOR[cls], alpha=0.5,
                   lw=0, marker=S.CLASS_MARKER[cls], zorder=3)
        dot(ax, means[v], ypos[v], cls, s=60, z=4)
    labeled = means.reindex(S.DELETERIOUS + S.NEUTRAL)
    a = auroc(labeled.to_dict(), sign)
    p_noise = wt_noise_p(used, col, sign)
    ax.set_yticks([top] + [ypos[v] for v in order],
                  [f"Wild type ({len(wt)} runs)"] + order, fontsize=8)
    ax.set_ylim(-0.6, top + 0.55)
    ax.set_xlim(*xlim)
    ax.set_title(f"{title}\n\n", fontsize=9.8)
    ax.text(0, 1.012, f"AUROC {a:.3f} · exact p = {exact_p(a):.3f} · "
            f"wild-type-noise p = {p_noise:.2f}", transform=ax.transAxes, fontsize=7.8,
            color=S.MUTED, va="bottom")
    ax.set_xlabel(xlabel, fontsize=8.5)
    strip_y(ax)
    S.vgrid(ax)
    return a, p_noise


def fig_md_vs_wt_noise(t156: pd.DataFrame, stack: pd.DataFrame) -> dict:
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 5.3))
    fig.subplots_adjust(left=0.135, right=0.985, top=S.top(fig, 1.95), bottom=0.6 / 5.3,
                        wspace=0.42)
    S.header(fig, "Variant runs against the wild-type noise floor",
             "One small dot per genuine 200 ns run, large dot = variant mean. Gray squares = "
             "the wild-type runs; the wash spans their range.\nWild-type-noise p: how often "
             "wild-type runs, dealt out as fake variants with the same run counts, rank as well.")
    legend(fig, class_items(["WT", "deleterious", "neutral", "zn_ctrl"]), top_in=1.0)
    a1, p1 = noise_panel(axes[0], t156, "f156_mean", -1, "a  T156–RNA hydrogen bonds",
                         "H-bonds per frame (1 ns sampling)", (0.33, 0.75))
    a2, p2 = noise_panel(axes[1], stack, "aromatic_any", -1, "b  Aromatic stacking on a base",
                         "fraction of frames with a ring stacked", (-0.04, 1.0))
    S.save(fig, "figM3_wt_noise")
    rows = []
    for name, runs, col in [("T156 H-bonds", t156, "f156_mean"),
                            ("Aromatic stacking (any ring)", stack, "aromatic_any")]:
        used = runs[runs["used"]]
        for v in [v for v in S.DISPLAY_ORDER if v != "WT"]:
            vals = used.loc[used["genuine"] & (used["label"] == v), col]
            rows.append({"readout": name, "variant": v, "class": S.VARIANT_CLASS[v],
                         "genuine_runs": len(vals), "mean": vals.mean()})
        wt = used.loc[used["wt_protein"], col]
        rows.append({"readout": name, "variant": "wild-type set", "class": "WT",
                     "genuine_runs": len(wt), "mean": wt.mean()})
    table("Per-variant means on genuine runs, and the wild-type set", pd.DataFrame(rows))
    return {"t156": (a1, p1), "stacking": (a2, p2)}


def stacking_table(stack: pd.DataFrame) -> pd.DataFrame:
    used = stack[stack["used"]]
    gen = used[used["genuine"]].groupby("label")[STACKING].mean()
    allr = used.groupby("label")[STACKING].mean()
    t = pd.DataFrame([{"readout": c, "AUROC_genuine": auroc(gen[c].to_dict(), -1),
                       "AUROC_published": auroc(allr[c].to_dict(), -1)} for c in STACKING])
    t["exact_p_genuine"] = t["AUROC_genuine"].map(exact_p)
    t["wt_stacked_runs"] = [f"{(used.loc[used['wt_protein'], c] > 0.5).sum()}"
                            f"/{used['wt_protein'].sum()}" for c in STACKING]
    table("Aromatic stacking readouts (13): genuine runs vs as first computed", t)
    return t


# ══════════════════════════════════════════════════════════════════════════
#  Conclusions
# ══════════════════════════════════════════════════════════════════════════
def fig_n8_ceiling(af3_iptm: float, md_best: float) -> tuple[int, int]:
    p53 = tail_p(5, 3)
    ks = np.arange(8, 16)
    bonf = 0.05 / N_FEATURES_SCREENED
    ref = 0.80

    fig = plt.figure(figsize=(8.6, 4.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], left=0.085, right=0.975,
                          top=S.top(fig, 1.25), bottom=0.62 / 4.9, wspace=0.3)
    S.header(fig, "Eight labeled variants cannot carry a significant result",
             f"Exact Mann–Whitney null for 5 deleterious vs 3 neutral: 56 equally likely "
             f"rankings. Even a perfect ranking gives p = {p53[15]:.3f},\nwell above the "
             f"Bonferroni bar for the {N_FEATURES_SCREENED} readouts tested "
             f"(0.05 / {N_FEATURES_SCREENED} = {bonf:.4f}).")

    ax = fig.add_subplot(gs[0])
    ax.axhspan(1e-4, bonf, color=S.WASH, lw=0, zorder=0)
    ax.axhline(0.05, color=S.INK2, lw=0.9)
    ax.axhline(bonf, color=S.INK, lw=0.9)
    ax.text(0.508, 0.05 * 1.12, "p = 0.05", fontsize=7.6, color=S.INK2, va="bottom")
    ax.text(0.508, bonf * 1.12, f"Bonferroni, {N_FEATURES_SCREENED} readouts", fontsize=7.6,
            color=S.INK, va="bottom")
    ax.text(0.508, bonf / 1.3, "unreachable at n = 8", fontsize=7.6, color=S.MUTED, va="top")
    ax.plot(ks / 15, p53[ks], color=S.INK2, lw=1.4, zorder=2, drawstyle="steps-mid", alpha=0.35)
    ax.scatter(ks / 15, p53[ks], s=22, color=S.INK2, zorder=3, edgecolors=S.SURFACE, lw=1)
    for val, c in [(af3_iptm, S.AF3), (md_best, S.MD), (1.0, S.INK)]:
        ax.scatter(val, p53[int(round(val * 15))], s=60, color=c, zorder=4,
                   edgecolors=S.SURFACE, lw=1.3)
    ax.text(1.0 - 0.025, p53[15], f"perfect\np = {p53[15]:.3f}", ha="right", va="center",
            fontsize=7.4, color=S.INK)
    keys = [Line2D([], [], ls="", marker="o", markersize=7, markerfacecolor=c,
                   markeredgecolor=S.SURFACE,
                   label=f"{lab}: {val:.3f}, p = {p53[int(round(val * 15))]:.3f}")
            for val, lab, c in [(af3_iptm, "AF3 pair-ipTM", S.AF3),
                                (md_best, "MD, aromatic stacking", S.MD)]]
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
    ax.axhline(ref, color=S.MUTED, lw=1.2)
    ax.text(60.5, ref + 0.008, f"a solid effect, AUROC {ref:.2f}", ha="right", va="bottom",
            fontsize=7.8, color=S.INK)
    need = int(next(2 * k for k, c in zip(per_class, crit_b) if c <= ref))
    need_u = int(next(2 * k for k, c in zip(per_class, crit_u) if c <= ref))
    ax.scatter(need, ref, s=60, color=S.INK, edgecolors=S.SURFACE, lw=1.3, zorder=4)
    ax.annotate(f"≈ {need} variants ({need // 2} + {need // 2})\nto clear Bonferroni",
                xy=(need, ref), xytext=(need + 3, ref + 0.075), ha="left", va="bottom",
                fontsize=7.8, color=S.INK, linespacing=1.3,
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
    table(f"Minimum AUROC by design size (balanced; {N_FEATURES_SCREENED} readouts)",
          pd.DataFrame({"variants_total": 2 * per_class, "AUROC_p<0.05": crit_u,
                        f"AUROC_p<0.05/{N_FEATURES_SCREENED}": crit_b}),
          floatfmt=(".0f", ".3f", ".3f"))
    return need, need_u


def fig_forking_paths(af3: pd.DataFrame, t156: pd.DataFrame, screen: pd.DataFrame) -> None:
    def window(n_ns: str) -> float:  # the expanse_orig runs (genuine), resSeq 6 = T156
        h = pd.read_csv(MDA / f"hbonds_per_residue_common{n_ns}.csv", index_col=0)
        return auroc(h.loc[6].to_dict(), -1)

    used = t156[t156["used"] & t156["genuine"]]
    first_run = used[used["source"] == "expanse_orig"].groupby("label")["f156_mean"].mean()
    pooled = used.groupby("label")["f156_mean"].mean()
    in_screen = screen.set_index("feature").loc["hbond_RNA_F156", "AUROC_genuine"]
    top = af3[af3["model"] == 0].groupby("variant").mean(numeric_only=True)
    allm = af3.groupby("variant").mean(numeric_only=True)
    paths = [
        ("MD", "T156 · 1 run per variant, first 139 ns", window("139")),
        ("MD", "T156 · 1 run per variant, first 269 ns", window("269")),
        ("MD", "T156 · 1 run per variant, 200 ns, every 1 ns", auroc(first_run.to_dict(), -1)),
        ("MD", "T156 · 1–2 runs per variant, every 1 ns", auroc(pooled.to_dict(), -1)),
        ("MD", "T156 · 1–2 runs per variant, 50 frames per run", in_screen),
        ("AF3", "pair-ipTM · mean of all 40 models",
         auroc(allm["pair_ipTM_protein_rna"].to_dict(), -1)),
        ("AF3", "pair-ipTM · top-ranked model per seed",
         auroc(top["pair_ipTM_protein_rna"].to_dict(), -1)),
        ("AF3", "min pLDDT · mean of all 40 models", auroc(allm["pLDDT_min"].to_dict(), -1)),
        ("AF3", "min pLDDT · top-ranked model per seed", auroc(top["pLDDT_min"].to_dict(), -1)),
    ]
    crit = critical_auroc(5, 3, 0.05)
    n_md = sum(1 for m, *_ in paths if m == "MD")
    fig, ax = plt.subplots(figsize=(8.6, 4.9))
    fig.subplots_adjust(left=0.43, right=0.96, top=S.top(fig, 1.5), bottom=0.62 / 4.9)
    S.header(fig, "Reasonable analysis choices move the answer more than biology does",
             "Same eight labeled variants and genuine runs only. Both values that look "
             "significant rest on a fragile choice:\na 139 ns window, or the confidence of "
             "the floppy C-terminal residue.")
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
        right = val >= 0.5
        ax.text(val + (0.015 if right else -0.015), y, f"{val:.3f}", va="center",
                ha="left" if right else "right", fontsize=8.2, color=S.INK)
    ax.axhline(len(paths) - n_md - 0.5, color=S.GRID, lw=0.8)
    ax.set_yticks(ys, [lab for _, lab, _ in paths])
    strip_y(ax)
    ax.set_xlim(0.15, 1.02)
    ax.set_ylim(-0.6, len(paths) - 0.4)
    ax.set_xlabel("AUROC, 5 deleterious vs 3 neutral")
    S.vgrid(ax)
    S.save(fig, "figC2_forking_paths")
    table("AUROC under alternative analysis choices (genuine runs)",
          pd.DataFrame(paths, columns=["method", "analysis", "AUROC"]))


def main() -> None:
    S.apply()
    af3 = load_af3()
    t156 = load_t156_runs()
    stack = load_stacking_runs()
    used = t156[t156["used"]]
    assert used["wt_protein"].sum() == 19 and (used["genuine"] & (used["label"] != "WT")).sum() == 14
    af3_iptm = auroc(af3.groupby("variant")["pair_ipTM_protein_rna"].mean().to_dict(), -1)
    assert abs(af3_iptm - 11 / 15) < 1e-9

    fig_structure()
    fig_af3_confidence(af3)
    fig_af3_rna_control(af3)
    fig_md_sampling(t156)
    screen = screen_table()
    fig_auroc_ranking(screen)
    noise = fig_md_vs_wt_noise(t156, stack)
    st = stacking_table(stack)
    md_best = st["AUROC_genuine"].max()
    assert abs(noise["t156"][0] - 4 / 15) < 1e-9 and abs(md_best - 14 / 15) < 1e-9
    need, need_u = fig_n8_ceiling(af3_iptm, md_best)
    fig_forking_paths(af3, t156, screen)
    (FIGDATA / "tables.md").write_text(
        "# Numbers behind the figures\n\nGenerated by figures/src/make_figures.py.\n\n"
        + "\n".join(TABLES))
    print(f"AF3 pair-ipTM {af3_iptm:.3f}; T156 genuine {noise['t156'][0]:.3f} "
          f"(WT-noise p {noise['t156'][1]:.2f}); stacking genuine {noise['stacking'][0]:.3f} "
          f"(WT-noise p {noise['stacking'][1]:.3f}); best screen MD "
          f"{screen[screen['method'] == 'MD']['AUROC_genuine'].max():.3f}; variants needed for "
          f"AUROC 0.80: {need_u} uncorrected, {need} after Bonferroni({N_FEATURES_SCREENED})")


if __name__ == "__main__":
    main()
