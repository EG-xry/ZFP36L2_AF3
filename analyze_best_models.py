#!/usr/bin/env python3
"""
ZFP36L2 AF3 — variant comparison: model-0 only (n=8) AND all models (n=40).

For each variant, computes key interface / confidence metrics and reports
mean ± SD.  Two parallel analyses are produced:

  1. model_0 only (top-ranked AF3 model per run, 8 per variant)
  2. all models   (all 5 models × 8 runs = 40 per variant)

Statistical comparisons use Welch's two-sample t-test (unequal variance).

Outputs (all under results/):
    variant_avg.csv / .md                model-0 averaged metrics
    variant_vs_wt.csv / .md              model-0 t-test table
    plots/pair_iptm_barplot.png  …       model-0 bar plots

    variant_avg_all.csv / .md            all-models averaged metrics
    variant_vs_wt_all.csv / .md          all-models t-test table
    plots/pair_iptm_barplot_all.png …    all-models bar plots
    plots/composite_deviation_with_plddt*.png
                                        alternate composite plots including pLDDT
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
from scipy import stats

# Import core analysis functions from the main module
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_complexes as ac  # noqa: E402

BEST_MODEL_INDEX = 0  # model_0 is the top-ranked model per AF3 run

# ── Column definitions ────────────────────────────────────────────────────
GROUP_COLS = ["Protein", "RNA", "ProteinExpClass", "RNAExpClass"]

OUTPUT_COLS = [
    "Protein", "RNA", "ProteinExpClass", "RNAExpClass",
    "pair_ipTM_protein_rna",
    "pair_pae_min_cross_mean",
    "pLDDT", "pLDDT_ZF1", "pLDDT_ZF2",
    "pLDDT_min",
    "RMSD_whole",
]
METRIC_COLS = [c for c in OUTPUT_COLS if c not in GROUP_COLS]

# Metrics for statistical testing (no RMSD — see note in methodology)
TEST_METRICS = [
    ("pair_ipTM_protein_rna", "pair ipTM (protein–RNA)", True),
    ("pair_pae_min_cross_mean", "PAE cross-mean (Å)",     False),
    ("pLDDT",                   "pLDDT",                   True),
    ("pLDDT_min",               "pLDDT min (Cα)",          True),
]

# Protein order: WT → deleterious → neutral
PROTEIN_ORDER = [
    "WT", "C206F", "C208A", "G204D", "R160W", "R160E", "R184C", "Y154H", "H212A",
    "D210E", "S185N", "T187S",
]

CLASS_COLORS = {"WT": "#89CFF0", "neutral": "#0D3B66", "deleterious": "#F4722B"}
EXCLUDE_FROM_FIGURES = {"H212A"}
ASSUMED_N_MODEL0 = 8
ASSUMED_N_ALL = 40


# ── Reusable helpers ──────────────────────────────────────────────────────
def _color_for(protein: str, class_map: dict) -> str:
    return CLASS_COLORS.get(class_map.get(protein, ""), "#888888")


def derive_protein_order(df_per_model: pd.DataFrame) -> list[str]:
    present = set(df_per_model["Protein"].dropna().astype(str))
    ordered = [p for p in PROTEIN_ORDER if p in present]
    extras = sorted(p for p in present if p not in PROTEIN_ORDER)
    return ordered + extras


def filter_for_figures(df_per_model: pd.DataFrame) -> pd.DataFrame:
    return df_per_model[~df_per_model["Protein"].isin(EXCLUDE_FROM_FIGURES)].copy()


def _assume_n(values: np.ndarray, target_n: int, seed: int) -> np.ndarray:
    """
    Enforce a fixed sample size for significance testing.
    - If len(values) < target_n: bootstrap upsample with replacement.
    - If len(values) > target_n: downsample without replacement.
    """
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return arr
    if len(arr) == target_n:
        return arr
    rng = np.random.default_rng(seed)
    if len(arr) < target_n:
        idx = rng.choice(len(arr), size=target_n, replace=True)
        return arr[idx]
    idx = rng.choice(len(arr), size=target_n, replace=False)
    return arr[idx]


def make_barplot(
    df_per_model: pd.DataFrame,
    metric_col: str,
    ylabel: str,
    higher_is_better: bool,
    protein_order: list,
    class_map: dict,
    sig_proteins: set,
    out_path: Path,
    subtitle: str = "",
) -> None:
    agg = (
        df_per_model.groupby("Protein")[metric_col]
        .agg(["mean", "std"])
        .reindex(protein_order)
    )
    colors = [_color_for(p, class_map) for p in protein_order]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = np.arange(len(protein_order))
    ax.bar(x, agg["mean"], yerr=agg["std"], capsize=4,
           color=colors, edgecolor="black", linewidth=0.6, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels(protein_order, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=13)
    ax.tick_params(axis="y", labelsize=12)
    direction = "↑ higher = better" if higher_is_better else "↓ lower = better"
    title = f"{ylabel}  ({direction})"
    if subtitle:
        title += f"\n{subtitle}"
    ax.set_title(title, fontsize=13)
    ax.grid(axis="y", alpha=0.3, zorder=0)

    # Significance stars at uniform height (top gridline)
    if sig_proteins:
        y_top = ax.get_ylim()[1]
        for i, protein in enumerate(protein_order):
            if protein in sig_proteins:
                ax.text(i, y_top, "*", ha="center", va="top",
                        fontsize=16, fontweight="bold", color="black")

    legend_elements = [
        Patch(facecolor=CLASS_COLORS["WT"], edgecolor="black", label="WT (light blue)"),
        Patch(facecolor=CLASS_COLORS["deleterious"], edgecolor="black", label="Deleterious (orange)"),
        Patch(facecolor=CLASS_COLORS["neutral"], edgecolor="black", label="Neutral (dark blue)"),
    ]
    ax.legend(handles=legend_elements, loc="best", fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


# ── Composite "variant health" master graph ─────────────────────────────
# Metrics used in the composite: (column, display_label, higher_is_better)
COMPOSITE_METRICS_BASE = [
    ("pair_ipTM_protein_rna", "ipTM\n(protein–RNA)", True),
    ("pair_pae_min_cross_mean", "PAE\n(cross-mean)",  False),
]
COMPOSITE_METRICS_WITH_PLDDT = COMPOSITE_METRICS_BASE + [
    ("pLDDT", "pLDDT", True),
]


def make_composite_graph(
    df_per_model: pd.DataFrame,
    protein_order: list,
    class_map: dict,
    out_path: Path,
    composite_metrics: list[tuple[str, str, bool]],
    subtitle: str = "",
) -> None:
    """
    Master figure combining ipTM and PAE into one view.

    Left panel : heatmap of z-scores (variant vs WT) for the selected metrics.
                 Red = worse than WT, blue = better than WT.
    Right panel: composite "deviation score" — mean z-score across metrics,
                 ranking variants by overall departure from WT.
    """
    variants = [p for p in protein_order if p != "WT"]
    wt_data = df_per_model[df_per_model["Protein"] == "WT"]

    # Compute z-scores: (variant_mean − WT_mean) / WT_sd
    # Flip sign so that POSITIVE always means WORSE
    z_matrix = np.zeros((len(variants), len(composite_metrics)))
    sig_matrix = np.full((len(variants), len(composite_metrics)), False)  # significance flags
    # Per-model z-scores for computing composite SE
    z_per_model: dict[str, list[np.ndarray]] = {p: [] for p in variants}
    labels_row = []
    labels_col = []

    for j, (col, label, higher_better) in enumerate(composite_metrics):
        wt_vals = wt_data[col].dropna().values
        wt_mean = np.mean(wt_vals)
        wt_std = np.std(wt_vals, ddof=1)
        if wt_std == 0:
            wt_std = 1e-9  # avoid division by zero
        labels_col.append(label)
        for i, protein in enumerate(variants):
            var_vals = df_per_model[df_per_model["Protein"] == protein][col].dropna().values
            var_mean = np.mean(var_vals)
            z = (var_mean - wt_mean) / wt_std
            if higher_better:
                z = -z
            z_matrix[i, j] = z
            # Welch's t-test for significance
            if len(wt_vals) >= 2 and len(var_vals) >= 2:
                _, p_val = stats.ttest_ind(var_vals, wt_vals, equal_var=False)
                sig_matrix[i, j] = p_val < 0.05
            # Per-model z-scores (one per model)
            z_ind = (var_vals - wt_mean) / wt_std
            if higher_better:
                z_ind = -z_ind
            z_per_model[protein].append(z_ind)
        if j == 0:
            labels_row = variants

    # Composite score: mean of z-scores (positive = worse overall)
    composite = np.mean(z_matrix, axis=1)

    # Composite SE: for each model, average its z-scores across metrics,
    # then compute SD / sqrt(n) of those per-model composites
    composite_se = np.zeros(len(variants))
    for i, protein in enumerate(variants):
        # z_per_model[protein] is a list of n_metrics arrays, each length n_models.
        # Stack into (n_metrics, n_models) and average over metrics -> per-model composite.
        stacked = np.vstack(z_per_model[protein])  # (n_metrics, n_models)
        per_model_composite = np.mean(stacked, axis=0)  # (n_models,)
        composite_se[i] = np.std(per_model_composite, ddof=1) / np.sqrt(len(per_model_composite))

    # ── Build figure ────────────────────────────────────────────────────
    fig = plt.figure(figsize=(11, 5.5))
    gs = GridSpec(1, 3, width_ratios=[3, 0.15, 1.8], wspace=0.05)

    # --- Left: heatmap ---
    ax_heat = fig.add_subplot(gs[0])
    vmax = max(abs(z_matrix.min()), abs(z_matrix.max()), 1.0)
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "bwr_custom", ["#0D3B66", "#FFFFFF", "#F4722B"], N=256,
    )
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax_heat.imshow(z_matrix, cmap=cmap, norm=norm, aspect="auto")

    ax_heat.set_xticks(np.arange(len(labels_col)))
    ax_heat.set_xticklabels(labels_col, fontsize=11, ha="center")
    ax_heat.set_yticks(np.arange(len(labels_row)))
    ax_heat.set_yticklabels(labels_row, fontsize=12)
    ax_heat.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Remove internal grid lines but keep the outer frame
    ax_heat.set_xticks([x - 0.5 for x in range(1, len(labels_col))], minor=True)
    ax_heat.set_yticks([y - 0.5 for y in range(1, len(labels_row))], minor=True)
    ax_heat.tick_params(which="minor", length=0)  # hide minor tick marks
    ax_heat.grid(False)  # no grid at all on heatmap

    # Annotate cells with z-score values (star if Welch t-test p < 0.05)
    for i in range(len(labels_row)):
        for j in range(len(labels_col)):
            val = z_matrix[i, j]
            txt_color = "white" if abs(val) > vmax * 0.6 else "black"
            star = " *" if sig_matrix[i, j] else ""
            ax_heat.text(j, i, f"{val:+.2f}{star}", ha="center", va="center",
                         fontsize=10, fontweight="bold", color=txt_color)

    # Colorbar
    ax_cb = fig.add_subplot(gs[1])
    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("z-score (+ = worse than WT)", fontsize=10)
    cb.ax.tick_params(labelsize=10)

    # --- Right: composite deviation bar ---
    ax_bar = fig.add_subplot(gs[2])
    y = np.arange(len(variants))
    bar_colors = [_color_for(p, class_map) for p in variants]
    ax_bar.barh(y, composite, xerr=composite_se, capsize=4,
                color=bar_colors, edgecolor="black", linewidth=0.6,
                error_kw={"linewidth": 1.2, "color": "black"})
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels([""] * len(variants))  # labels already on heatmap
    ax_bar.set_xlabel("Composite deviation\n(mean z-score, + = worse)", fontsize=10)
    ax_bar.invert_yaxis()
    ax_bar.tick_params(axis="x", labelsize=11)
    ax_bar.axvline(0, color="black", linewidth=0.8, linestyle="--", alpha=0.5)
    ax_bar.grid(axis="x", alpha=0.3)

    # Title — place high enough to clear column headers
    title = "Variant Deviation from WT — Composite View"
    if subtitle:
        title += f"\n{subtitle}"
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.08)

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def welch_ttest_table(
    df_per_model: pd.DataFrame,
    protein_order: list,
    class_map: dict,
    assumed_n: int,
    analysis_tag: str,
) -> pd.DataFrame:
    wt = df_per_model[df_per_model["Protein"] == "WT"]
    rows = []
    for protein in protein_order:
        if protein == "WT":
            continue
        variant = df_per_model[df_per_model["Protein"] == protein]
        exp_class = class_map.get(protein, "")
        for col, label, higher_better in TEST_METRICS:
            wt_raw = wt[col].dropna().values
            var_raw = variant[col].dropna().values
            if len(wt_raw) < 2 or len(var_raw) < 2:
                continue
            wt_seed = abs(hash((analysis_tag, "WT", col, assumed_n))) % (2**32)
            var_seed = abs(hash((analysis_tag, protein, col, assumed_n))) % (2**32)
            wt_vals = _assume_n(wt_raw, assumed_n, wt_seed)
            var_vals = _assume_n(var_raw, assumed_n, var_seed)
            if len(wt_vals) < 2 or len(var_vals) < 2:
                continue
            t_stat, p_val = stats.ttest_ind(var_vals, wt_vals, equal_var=False)
            wt_mean = float(np.mean(wt_vals))
            var_mean = float(np.mean(var_vals))
            sig = p_val < 0.05
            if var_mean > wt_mean:
                direction = "higher"
            elif var_mean < wt_mean:
                direction = "lower"
            else:
                direction = "equal"
            if higher_better:
                better = var_mean > wt_mean
            else:
                better = var_mean < wt_mean
            rows.append({
                "Protein": protein,
                "ProteinExpClass": exp_class,
                "Metric": label,
                "WT_mean": round(wt_mean, 4),
                "Variant_mean": round(var_mean, 4),
                "Difference": round(var_mean - wt_mean, 4),
                "t_statistic": round(t_stat, 4),
                "p_value": p_val,
                "Significant (p<0.05)": "Yes" if sig else "No",
                "Direction vs WT": direction,
                "Better than WT": ("Yes *" if better else "No") if sig else ("-"),
                "n_assumed_each_group": assumed_n,
            })
    return pd.DataFrame(rows)


def run_analysis(
    df_all: pd.DataFrame,
    output_dir: Path,
    suffix: str,
    description: str,
) -> None:
    """
    Produce averaged table, significance table, and bar plots for a given
    set of per-model rows.  ``suffix`` is appended to every output filename
    (empty string for model-0, '_all' for all models).
    """
    tag = f" ({description})" if description else ""
    sfx = f"_{suffix}" if suffix else ""
    assumed_n = ASSUMED_N_ALL if suffix == "all" else ASSUMED_N_MODEL0

    class_map = dict(zip(df_all["Protein"], df_all["ProteinExpClass"]))
    protein_order = derive_protein_order(df_all)
    n_by_protein = (
        df_all.groupby("Protein").size().reindex(protein_order).dropna().astype(int)
    )
    if len(n_by_protein):
        n_summary = ", ".join(f"{protein}={n}" for protein, n in n_by_protein.items())
    else:
        n_summary = "n/a"

    # ── 1. Averaged table ─────────────────────────────────────────────────
    df_avg = df_all.groupby(GROUP_COLS, dropna=False)[METRIC_COLS].mean().reset_index()
    df_sd = df_all.groupby(GROUP_COLS, dropna=False)[METRIC_COLS].std(ddof=1).add_suffix("_sd").reset_index()
    df_out = df_avg.merge(df_sd, on=GROUP_COLS, how="left")
    n_runs = df_all.groupby(GROUP_COLS).size().rename("n_models").reset_index()
    df_out = df_out.merge(n_runs, on=GROUP_COLS, how="left")
    df_out.sort_values(["ProteinExpClass", "Protein"], inplace=True)

    avg_csv = output_dir / f"variant_avg{sfx}.csv"
    avg_md = output_dir / f"variant_avg{sfx}.md"
    df_out.to_csv(avg_csv, index=False)
    avg_md.write_text("\n".join([
        f"# ZFP36L2 AF3 — Averaged metrics{tag}",
        "",
        f"- Variants analysed: {len(df_out)}",
        f"- Models per variant: {n_summary}",
        f"- {description}",
        "- `_sd` columns are standard deviations (ddof = 1)",
        "",
        "## Metric definitions",
        "",
        "| Metric | Description |",
        "|--------|-------------|",
        "| pair_ipTM_protein_rna | Pairwise protein↔RNA ipTM (interface confidence; higher = better) |",
        "| pair_pae_min_cross_mean | Mean of PAE(protein→RNA) and PAE(RNA→protein); lower = better |",
        "| pLDDT | Mean per-residue Cα confidence (0–100) |",
        "| pLDDT_ZF1 / ZF2 | Mean Cα pLDDT for zinc finger 1 / 2 |",
        "| pLDDT_min | Minimum per-residue Cα pLDDT across the protein domain (matches AF3 server coloring) |",
        "| RMSD_whole | Backbone RMSD vs WT reference (Å) |",
        "",
        df_out.to_markdown(index=False),
    ]) + "\n")
    print(f"Wrote {avg_csv}")
    print(f"Wrote {avg_md}")

    # ── 2. Welch t-test significance table ────────────────────────────────
    df_sig = welch_ttest_table(
        df_all,
        protein_order,
        class_map,
        assumed_n=assumed_n,
        analysis_tag=suffix or "model0",
    )

    sig_csv = output_dir / f"variant_vs_wt{sfx}.csv"
    sig_md = output_dir / f"variant_vs_wt{sfx}.md"
    df_sig.to_csv(sig_csv, index=False)
    sig_md.write_text("\n".join([
        f"# ZFP36L2 AF3 — Welch's t-test: variant vs WT{tag}",
        "",
        "## Methodology",
        "",
        "For each metric, an independent two-sample Welch's t-test (unequal",
        "variance) compares each variant's available values against WT's",
        "available values for the same analysis set. Significance: α = 0.05.",
        f"For hypothesis-testing only, each group is normalized to assumed n={assumed_n}",
        "(bootstrap upsampling when fewer values exist; deterministic downsampling when more exist).",
        "",
        "## Metric interpretation",
        "",
        "| Metric | Higher is … | Meaning |",
        "|--------|-------------|---------|",
        "| pair ipTM (protein–RNA) | **Better** | Pairwise protein↔RNA predicted TM-score; higher = more confident interface |",
        "| PAE cross-mean | **Worse** | Predicted aligned error between protein ↔ RNA (Å); higher = less confident |",
        "| pLDDT | **Better** | Mean per-residue Cα confidence (0–100); higher = more confident structure |",
        "| pLDDT min (Cα) | **Better** | Lowest per-residue Cα pLDDT in the protein domain; higher = fewer disordered regions |",
        "",
        "RMSD is excluded from this test because WT RMSD is ~0 by definition",
        "(self-aligned reference), making the comparison trivially significant",
        "for every variant.",
        "",
        "## Results",
        "",
        '- **"Better than WT = Yes *"** means the variant is *significantly* better than WT for that metric (p < 0.05).',
        '- **"-"** means the difference is not statistically significant.',
        "",
        df_sig.to_markdown(index=False),
    ]) + "\n")
    print(f"Wrote {sig_csv}")
    print(f"Wrote {sig_md}")

    # ── 3. Bar plots ──────────────────────────────────────────────────────
    df_fig = filter_for_figures(df_all)
    protein_order_fig = derive_protein_order(df_fig)
    description_base = description.split(";")[0].strip() if description else ""
    subtitle_fig = description_base

    sig_lookup: dict[str, set] = {}
    for _, row in df_sig.iterrows():
        if row["Significant (p<0.05)"] == "Yes":
            sig_lookup.setdefault(row["Metric"], set()).add(row["Protein"])

    for col, label, higher_better in TEST_METRICS:
        safe_name = col.lower().replace("pair_pae_min_cross_mean", "pae_cross_mean")
        fname = f"{safe_name}_barplot{sfx}.png"
        make_barplot(
            df_fig, col, label, higher_better, protein_order_fig,
            class_map, sig_lookup.get(label, set()),
            output_dir / "plots" / fname, subtitle=subtitle_fig,
        )
        print(f"Wrote {output_dir / 'plots' / fname}")

    # RMSD bar plot (no stars — WT is self-aligned, so comparison is meaningless)
    rmsd_fname = f"rmsd_whole_barplot{sfx}.png"
    make_barplot(
        df_fig, "RMSD_whole", "RMSD (Å)", False, protein_order_fig,
        class_map, set(), output_dir / "plots" / rmsd_fname,
        subtitle=subtitle_fig,
    )
    print(f"Wrote {output_dir / 'plots' / rmsd_fname}")

    # ── 4. Composite master graph (ipTM + PAE) ──────────────────────────
    composite_fname = f"composite_deviation{sfx}.png"
    make_composite_graph(
        df_fig, protein_order_fig, class_map,
        output_dir / "plots" / composite_fname,
        COMPOSITE_METRICS_BASE,
        subtitle=subtitle_fig,
    )
    print(f"Wrote {output_dir / 'plots' / composite_fname}")

    composite_plddt_fname = f"composite_deviation_with_plddt{sfx}.png"
    make_composite_graph(
        df_fig, protein_order_fig, class_map,
        output_dir / "plots" / composite_plddt_fname,
        COMPOSITE_METRICS_WITH_PLDDT,
        subtitle=subtitle_fig,
    )
    print(f"Wrote {output_dir / 'plots' / composite_plddt_fname}")


# ── Main ──────────────────────────────────────────────────────────────────
def main() -> None:
    input_root = ac.DEFAULT_DATA_DIR
    output_dir = ac.BASE_DIR / "results"
    ac.init_output_dirs(output_dir)

    # ── Discover all structures ───────────────────────────────────────────
    all_complexes = ac.list_complexes(input_root)

    # ── Reference structures & confidence data ────────────────────────────
    refs = ac.build_reference_map(all_complexes)
    any_ref_path = next(iter(next(iter(refs.values())).values()))
    canonical_seq = ac.load_canonical_sequence(any_ref_path)
    confidence_metrics = ac.load_confidence_metrics(input_root)
    full_data_plddts = ac.load_full_data_plddts(input_root)
    ref_cache: dict = {}

    # ── Analyse ALL structures once ───────────────────────────────────────
    print(f"Analysing {len(all_complexes)} structures …")
    all_records = []
    for meta in all_complexes:
        ref_by_model = refs[meta.reference_label]
        ref_path = (
            ref_by_model.get(meta.model_index)
            or ref_by_model.get(0)
            or next(iter(ref_by_model.values()))
        )
        cache_key = f"{meta.reference_label}:{Path(ref_path).name}"
        if cache_key not in ref_cache:
            ref_cache[cache_key] = ac.mda.Universe(ref_path)
        record = ac.analyze_complex(
            meta, ref_cache[cache_key], canonical_seq,
            confidence_metrics, full_data_plddts,
        )
        all_records.append(record)

    df_everything = pd.DataFrame(all_records)

    # ── Split into model-0 only vs all ────────────────────────────────────
    df_model0 = df_everything[
        df_everything["Model"] % ac.MODELS_PER_RUN == BEST_MODEL_INDEX
    ].copy()

    n_model0_by_protein = (
        df_model0.groupby("Protein").size().reindex(derive_protein_order(df_model0)).dropna().astype(int)
    )
    n_all_by_protein = (
        df_everything.groupby("Protein").size().reindex(derive_protein_order(df_everything)).dropna().astype(int)
    )
    n_model0_summary = ", ".join(f"{p}={n}" for p, n in n_model0_by_protein.items())
    n_all_summary = ", ".join(f"{p}={n}" for p, n in n_all_by_protein.items())

    # ── Run both analyses ─────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print("  Analysis 1: model_0 only")
    print(f"  Counts by variant: {n_model0_summary}")
    print(f"{'='*60}")
    run_analysis(df_model0, output_dir, "",
                 f"model_0 only; counts by variant: {n_model0_summary}")

    print(f"\n{'='*60}")
    print("  Analysis 2: all models")
    print(f"  Counts by variant: {n_all_summary}")
    print(f"{'='*60}")
    run_analysis(df_everything, output_dir, "all",
                 f"all models; counts by variant: {n_all_summary}")

    print("\nDone.")


if __name__ == "__main__":
    main()
