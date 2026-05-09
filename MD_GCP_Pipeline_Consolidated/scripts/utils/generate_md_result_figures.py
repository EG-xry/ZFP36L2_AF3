#!/usr/bin/env python3
"""Generate focused MD result figures and summary stats from compiled bundle tables.

Outputs:
- figures/figure06_md_variant_delta_vs_wt.png
- figures/figure07_md_class_delta_summary.png
- tables/md_neutral_vs_deleterious_vs_wt_stats.csv
- MD_RESULTS_SUMMARY.md
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


CLASS_COLORS = {
    "WT": "#2E7D32",
    "neutral": "#1565C0",
    "deleterious": "#C62828",
}


def _ci95(series: pd.Series) -> tuple[float, float]:
    vals = series.dropna().astype(float).to_numpy()
    n = len(vals)
    if n == 0:
        return np.nan, np.nan
    mean = float(np.mean(vals))
    if n == 1:
        return mean, mean
    sem = stats.sem(vals)
    ci = float(stats.t.ppf(0.975, df=n - 1) * sem)
    return mean - ci, mean + ci


def _mean_sd_n(series: pd.Series) -> tuple[float, float, int]:
    vals = series.dropna().astype(float).to_numpy()
    n = len(vals)
    if n == 0:
        return np.nan, np.nan, 0
    mean = float(np.mean(vals))
    sd = float(np.std(vals, ddof=1)) if n > 1 else 0.0
    return mean, sd, n


def _variant_order(df: pd.DataFrame) -> list[str]:
    wt = ["WT"] if "WT" in set(df["variant"]) else []
    neutral = sorted(df.loc[df["classification"] == "neutral", "variant"].dropna().astype(str).unique())
    deleterious = sorted(df.loc[df["classification"] == "deleterious", "variant"].dropna().astype(str).unique())
    extras = sorted(set(df["variant"].dropna().astype(str)) - set(wt) - set(neutral) - set(deleterious))
    return wt + neutral + deleterious + extras


def make_variant_delta_figure(df: pd.DataFrame, out_path: Path) -> None:
    order = _variant_order(df)
    d = df.set_index("variant").reindex(order).reset_index()
    x = np.arange(len(d))
    colors = [CLASS_COLORS.get(c, "#666666") for c in d["classification"]]
    markers = ["^" if s == "fallback_short_rerun" else "o" for s in d["run_source_group"]]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    metrics = [
        ("delta_ca_vs_wt", "Protein drift delta vs WT (CA RMSD A)"),
        ("delta_rna_vs_wt", "RNA drift delta vs WT (P RMSD A)"),
    ]

    for ax, (col, ylabel) in zip(axes, metrics):
        vals = d[col].astype(float).to_numpy()
        ax.bar(x, vals, color=colors, edgecolor="black", linewidth=0.6, alpha=0.75)
        for xi, yi, marker in zip(x, vals, markers):
            ax.scatter(xi, yi, marker=marker, s=55, color="black", zorder=3)
        ax.axhline(0.0, color="black", linestyle="--", linewidth=0.9, alpha=0.8)
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=0.3)

    axes[1].set_xticks(x)
    axes[1].set_xticklabels(d["variant"], rotation=30, ha="right")
    axes[0].set_title("MD drift deltas vs WT baseline (mixed dataset: 8 main + 2 fallback)")

    legend_handles = [
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=CLASS_COLORS["WT"], markeredgecolor="black", markersize=10, label="WT"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=CLASS_COLORS["neutral"], markeredgecolor="black", markersize=10, label="Neutral"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor=CLASS_COLORS["deleterious"], markeredgecolor="black", markersize=10, label="Deleterious"),
        plt.Line2D([0], [0], marker="o", color="black", linestyle="None", markersize=7, label="Main batch"),
        plt.Line2D([0], [0], marker="^", color="black", linestyle="None", markersize=7, label="Fallback short rerun"),
    ]
    axes[0].legend(handles=legend_handles, loc="best", fontsize=9, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def make_class_summary_figure(df: pd.DataFrame, out_path: Path, pvals: dict[str, float]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    class_order = ["deleterious", "neutral"]
    x = np.arange(len(class_order))
    rng = np.random.default_rng(42)

    panels = [
        ("delta_ca_vs_wt", "Protein delta vs WT (CA RMSD A)"),
        ("delta_rna_vs_wt", "RNA delta vs WT (P RMSD A)"),
    ]

    for ax, (col, title) in zip(axes, panels):
        for i, cls in enumerate(class_order):
            sub = df[df["classification"] == cls][col].astype(float).dropna().to_numpy()
            jitter = rng.uniform(-0.08, 0.08, size=len(sub))
            ax.scatter(
                np.full(len(sub), i) + jitter,
                sub,
                s=60,
                color=CLASS_COLORS[cls],
                edgecolor="black",
                linewidth=0.5,
                alpha=0.9,
                zorder=3,
            )
            if len(sub):
                mean = float(np.mean(sub))
                lo, hi = _ci95(pd.Series(sub))
                ax.errorbar(i, mean, yerr=[[mean - lo], [hi - mean]], fmt="o", color="black", capsize=5, zorder=4)
                ax.text(i, mean, f"  mean={mean:+.2f}", fontsize=8, va="center")

        ax.axhline(0.0, color="black", linestyle="--", linewidth=0.9, alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"Deleterious\n(n={len(df[df['classification']=='deleterious'])})", f"Neutral\n(n={len(df[df['classification']=='neutral'])})"])
        ax.set_title(title)
        ax.grid(axis="y", alpha=0.3)

    axes[0].set_ylabel("Delta vs WT baseline")
    fig.suptitle(
        "Class-level MD delta summary\n"
        f"Welch p (CA)={pvals['welch_ca']:.3f}, Welch p (RNA)={pvals['welch_rna']:.3f}",
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_stats_table(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cls in ["WT", "deleterious", "neutral"]:
        sub = df[df["classification"] == cls]
        for col in [
            "ca_rmsd_a_final_vs_start",
            "p_rmsd_b_final_vs_start",
            "delta_ca_vs_wt",
            "delta_rna_vs_wt",
        ]:
            mean, sd, n = _mean_sd_n(sub[col])
            rows.append(
                {
                    "group": cls,
                    "metric": col,
                    "n": n,
                    "mean": mean,
                    "sd": sd,
                }
            )

    # Between-class tests (deleterious vs neutral)
    d = df[df["classification"] == "deleterious"]
    n = df[df["classification"] == "neutral"]
    for col in ["delta_ca_vs_wt", "delta_rna_vs_wt"]:
        dvals = d[col].astype(float).dropna().to_numpy()
        nvals = n[col].astype(float).dropna().to_numpy()
        t_stat, welch_p = stats.ttest_ind(dvals, nvals, equal_var=False)
        u_stat, mwu_p = stats.mannwhitneyu(dvals, nvals, alternative="two-sided")
        rows.append(
            {
                "group": "deleterious_vs_neutral",
                "metric": col,
                "n": len(dvals) + len(nvals),
                "mean": np.nan,
                "sd": np.nan,
                "test": "welch_t",
                "statistic": t_stat,
                "p_value": welch_p,
            }
        )
        rows.append(
            {
                "group": "deleterious_vs_neutral",
                "metric": col,
                "n": len(dvals) + len(nvals),
                "mean": np.nan,
                "sd": np.nan,
                "test": "mann_whitney_u",
                "statistic": u_stat,
                "p_value": mwu_p,
            }
        )

    # One-sample tests vs WT baseline (delta == 0)
    for cls in ["deleterious", "neutral"]:
        sub = df[df["classification"] == cls]
        for col in ["delta_ca_vs_wt", "delta_rna_vs_wt"]:
            vals = sub[col].astype(float).dropna().to_numpy()
            t_stat, p_val = stats.ttest_1samp(vals, 0.0)
            rows.append(
                {
                    "group": f"{cls}_vs_wt_baseline",
                    "metric": col,
                    "n": len(vals),
                    "mean": np.nan,
                    "sd": np.nan,
                    "test": "one_sample_t_vs_zero",
                    "statistic": t_stat,
                    "p_value": p_val,
                }
            )

    out = pd.DataFrame(rows)
    cols = ["group", "metric", "n", "mean", "sd", "test", "statistic", "p_value"]
    for c in cols:
        if c not in out.columns:
            out[c] = np.nan
    return out[cols]


def write_summary_md(path: Path, df: pd.DataFrame, stats_table: pd.DataFrame) -> None:
    wt = df[df["classification"] == "WT"].iloc[0]
    d = df[df["classification"] == "deleterious"]
    n = df[df["classification"] == "neutral"]

    d_ca = d["delta_ca_vs_wt"].mean()
    n_ca = n["delta_ca_vs_wt"].mean()
    d_rna = d["delta_rna_vs_wt"].mean()
    n_rna = n["delta_rna_vs_wt"].mean()

    welch_ca = stats_table[
        (stats_table["group"] == "deleterious_vs_neutral")
        & (stats_table["metric"] == "delta_ca_vs_wt")
        & (stats_table["test"] == "welch_t")
    ]["p_value"].iloc[0]
    welch_rna = stats_table[
        (stats_table["group"] == "deleterious_vs_neutral")
        & (stats_table["metric"] == "delta_rna_vs_wt")
        & (stats_table["test"] == "welch_t")
    ]["p_value"].iloc[0]

    lines = [
        "# MD Results Summary (Mixed 8+2 Dataset)",
        "",
        "## Data used",
        "",
        "- 8 variants from main cloud batch output.",
        "- 2 variants from fallback short rerun: S185N and C208A.",
        "- WT baseline is the WT run in the same synced cloud output set.",
        "",
        "## WT baseline",
        "",
        f"- WT protein drift (CA RMSD): `{wt['ca_rmsd_a_final_vs_start']:.4f} A`",
        f"- WT RNA drift (P RMSD): `{wt['p_rmsd_b_final_vs_start']:.4f} A`",
        "",
        "## Mean deltas vs WT (negative = lower drift than WT)",
        "",
        f"- Deleterious, protein delta: `{d_ca:+.4f} A`",
        f"- Neutral, protein delta: `{n_ca:+.4f} A`",
        f"- Deleterious, RNA delta: `{d_rna:+.4f} A`",
        f"- Neutral, RNA delta: `{n_rna:+.4f} A`",
        "",
        "## Neutral vs Deleterious class comparison",
        "",
        f"- Welch t-test (protein delta): `p={welch_ca:.4f}`",
        f"- Welch t-test (RNA delta): `p={welch_rna:.4f}`",
        "",
        "Interpretation: with this small mixed dataset, there is no statistically strong class separation between neutral and deleterious in MD drift deltas.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate focused MD result figures from compiled tables.")
    parser.add_argument(
        "--bundle-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "presentation" / "compiled_20260429_mixed8plus2",
        help="Compiled bundle directory containing tables/ and figures/.",
    )
    args = parser.parse_args()

    bundle_dir = args.bundle_dir
    tables_dir = bundle_dir / "tables"
    figures_dir = bundle_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    combined = pd.read_csv(tables_dir / "combined_af3_md_selected.csv")
    d = combined[combined["qc_status"] == "complete_minimal_md"].copy()
    d = d[d["classification"].isin(["WT", "neutral", "deleterious"])].copy()

    wt_row = d[d["classification"] == "WT"]
    if wt_row.empty:
        raise RuntimeError("WT row not found in combined table.")
    wt = wt_row.iloc[0]
    d["delta_ca_vs_wt"] = d["ca_rmsd_a_final_vs_start"].astype(float) - float(wt["ca_rmsd_a_final_vs_start"])
    d["delta_rna_vs_wt"] = d["p_rmsd_b_final_vs_start"].astype(float) - float(wt["p_rmsd_b_final_vs_start"])
    if "run_source_group" not in d.columns:
        d["run_source_group"] = "unknown"

    d.to_csv(tables_dir / "md_variant_deltas_vs_wt.csv", index=False)

    # Stats table
    stats_table = build_stats_table(d)
    stats_csv = tables_dir / "md_neutral_vs_deleterious_vs_wt_stats.csv"
    stats_table.to_csv(stats_csv, index=False)

    # Figure 06: per-variant deltas vs WT
    make_variant_delta_figure(d, figures_dir / "figure06_md_variant_delta_vs_wt.png")

    # Figure 07: class-level deltas
    deltas = d[d["classification"].isin(["deleterious", "neutral"])].copy()
    welch_ca = stats.ttest_ind(
        deltas[deltas["classification"] == "deleterious"]["delta_ca_vs_wt"],
        deltas[deltas["classification"] == "neutral"]["delta_ca_vs_wt"],
        equal_var=False,
    ).pvalue
    welch_rna = stats.ttest_ind(
        deltas[deltas["classification"] == "deleterious"]["delta_rna_vs_wt"],
        deltas[deltas["classification"] == "neutral"]["delta_rna_vs_wt"],
        equal_var=False,
    ).pvalue
    make_class_summary_figure(
        deltas,
        figures_dir / "figure07_md_class_delta_summary.png",
        pvals={"welch_ca": float(welch_ca), "welch_rna": float(welch_rna)},
    )

    # Summary markdown
    write_summary_md(bundle_dir / "MD_RESULTS_SUMMARY.md", d, stats_table)

    print(f"Wrote: {figures_dir / 'figure06_md_variant_delta_vs_wt.png'}")
    print(f"Wrote: {figures_dir / 'figure07_md_class_delta_summary.png'}")
    print(f"Wrote: {stats_csv}")
    print(f"Wrote: {bundle_dir / 'MD_RESULTS_SUMMARY.md'}")


if __name__ == "__main__":
    main()
