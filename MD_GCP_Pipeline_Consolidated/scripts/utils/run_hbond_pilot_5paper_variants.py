#!/usr/bin/env python3
"""Pilot H-bond comparison for the 5 paper variants vs WT counterpart residues.

This is a short-timescale, non-GROMACS pilot on currently available MD outputs.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis


VARIANTS = ["Y154H", "R160W", "R184C", "G204D", "C206F"]
MUT_RESID_LOCAL = {
    "Y154H": 4,
    "R160W": 10,
    "R184C": 34,
    "G204D": 54,
    "C206F": 56,
}

DONORS_SEL = "(protein or nucleic) and (name N* or name O* or name S*)"
HYDROGENS_SEL = "(protein or nucleic) and name H*"
ACCEPTORS_SEL = "(protein or nucleic) and (name O* or name N* or name S*) and not name H*"


def _resname_from_input(inputs_root: Path, variant: str, resid_local: int) -> str:
    u = mda.Universe(str(inputs_root / variant / "input_af.pdb"))
    ag = u.select_atoms(f"protein and segid A and resid {resid_local}")
    return ag.residues[0].resname if len(ag.residues) else "UNK"


def _hbond_count_stats(variant_output_dir: Path, resid_local: int) -> dict[str, float | int]:
    u = mda.Universe(
        str(variant_output_dir / "final_frame.pdb"),
        str(variant_output_dir / "production.dcd"),
    )
    residue_sel = f"protein and segid A and resid {resid_local}"
    partner_sel = f"(protein or nucleic) and not ({residue_sel})"

    h = HydrogenBondAnalysis(
        u,
        donors_sel=DONORS_SEL,
        hydrogens_sel=HYDROGENS_SEL,
        acceptors_sel=ACCEPTORS_SEL,
        between=[residue_sel, partner_sel],
        d_a_cutoff=3.0,
        d_h_a_angle_cutoff=150,
    )
    h.run(step=1)

    frames = h.results.hbonds[:, 0].astype(int) if h.results.hbonds.size else np.array([], dtype=int)
    counts = np.zeros(len(u.trajectory), dtype=int)
    for frame in frames:
        if 0 <= frame < len(counts):
            counts[frame] += 1

    return {
        "n_frames": int(len(counts)),
        "mean_hbonds": float(np.mean(counts)),
        "p0": float(np.mean(counts == 0)),
        "p1": float(np.mean(counts == 1)),
        "p2": float(np.mean(counts == 2)),
        "p3plus": float(np.mean(counts >= 3)),
        "max_hbonds": int(np.max(counts) if len(counts) else 0),
    }


def run_analysis(consolidated_root: Path, bundle_dir: Path) -> pd.DataFrame:
    inputs_root = consolidated_root / "runs" / "variants"
    outputs_root = consolidated_root / "runs" / "cloud_outputs" / "variants"
    out_rows: list[dict[str, object]] = []

    for variant in VARIANTS:
        resid_local = MUT_RESID_LOCAL[variant]
        mut_dir = outputs_root / variant / "input_af_cleaned"
        wt_dir = outputs_root / "WT" / "input_af_cleaned"

        mut_stats = _hbond_count_stats(mut_dir, resid_local)
        wt_stats = _hbond_count_stats(wt_dir, resid_local)

        row: dict[str, object] = {
            "variant": variant,
            "mut_resid_local": resid_local,
            "wt_resname": _resname_from_input(inputs_root, "WT", resid_local),
            "mut_resname": _resname_from_input(inputs_root, variant, resid_local),
        }
        for key, val in mut_stats.items():
            row[f"mut_{key}"] = val
        for key, val in wt_stats.items():
            row[f"wt_{key}"] = val
        row["delta_mean_hbonds"] = float(row["mut_mean_hbonds"]) - float(row["wt_mean_hbonds"])
        out_rows.append(row)

    df = pd.DataFrame(out_rows)
    out_table = bundle_dir / "tables" / "md_hbond_pilot_vs_wt_5paper_variants.csv"
    out_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_table, index=False)
    return df


def write_figure(df: pd.DataFrame, bundle_dir: Path) -> Path:
    x = np.arange(len(df))
    width = 0.38

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.bar(x - width / 2, df["wt_mean_hbonds"], width=width, label="WT counterpart residue", color="#2E7D32", edgecolor="black", linewidth=0.6)
    ax.bar(x + width / 2, df["mut_mean_hbonds"], width=width, label="Mutant residue", color="#C62828", edgecolor="black", linewidth=0.6)

    ymax = float(
        max(
            float(df["wt_mean_hbonds"].max()),
            float(df["mut_mean_hbonds"].max()),
            1.0,
        )
    )
    ylim_top = ymax * 1.22 + 0.12
    ax.set_ylim(0.0, ylim_top)

    for i, d in enumerate(df["delta_mean_hbonds"]):
        y_base = float(max(df.loc[i, "wt_mean_hbonds"], df.loc[i, "mut_mean_hbonds"]))
        y_text = min(y_base + (0.05 * ymax) + 0.03, ylim_top * 0.94)
        ax.text(
            i,
            y_text,
            f"Δ={d:+.2f}",
            ha="center",
            va="bottom",
            fontsize=9,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "alpha": 0.7, "edgecolor": "none"},
            clip_on=True,
        )

    labels = [f"{v}\n{wt}->{mut}" for v, wt, mut in zip(df["variant"], df["wt_resname"], df["mut_resname"])]
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Mean H-bond count per frame")
    ax.set_title("Pilot: Residue-level H-bond counts vs WT (short MD)")
    ax.grid(axis="y", alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()

    out_fig = bundle_dir / "figures" / "figure08_hbond_pilot_5paper_variants.png"
    out_fig.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fig, dpi=300)
    plt.close(fig)
    return out_fig


def write_summary(df: pd.DataFrame, bundle_dir: Path) -> Path:
    lines = [
        "# H-bond Pilot Summary (5 Paper Variants)",
        "",
        "This is a pilot comparison on short available trajectories, not the paper's 1000 ns GROMACS setup.",
        "",
        "## Directional comparison vs WT counterpart residues",
        "",
    ]
    for _, r in df.iterrows():
        direction = "higher" if r["delta_mean_hbonds"] > 0 else ("lower" if r["delta_mean_hbonds"] < 0 else "equal")
        lines.append(
            f"- {r['variant']}: mutant mean `{r['mut_mean_hbonds']:.2f}` vs WT `{r['wt_mean_hbonds']:.2f}` (Δ `{r['delta_mean_hbonds']:+.2f}`, {direction})"
        )
    lines += [
        "",
        "## Caveats",
        "",
        "- Uses short trajectories from current OpenMM cloud outputs.",
        "- Uses geometric H-bond detection in MDAnalysis, not `gmx hbond`.",
        "- No replicate-level uncertainty estimates yet.",
        "",
    ]
    out_md = bundle_dir / "HBOND_PILOT_SUMMARY.md"
    out_md.write_text("\n".join(lines) + "\n")
    return out_md


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--consolidated-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
        help="Path to MD_GCP_Pipeline_Consolidated root.",
    )
    parser.add_argument(
        "--bundle-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "presentation" / "compiled_20260429_mixed8plus2",
        help="Compiled bundle directory where outputs will be written.",
    )
    args = parser.parse_args()

    df = run_analysis(args.consolidated_root, args.bundle_dir)
    fig_path = write_figure(df, args.bundle_dir)
    md_path = write_summary(df, args.bundle_dir)
    print(f"Wrote: {args.bundle_dir / 'tables' / 'md_hbond_pilot_vs_wt_5paper_variants.csv'}")
    print(f"Wrote: {fig_path}")
    print(f"Wrote: {md_path}")


if __name__ == "__main__":
    main()
