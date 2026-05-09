#!/usr/bin/env python3
"""
Tau-style CA–CA distance analysis for the prep_run trajectories.

Outputs match the AF-MI style notebook:
  * Pairwise distance matrix (Å) written as CSV
  * Contact probability map (distance threshold in Å)
  * Pairwise distance probability distribution (CSV + plot)
  * 2x2 diagnostic panel (matrix heatmap, histogram, contact map, scatter)
  * Distance vs sequence separation table (Å) with mean curve
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable, Sequence, Tuple

import matplotlib

# Force a non-interactive backend to avoid Qt warnings on headless nodes.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--topology",
        type=Path,
        default=Path("AlphaFold-IDP/prep_run/input_af.pdb"),
        help="Path to the coarse-grained topology (default: input_af.pdb).",
    )
    parser.add_argument(
        "--trajectories",
        type=Path,
        nargs="+",
        default=[
            Path("AlphaFold-IDP/prep_run/output_0.dcd"),
            Path("AlphaFold-IDP/prep_run/output_1.dcd"),
        ],
        help="One or more DCD trajectories to analyze.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=25,
        help="Frames per chunk when streaming the DCD files.",
    )
    parser.add_argument(
        "--contact-cutoff",
        type=float,
        default=12.0,
        help="CA-CA contact threshold in Å (default: 12 Å).",
    )
    parser.add_argument(
        "--hist-max",
        type=float,
        default=500.0,
        help="Upper bound (Å) for the pairwise distance histogram.",
    )
    parser.add_argument(
        "--hist-bins",
        type=int,
        default=250,
        help="Number of bins for the pairwise distance histogram.",
    )
    parser.add_argument(
        "--scatter-pairs",
        type=int,
        default=600,
        help="Number of residue pairs to sample per chunk for the scatter plot.",
    )
    parser.add_argument(
        "--scatter-frames",
        type=int,
        default=5,
        help="Number of frames to sample per chunk for the scatter plot.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=2025,
        help="Random seed for scatter sampling.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("AlphaFold-IDP/prep_run/analysis_distances"),
        help="Directory where analysis outputs will be written.",
    )
    return parser.parse_args()


def expand_paths(paths: Sequence[Path]) -> list[Path]:
    return [p.expanduser().resolve() for p in paths]


def iter_chunks(
    traj_paths: Sequence[Path],
    topology: md.Topology,
    chunk_size: int,
) -> Iterable[md.Trajectory]:
    for traj_path in traj_paths:
        for chunk in md.iterload(str(traj_path), top=topology, chunk=chunk_size):
            yield chunk


def sample_distances(
    distances_A: np.ndarray,
    separations: np.ndarray,
    scatter_pairs: int,
    scatter_frames: int,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """Randomly sample distances for the scatter plot."""
    frames = distances_A.shape[0]
    n_pairs = distances_A.shape[1]
    frame_sel = rng.choice(frames, size=min(frames, scatter_frames), replace=False)
    pair_sel = rng.choice(n_pairs, size=min(n_pairs, scatter_pairs), replace=False)
    sub = distances_A[np.ix_(frame_sel, pair_sel)]
    sampled_sep = np.tile(separations[pair_sel], len(frame_sel))
    sampled_dist = sub.reshape(-1)
    return sampled_sep, sampled_dist


def compute_statistics(args: argparse.Namespace) -> dict:
    top_path = args.topology.expanduser().resolve()
    traj_paths = expand_paths(args.trajectories)
    template = md.load(str(top_path), standard_names=False)
    topology = template.topology
    residue_labels = [f"{res.resSeq}-{res.name}" for res in topology.residues]
    n_residues = template.n_atoms

    upper_i, upper_j = np.triu_indices(n_residues, k=1)
    pairs = np.vstack([upper_i, upper_j]).astype(np.int32).T
    n_pairs = pairs.shape[0]
    separations = np.abs(pairs[:, 0] - pairs[:, 1])

    contact_cutoff_nm = args.contact_cutoff / 10.0
    hist_edges = np.linspace(0.0, args.hist_max, args.hist_bins + 1, dtype=np.float64)

    hist_counts = np.zeros(args.hist_bins, dtype=np.float64)
    dist_sum = np.zeros(n_pairs, dtype=np.float64)
    dist_sumsq = np.zeros(n_pairs, dtype=np.float64)
    contact_counts = np.zeros(n_pairs, dtype=np.float64)

    samples_sep: list[np.ndarray] = []
    samples_dist: list[np.ndarray] = []
    rng = np.random.default_rng(args.seed)

    total_frames = 0
    chunk_idx = 0
    for chunk in iter_chunks(traj_paths, topology, args.chunk_size):
        chunk_idx += 1
        distances_nm = md.compute_distances(chunk, pairs)
        if np.isnan(distances_nm).any():
            raise ValueError(f"NaNs encountered in chunk {chunk_idx}")

        distances_A = distances_nm * 10.0
        dist_sum += distances_nm.sum(axis=0)
        dist_sumsq += np.square(distances_nm).sum(axis=0)
        contact_counts += (distances_nm <= contact_cutoff_nm).sum(axis=0)
        hist_counts += np.histogram(distances_A, bins=hist_edges)[0]
        total_frames += chunk.n_frames

        sampled_sep, sampled_dist = sample_distances(
            distances_A, separations, args.scatter_pairs, args.scatter_frames, rng
        )
        samples_sep.append(sampled_sep)
        samples_dist.append(sampled_dist)

    if total_frames == 0:
        raise RuntimeError("No frames were processed from the supplied trajectories.")

    samples = {
        "separations": np.concatenate(samples_sep),
        "distances_A": np.concatenate(samples_dist),
    }

    stats = {
        "topology": top_path,
        "traj_paths": traj_paths,
        "residue_labels": residue_labels,
        "n_residues": n_residues,
        "pairs": pairs,
        "separations": separations,
        "hist_edges": hist_edges,
        "hist_counts": hist_counts,
        "dist_sum": dist_sum,
        "dist_sumsq": dist_sumsq,
        "contact_counts": contact_counts,
        "total_frames": total_frames,
        "samples": samples,
    }
    return stats


def summarize_results(stats: dict, args: argparse.Namespace, outdir: Path) -> None:
    residue_labels = stats["residue_labels"]
    n_residues = stats["n_residues"]
    pairs = stats["pairs"]
    separations = stats["separations"]
    total_frames = stats["total_frames"]

    dist_sum = stats["dist_sum"]
    dist_sumsq = stats["dist_sumsq"]
    contact_counts = stats["contact_counts"]
    mean_dist_nm = dist_sum / total_frames
    var_dist_nm = np.clip(dist_sumsq / total_frames - np.square(mean_dist_nm), 0.0, None)
    std_dist_nm = np.sqrt(var_dist_nm)
    contact_prob = contact_counts / total_frames

    mean_matrix_A = np.zeros((n_residues, n_residues), dtype=np.float64)
    contact_matrix = np.zeros_like(mean_matrix_A)

    mean_matrix_A[pairs[:, 0], pairs[:, 1]] = mean_dist_nm * 10.0
    mean_matrix_A[pairs[:, 1], pairs[:, 0]] = mean_matrix_A[pairs[:, 0], pairs[:, 1]]
    contact_matrix[pairs[:, 0], pairs[:, 1]] = contact_prob
    contact_matrix[pairs[:, 1], pairs[:, 0]] = contact_prob

    pd.DataFrame(mean_matrix_A, index=residue_labels, columns=residue_labels).to_csv(
        outdir / "pairwise_mean_distance_angstrom.csv"
    )
    pd.DataFrame(contact_matrix, index=residue_labels, columns=residue_labels).to_csv(
        outdir / "contact_map_probability.csv"
    )

    hist_edges = stats["hist_edges"]
    hist_counts = stats["hist_counts"]
    hist_centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
    hist_prob = hist_counts / hist_counts.sum()
    pd.DataFrame({"distance_A": hist_centers, "probability": hist_prob}).to_csv(
        outdir / "pairwise_distance_probability_distribution.csv", index=False
    )

    seq_measurements = np.bincount(separations, minlength=n_residues) * total_frames
    sep_sum = np.bincount(separations, weights=dist_sum, minlength=n_residues)
    sep_sumsq = np.bincount(separations, weights=dist_sumsq, minlength=n_residues)
    valid = (seq_measurements > 0) & (np.arange(n_residues) > 0)

    mean_sep_nm = np.zeros(n_residues)
    std_sep_nm = np.zeros(n_residues)
    mean_sep_nm[valid] = sep_sum[valid] / seq_measurements[valid]
    std_sep_nm[valid] = np.sqrt(
        np.clip(sep_sumsq[valid] / seq_measurements[valid] - np.square(mean_sep_nm[valid]), 0.0, None)
    )

    seq_sep = np.arange(1, n_residues, dtype=int)
    sep_df = pd.DataFrame(
        {
            "sequence_separation": seq_sep,
            "mean_distance_A": mean_sep_nm[1:] * 10.0,
            "std_distance_A": std_sep_nm[1:] * 10.0,
            "contact_probability": np.bincount(
                separations, weights=contact_counts, minlength=n_residues
            )[1:]
            / seq_measurements[1:],
        }
    )
    sep_df.to_csv(outdir / "distance_vs_sequence_separation.csv", index=False)

    metadata = {
        "topology": str(stats["topology"]),
        "trajectories": [str(p) for p in stats["traj_paths"]],
        "n_residues": int(n_residues),
        "n_pairs": int(pairs.shape[0]),
        "total_frames": int(total_frames),
        "contact_cutoff_A": args.contact_cutoff,
        "hist_max_A": args.hist_max,
        "hist_bins": args.hist_bins,
        "chunk_size": args.chunk_size,
        "scatter_pairs_per_chunk": args.scatter_pairs,
        "scatter_frames_per_chunk": args.scatter_frames,
        "seed": args.seed,
    }
    with open(outdir / "analysis_summary.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)

    plots = {
        "hist": (hist_centers, hist_prob),
        "mean_matrix_A": mean_matrix_A,
        "contact_matrix": contact_matrix,
        "seq_sep_df": sep_df,
    }
    generate_figures(plots, stats["samples"], args, outdir)


def generate_figures(
    plots: dict,
    samples: dict,
    args: argparse.Namespace,
    outdir: Path,
) -> None:
    mean_matrix_A = plots["mean_matrix_A"]
    contact_matrix = plots["contact_matrix"]
    hist_centers, hist_prob = plots["hist"]
    sep_df = plots["seq_sep_df"]
    sample_sep = samples["separations"]
    sample_dist = samples["distances_A"]

    max_dist = np.percentile(mean_matrix_A[np.triu_indices_from(mean_matrix_A, k=1)], 95)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    ax1 = axes[0, 0]
    im1 = ax1.imshow(
        mean_matrix_A,
        cmap="viridis",
        origin="lower",
        interpolation="nearest",
        vmin=0.0,
        vmax=max_dist,
    )
    ax1.set_title("Mean CA-CA Distance (Å)", fontsize=14, fontweight="bold")
    ax1.set_xlabel("Residue j")
    ax1.set_ylabel("Residue i")
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, label="Å")

    ax2 = axes[0, 1]
    width = hist_centers[1] - hist_centers[0]
    ax2.bar(hist_centers, hist_prob, width=width, color="steelblue", alpha=0.8, edgecolor="black")
    mean_A = np.sum(hist_centers * hist_prob)
    cdf = np.cumsum(hist_prob)
    median_A = np.interp(0.5, cdf, hist_centers)
    ax2.axvline(mean_A, color="red", linestyle="--", label=f"Mean: {mean_A:.1f} Å")
    ax2.axvline(median_A, color="green", linestyle="-.", label=f"Median: {median_A:.1f} Å")
    ax2.set_xlabel("Distance (Å)")
    ax2.set_ylabel("Probability")
    ax2.set_title("Pairwise Distance Distribution", fontsize=14, fontweight="bold")
    ax2.legend()

    ax3 = axes[1, 0]
    im3 = ax3.imshow(
        contact_matrix,
        cmap="Reds",
        interpolation="nearest",
        origin="lower",
        vmin=0.0,
        vmax=1.0,
    )
    ax3.set_xlabel("Residue j")
    ax3.set_ylabel("Residue i")
    ax3.set_title(f"Contact Map (≤ {args.contact_cutoff:.1f} Å)", fontsize=14, fontweight="bold")
    plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04, label="P(contact)")

    ax4 = axes[1, 1]
    ax4.scatter(
        sample_sep,
        sample_dist,
        s=4,
        alpha=0.08,
        color="purple",
        rasterized=True,
    )
    ax4.plot(
        sep_df["sequence_separation"],
        sep_df["mean_distance_A"],
        color="red",
        linewidth=2,
        label="Mean distance",
    )
    ax4.set_xlabel("Sequence separation |i - j|")
    ax4.set_ylabel("Distance (Å)")
    ax4.set_yscale("log")
    ax4.set_title("Distance vs Sequence Separation", fontsize=14, fontweight="bold")
    ax4.grid(True, alpha=0.3)
    ax4.legend()

    fig.tight_layout()
    fig.savefig(outdir / "pairwise_distance_panels.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig2, ax = plt.subplots(figsize=(10, 6))
    ax.plot(hist_centers, hist_prob, color="navy", linewidth=2, label="Probability")
    ax.fill_between(hist_centers, hist_prob, color="navy", alpha=0.25)
    ax.set_xlabel("CA-CA Distance (Å)")
    ax.set_ylabel("Probability Density")
    ax.set_title("Pairwise Distance Probability Distribution", fontsize=16, fontweight="bold")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig2.tight_layout()
    fig2.savefig(outdir / "distance_probability_distribution.png", dpi=300, bbox_inches="tight")
    plt.close(fig2)


def main() -> None:
    args = parse_args()
    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    stats = compute_statistics(args)
    summarize_results(stats, args, outdir)
    print(f"Wrote Tau-style analysis outputs to {outdir}")


if __name__ == "__main__":
    main()

