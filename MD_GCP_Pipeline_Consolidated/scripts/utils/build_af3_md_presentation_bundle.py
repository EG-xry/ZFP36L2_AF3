#!/usr/bin/env python3
"""
Build a presentation-ready AF3 + MD summary bundle.

Outputs include:
- merged tables for AF3 static and MD minimal runs
- new AF3-focused figures (run stability, class separation, AF3↔MD bridge)
- copied key legacy figures from prior analyses
- concise markdown reports of findings and next AF3-only steps
"""
from __future__ import annotations

import argparse
import math
import shutil
from datetime import date
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import rankdata

# -----------------------------------------------------------------------------
# Project constants
# -----------------------------------------------------------------------------
CLASS_MAP: Dict[str, str] = {
    "WT": "WT",
    "Y154H": "deleterious",
    "R160E": "deleterious",
    "R160W": "deleterious",
    "R184C": "deleterious",
    "G204D": "deleterious",
    "C206F": "deleterious",
    "C208A": "deleterious",  # control outside paper, intentionally retained
    "S185N": "neutral",
    "T187S": "neutral",
    "D210E": "neutral",
}

SELECTED_VARIANTS: List[str] = [
    "WT",
    "Y154H",
    "R160E",
    "R160W",
    "R184C",
    "S185N",
    "T187S",
    "G204D",
    "C206F",
    "C208A",
    "D210E",
]

EXCLUDED_VARIANTS: List[str] = ["H212A"]

CLASS_ORDER = ["WT", "neutral", "deleterious"]
CLASS_COLORS = {
    "WT": "#2E7D32",
    "neutral": "#1565C0",
    "deleterious": "#C62828",
}


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------
def variant_order_for_plots(df: pd.DataFrame) -> List[str]:
    present = set(df["variant"].astype(str))
    order = ["WT"]
    order += sorted(v for v in present if CLASS_MAP.get(v) == "neutral")
    order += sorted(v for v in present if CLASS_MAP.get(v) == "deleterious")
    extras = sorted(v for v in present if v not in order)
    return order + extras


def class_color(variant: str) -> str:
    return CLASS_COLORS.get(CLASS_MAP.get(variant, ""), "#616161")


def extract_variant_from_complex_base(base: str) -> str:
    text = base.lower()
    import re

    hits = re.findall(r"(?<![a-z0-9])[a-z]\d{3}[a-z](?![a-z0-9])", text)
    if hits:
        return hits[0].upper()
    if text.startswith("fold_zfp36l2_wt_"):
        return "WT"
    return "UNKNOWN"


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def safe_float(x) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def roc_auc_from_scores(y_true: np.ndarray, scores: np.ndarray) -> float:
    """AUC via Mann-Whitney U / rank formulation."""
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)
    pos = scores[y_true == 1]
    neg = scores[y_true == 0]
    n_pos = len(pos)
    n_neg = len(neg)
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = rankdata(np.concatenate([pos, neg]), method="average")
    rank_sum_pos = np.sum(ranks[:n_pos])
    auc = (rank_sum_pos - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return float(auc)


def auc_bootstrap_ci(
    y_true: np.ndarray,
    scores: np.ndarray,
    n_boot: int = 5000,
    seed: int = 7,
) -> Tuple[float, float]:
    """Stratified bootstrap CI for AUC."""
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    if len(pos_idx) < 2 or len(neg_idx) < 2:
        return float("nan"), float("nan")

    vals = []
    for _ in range(n_boot):
        b_pos = rng.choice(pos_idx, size=len(pos_idx), replace=True)
        b_neg = rng.choice(neg_idx, size=len(neg_idx), replace=True)
        idx = np.concatenate([b_pos, b_neg])
        auc = roc_auc_from_scores(y_true[idx], scores[idx])
        if math.isfinite(auc):
            vals.append(auc)
    if not vals:
        return float("nan"), float("nan")
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def auc_permutation_pvalue(
    y_true: np.ndarray,
    scores: np.ndarray,
    n_perm: int = 20000,
    seed: int = 17,
) -> float:
    """One-sided permutation p-value: P(AUC_perm >= AUC_obs)."""
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true).astype(int)
    scores = np.asarray(scores).astype(float)
    obs = roc_auc_from_scores(y_true, scores)
    if not math.isfinite(obs):
        return float("nan")
    count = 0
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        auc = roc_auc_from_scores(y_perm, scores)
        if auc >= obs:
            count += 1
    # +1 correction
    return float((count + 1) / (n_perm + 1))


def build_af3_composite_score(df: pd.DataFrame) -> pd.DataFrame:
    """Higher score = further from WT in an unfavorable direction."""
    out = df.copy()
    wt = out.loc[out["variant"] == "WT"].iloc[0]

    wt_iptm = safe_float(wt["pair_ipTM_protein_rna"])
    wt_iptm_sd = max(safe_float(wt["pair_ipTM_protein_rna_sd"]), 1e-6)
    wt_pae = safe_float(wt["pair_pae_min_cross_mean"])
    wt_pae_sd = max(safe_float(wt["pair_pae_min_cross_mean_sd"]), 1e-6)
    wt_plddt_min = safe_float(wt["pLDDT_min"])
    wt_plddt_min_sd = max(safe_float(wt["pLDDT_min_sd"]), 1e-6)

    out["z_bad_iptm"] = (wt_iptm - out["pair_ipTM_protein_rna"]) / wt_iptm_sd
    out["z_bad_pae"] = (out["pair_pae_min_cross_mean"] - wt_pae) / wt_pae_sd
    out["z_bad_plddt_min"] = (wt_plddt_min - out["pLDDT_min"]) / wt_plddt_min_sd
    out["af3_composite_badness"] = out[["z_bad_iptm", "z_bad_pae", "z_bad_plddt_min"]].mean(axis=1)
    return out


# -----------------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------------
def load_per_model_pair_metrics(workspace_root: Path) -> pd.DataFrame:
    """Load AF3 per-model pair metrics from summary_confidences JSON files."""
    import json
    import re

    data_root = workspace_root / "Data"
    json_files = sorted(data_root.glob("**/*_summary_confidences_*.json"))
    if not json_files:
        raise RuntimeError(f"No summary confidence JSON files found under {data_root}")

    # Stable run indexing by folder path avoids '(1)' vs no-suffix collisions.
    run_folder_rank: Dict[Tuple[str, str], int] = {}
    per_variant_folders: Dict[str, List[str]] = {}
    parsed = []
    for p in json_files:
        m = re.match(r"(.+)_summary_confidences_(\d+)\.json$", p.name)
        if not m:
            continue
        base = m.group(1)
        local_model = int(m.group(2))
        variant = extract_variant_from_complex_base(base)
        if variant not in SELECTED_VARIANTS:
            continue
        folder = str(p.parent.resolve())
        per_variant_folders.setdefault(variant, [])
        if folder not in per_variant_folders[variant]:
            per_variant_folders[variant].append(folder)
        parsed.append((p, base, variant, local_model, folder))

    for variant, folders in per_variant_folders.items():
        for i, folder in enumerate(sorted(folders), start=1):
            run_folder_rank[(variant, folder)] = i

    rows = []
    for p, base, variant, local_model, folder in parsed:
        try:
            data = json.loads(p.read_text())
        except Exception:
            continue

        chain_ptm = data.get("chain_ptm", None)
        chain_pair_iptm = data.get("chain_pair_iptm", None)
        pae_matrix = data.get("chain_pair_pae_min", None)

        valid = []
        if isinstance(chain_ptm, list):
            valid = [i for i, v in enumerate(chain_ptm) if v is not None]
        if not valid:
            continue

        pidx = max(valid, key=lambda i: float(chain_ptm[i]))
        candidates = [i for i in valid if i != pidx]
        ridx = None

        def _pae_ok(i: int, j: int) -> bool:
            try:
                return (
                    isinstance(pae_matrix, list)
                    and pae_matrix[i][j] is not None
                    and pae_matrix[j][i] is not None
                )
            except Exception:
                return False

        good = [j for j in candidates if _pae_ok(pidx, j)]
        if len(good) == 1:
            ridx = good[0]
        elif len(good) > 1 and isinstance(chain_pair_iptm, list):
            ridx = max(
                good,
                key=lambda j: float(
                    chain_pair_iptm[pidx][j] if chain_pair_iptm[pidx][j] is not None else -1.0
                ),
            )
        elif candidates:
            ridx = candidates[0]

        pair_iptm = float("nan")
        pair_pae = float("nan")
        if ridx is not None and isinstance(chain_pair_iptm, list):
            try:
                v = chain_pair_iptm[pidx][ridx]
                pair_iptm = float(v) if v is not None else float("nan")
            except Exception:
                pass
        if ridx is not None and isinstance(pae_matrix, list):
            try:
                v1 = pae_matrix[pidx][ridx]
                v2 = pae_matrix[ridx][pidx]
                vals = [float(v) for v in [v1, v2] if v is not None]
                if vals:
                    pair_pae = float(np.mean(vals))
            except Exception:
                pass

        run_idx = run_folder_rank[(variant, folder)]
        model_global_idx = (run_idx - 1) * 5 + local_model

        rows.append(
            {
                "complex_base": base,
                "variant": variant,
                "classification": CLASS_MAP[variant],
                "run_idx": run_idx,
                "model_global_idx": int(model_global_idx),
                "model_local_idx": int(local_model),
                "iptm": safe_float(data.get("iptm")),
                "ptm": safe_float(data.get("ptm")),
                "pair_ipTM_protein_rna": pair_iptm,
                "pair_pae_min_cross_mean": pair_pae,
                "ranking_score": safe_float(data.get("ranking_score")),
            }
        )

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No AF3 per-model metrics parsed from Data/* summary_confidences JSON files")

    df = df[df["variant"].isin(SELECTED_VARIANTS)].copy()
    return df


def load_af3_variant_summary(workspace_root: Path) -> pd.DataFrame:
    src = workspace_root / "results" / "variant_avg_all.csv"
    df = pd.read_csv(src)
    df = df.rename(columns={"Protein": "variant", "ProteinExpClass": "classification_source"})
    df = df[df["variant"].isin(SELECTED_VARIANTS)].copy()
    df["classification"] = df["variant"].map(CLASS_MAP)
    df.sort_values("variant", inplace=True)
    return df


def load_md_tables(consolidated_root: Path) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    md_qc = pd.read_csv(consolidated_root / "manifests" / "cloud_variant_analysis.csv")

    md_ca = pd.read_csv(consolidated_root / "manifests" / "cloud_variant_structural_drift.csv")
    md_rna = pd.read_csv(consolidated_root / "manifests" / "cloud_variant_structural_drift_rna.csv")

    for df in (md_qc, md_ca, md_rna):
        if "classification" in df.columns:
            df["classification"] = df["variant"].map(CLASS_MAP).fillna(df["classification"])

    return md_qc, md_ca, md_rna


def _read_env(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path.exists():
        return out
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.strip()
        if not line or "=" not in line:
            continue
        key, val = line.split("=", 1)
        out[key.strip()] = val.strip()
    return out


def _production_profile(final_step: float | int | str | None) -> str:
    try:
        step = int(float(final_step))  # type: ignore[arg-type]
    except Exception:
        return "unknown"
    if step >= 100000:
        return "extended_0.23ns_like"
    if step >= 50000:
        return "short_0.10ns_like"
    if step >= 25000:
        return "short_0.05ns_like"
    return "very_short_or_incomplete"


def load_md_run_provenance(consolidated_root: Path, md_qc: pd.DataFrame) -> pd.DataFrame:
    root = consolidated_root / "runs" / "cloud_outputs" / "variants"
    qc_lookup = md_qc.set_index("variant").to_dict(orient="index")

    rows: List[Dict[str, object]] = []
    for variant in SELECTED_VARIANTS:
        env_path = root / variant / "task_metadata.env"
        meta = _read_env(env_path)
        rel_input = meta.get("REL_INPUT_PDB", "")
        task_count = meta.get("BATCH_TASK_COUNT", "")
        final_step = qc_lookup.get(variant, {}).get("prod_final_step", "")

        source_group = "no_md_run"
        if env_path.exists():
            source_group = "main_batch"
            if "variants_failed_rerun/" in rel_input or task_count == "2":
                source_group = "fallback_short_rerun"

        rows.append(
            {
                "variant": variant,
                "classification": CLASS_MAP.get(variant, ""),
                "task_metadata_present": env_path.exists(),
                "rel_input_pdb": rel_input if rel_input else "no_md_run",
                "batch_task_index": meta.get("BATCH_TASK_INDEX", "") if env_path.exists() else "unknown",
                "batch_task_count": task_count if task_count else ("unknown" if not env_path.exists() else ""),
                "run_source_group": source_group,
                "prod_final_step": final_step,
                "production_profile": _production_profile(final_step),
            }
        )

    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def plot_af3_variant_metric_bars(df_variant: pd.DataFrame, out_path: Path) -> None:
    order = variant_order_for_plots(df_variant)
    d = df_variant.set_index("variant").reindex(order).reset_index()
    x = np.arange(len(d))
    colors = [class_color(v) for v in d["variant"]]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    axes[0].bar(
        x,
        d["pair_ipTM_protein_rna"],
        yerr=d["pair_ipTM_protein_rna_sd"],
        color=colors,
        edgecolor="black",
        linewidth=0.6,
        capsize=3,
    )
    axes[0].set_ylabel("pair ipTM")
    axes[0].set_title("AF3 static: Protein-RNA interface confidence by variant")
    axes[0].grid(axis="y", alpha=0.3)

    axes[1].bar(
        x,
        d["pair_pae_min_cross_mean"],
        yerr=d["pair_pae_min_cross_mean_sd"],
        color=colors,
        edgecolor="black",
        linewidth=0.6,
        capsize=3,
    )
    axes[1].set_ylabel("PAE cross-mean (A)")
    axes[1].set_xlabel("Variant")
    axes[1].grid(axis="y", alpha=0.3)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(d["variant"], rotation=30, ha="right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_af3_model0_run_stability(df_models: pd.DataFrame, out_path: Path) -> None:
    d = df_models[df_models["model_local_idx"] == 0].copy()
    order = variant_order_for_plots(d.rename(columns={"variant": "variant"}))
    d["variant"] = pd.Categorical(d["variant"], categories=order, ordered=True)
    d.sort_values(["variant", "run_idx"], inplace=True)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    # per-variant run traces: pair ipTM
    for variant, sub in d.groupby("variant", observed=True):
        xv = np.full(len(sub), order.index(str(variant)), dtype=float)
        jitter = (sub["run_idx"].to_numpy() - sub["run_idx"].mean()) * 0.02
        axes[0].scatter(xv + jitter, sub["pair_ipTM_protein_rna"], color=class_color(str(variant)), s=24, alpha=0.85)
    med0 = d.groupby("variant", observed=True)["pair_ipTM_protein_rna"].median().reindex(order)
    axes[0].plot(np.arange(len(order)), med0.values, color="black", linewidth=1.2, label="Per-variant median")
    axes[0].set_ylabel("pair ipTM (model0 across runs)")
    axes[0].grid(axis="y", alpha=0.3)
    axes[0].legend(loc="best", fontsize=9)

    for variant, sub in d.groupby("variant", observed=True):
        xv = np.full(len(sub), order.index(str(variant)), dtype=float)
        jitter = (sub["run_idx"].to_numpy() - sub["run_idx"].mean()) * 0.02
        axes[1].scatter(xv + jitter, sub["pair_pae_min_cross_mean"], color=class_color(str(variant)), s=24, alpha=0.85)
    med1 = d.groupby("variant", observed=True)["pair_pae_min_cross_mean"].median().reindex(order)
    axes[1].plot(np.arange(len(order)), med1.values, color="black", linewidth=1.2, label="Per-variant median")
    axes[1].set_ylabel("PAE cross-mean (A), model0")
    axes[1].set_xlabel("Variant")
    axes[1].grid(axis="y", alpha=0.3)

    axes[1].set_xticks(np.arange(len(order)))
    axes[1].set_xticklabels(order, rotation=30, ha="right")

    fig.suptitle("AF3 run-to-run stability (top-ranked model per run)", y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_af3_class_separation(df_score_mut: pd.DataFrame, out_path: Path) -> None:
    classes = ["neutral", "deleterious"]
    fig, ax = plt.subplots(figsize=(7, 5))

    box_data = [
        df_score_mut.loc[df_score_mut["classification"] == c, "af3_composite_badness"].to_numpy()
        for c in classes
    ]
    b = ax.boxplot(
        box_data,
        tick_labels=classes,
        patch_artist=True,
        showmeans=True,
        meanprops={
            "marker": "^",
            "markerfacecolor": "#2E7D32",
            "markeredgecolor": "black",
            "markersize": 8,
        },
    )
    for patch, c in zip(b["boxes"], classes):
        patch.set_facecolor(CLASS_COLORS[c])
        patch.set_alpha(0.35)

    # Estimate box vertical spans so labels can avoid sitting on top of the box body.
    class_box_y: Dict[str, Tuple[float, float]] = {}
    for c, ys in zip(classes, box_data):
        q1 = float(np.percentile(ys, 25))
        q3 = float(np.percentile(ys, 75))
        class_box_y[c] = (q1, q3)

    rng = np.random.default_rng(13)
    placed_xy: List[Tuple[float, float]] = []
    for i, c in enumerate(classes, start=1):
        sub = (
            df_score_mut.loc[df_score_mut["classification"] == c, ["variant", "af3_composite_badness"]]
            .sort_values("af3_composite_badness")
            .reset_index(drop=True)
        )
        ys = sub["af3_composite_badness"].to_numpy()
        xs = i + rng.normal(0, 0.04, size=len(ys))
        ax.scatter(xs, ys, s=45, color=CLASS_COLORS[c], edgecolor="black", linewidth=0.4, zorder=3)
        q1, q3 = class_box_y[c]
        box_x_min, box_x_max = float(i - 0.25), float(i + 0.25)
        for x, y, name in zip(xs, ys, sub["variant"]):
            # Keep labels close to points, but choose nearby offsets that avoid
            # the dashed baseline and the core box area.
            candidates = [
                (0.026, 0.028),
                (-0.026, 0.028),
                (0.026, -0.028),
                (-0.026, -0.028),
                (0.034, 0.000),
                (-0.034, 0.000),
                (0.000, 0.040),
                (0.000, -0.040),
            ]
            # Manual preference: S185N label should sit above/left when possible.
            if str(name).upper() == "S185N":
                candidates = [
                    (-0.020, 0.040),
                    (0.000, 0.045),
                    (-0.030, 0.028),
                    (-0.028, -0.024),
                    (0.026, 0.028),
                    (0.026, -0.028),
                ]

            chosen = None
            for dx, dy in candidates:
                tx = float(x + dx)
                ty = float(y + dy)
                # Avoid placing text directly on dashed y=0 line.
                if abs(ty) < 0.08:
                    continue
                # Avoid core box body overlap.
                in_box_body = (box_x_min <= tx <= box_x_max) and ((q1 - 0.02) <= ty <= (q3 + 0.02))
                if in_box_body:
                    continue
                # Avoid label-label collisions.
                collides = any(abs(tx - px) < 0.06 and abs(ty - py) < 0.06 for px, py in placed_xy)
                if collides:
                    continue
                chosen = (tx, ty, dx)
                break

            if chosen is None:
                # Fallback: keep it near point and shift away from y=0.
                dx, dy = (0.028, 0.028 if y >= 0 else -0.028)
                tx = float(x + dx)
                ty = float(y + dy)
                if abs(ty) < 0.08:
                    ty += 0.10 if ty >= 0 else -0.10
                chosen = (tx, ty, dx)

            tx, ty, dx = chosen
            placed_xy.append((tx, ty))
            ha = "left" if dx >= 0 else "right"
            ax.text(
                tx,
                ty,
                name,
                fontsize=8,
                va="center",
                ha=ha,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.75, "pad": 0.2},
                zorder=4,
            )

    ax.axhline(0.0, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_ylabel("AF3 composite badness (WT-referenced z-score)")
    ax.set_title("Static AF3 class separation (variant-level means)")
    mean_handle = Line2D(
        [0],
        [0],
        marker="^",
        linestyle="None",
        markerfacecolor="#2E7D32",
        markeredgecolor="black",
        markersize=8,
        label="Class mean",
    )
    ax.legend(handles=[mean_handle], loc="upper left", fontsize=8, framealpha=0.9)
    ax.grid(axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_af3_md_bridge(df_bridge: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)

    for _, row in df_bridge.iterrows():
        c = class_color(str(row["variant"]))
        axes[0].scatter(row["af3_composite_badness"], row["ca_rmsd_a_final_vs_start"], s=60, color=c, edgecolor="black", linewidth=0.5)
        axes[0].text(row["af3_composite_badness"] + 0.02, row["ca_rmsd_a_final_vs_start"], row["variant"], fontsize=8)

        axes[1].scatter(row["af3_composite_badness"], row["p_rmsd_b_final_vs_start"], s=60, color=c, edgecolor="black", linewidth=0.5)
        axes[1].text(row["af3_composite_badness"] + 0.02, row["p_rmsd_b_final_vs_start"], row["variant"], fontsize=8)

    axes[0].set_ylabel("MD protein drift (CA RMSD A)")
    axes[1].set_ylabel("MD RNA drift (P RMSD A)")
    for ax in axes:
        ax.set_xlabel("AF3 composite badness")
        ax.grid(alpha=0.3)

    fig.suptitle("Bridge: AF3 static score vs short MD drift", y=0.99)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_md_qc_overview(df_md: pd.DataFrame, out_path: Path) -> None:
    d = df_md.copy()
    d["classification"] = d["variant"].map(CLASS_MAP)
    order = variant_order_for_plots(d)
    d = d.set_index("variant").reindex(order).reset_index()
    x = np.arange(len(d))
    colors = [class_color(v) for v in d["variant"]]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

    axes[0].bar(x, d["prod_speed_mean_nsday"], color=colors, edgecolor="black", linewidth=0.6)
    axes[0].set_ylabel("Speed (ns/day)")
    axes[0].set_title("Cloud minimal-MD QC snapshot")
    axes[0].grid(axis="y", alpha=0.3)

    axes[1].bar(x, d["prod_temp_mean_k"], color=colors, edgecolor="black", linewidth=0.6)
    axes[1].axhline(300.0, color="black", linestyle="--", linewidth=0.8, alpha=0.7)
    axes[1].set_ylabel("Production mean T (K)")
    axes[1].set_xlabel("Variant")
    axes[1].grid(axis="y", alpha=0.3)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(d["variant"], rotation=30, ha="right")

    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


# -----------------------------------------------------------------------------
# Reporting
# -----------------------------------------------------------------------------
def write_report(
    report_path: Path,
    bundle_rel: str,
    auc: float,
    auc_lo: float,
    auc_hi: float,
    auc_p: float,
    df_score: pd.DataFrame,
    df_md: pd.DataFrame,
    df_md_provenance: pd.DataFrame,
) -> None:
    top_bad = df_score.sort_values("af3_composite_badness", ascending=False).head(3)
    top_good = df_score.sort_values("af3_composite_badness", ascending=True).head(3)

    done = int((df_md["qc_status"] == "complete_minimal_md").sum()) if "qc_status" in df_md.columns else 0
    total = int(len(df_md))
    n_neutral = int((df_score["classification"] == "neutral").sum())
    n_del = int((df_score["classification"] == "deleterious").sum())
    df_md_present = df_md_provenance[df_md_provenance["task_metadata_present"] == True].copy()
    n_main = int((df_md_present["run_source_group"] == "main_batch").sum())
    n_fallback = int((df_md_present["run_source_group"] == "fallback_short_rerun").sum())
    fallback_variants = (
        df_md_present.loc[df_md_present["run_source_group"] == "fallback_short_rerun", "variant"]
        .dropna()
        .astype(str)
        .tolist()
    )
    fallback_text = ", ".join(fallback_variants) if fallback_variants else "none"

    lines = [
        "# AF3 + MD Compiled Results",
        "",
        f"Bundle folder: `{bundle_rel}`",
        f"Date: `{date.today().isoformat()}`",
        "",
        "## Scope",
        "",
        "- Included variants: WT + 10 variants (includes controls C208A and R160E).",
        "- Excluded variant: H212A.",
        "- AF3 static source: existing multi-run outputs in `Data/*` and summary tables in `results/*`.",
        "- MD source: cloud minimal-MD outputs in `MD_GCP_Pipeline_Consolidated/manifests/*`.",
        f"- MD provenance split used here: `{n_main}` from main batch + `{n_fallback}` from fallback short rerun (`{fallback_text}`).",
        "",
        "## Key Findings",
        "",
        f"- AF3 composite deleterious-vs-neutral AUC (variant-level means): `{auc:.3f}` (95% bootstrap CI `{auc_lo:.3f}` to `{auc_hi:.3f}`), permutation p-value `{auc_p:.4f}`.",
        f"- Minimal-MD completion status in current table: `{done}/{total}` variants marked `complete_minimal_md`.",
        "- AF3 composite score is WT-referenced and directional (higher = worse than WT on interface confidence).",
        "",
        "### Most WT-deviated variants by AF3 composite",
        "",
    ]

    for _, row in top_bad.iterrows():
        lines.append(
            f"- {row['variant']}: score `{row['af3_composite_badness']:.2f}` ({row['classification']})"
        )

    lines += [
        "",
        "### Closest to WT by AF3 composite",
        "",
    ]

    for _, row in top_good.iterrows():
        lines.append(
            f"- {row['variant']}: score `{row['af3_composite_badness']:.2f}` ({row['classification']})"
        )

    lines += [
        "",
        "## Files",
        "",
        "- `tables/af3_variant_summary_selected.csv`: canonical AF3 static means/SD (all models) for selected variants.",
        "- `tables/af3_per_model_pair_metrics_selected.csv`: per-model pair metrics parsed from AF3 summary JSONs.",
        "- `tables/combined_af3_md_selected.csv`: AF3 + MD merged table for cross-modal comparisons.",
        "- `tables/md_run_provenance_selected.csv`: per-variant MD run-source provenance (main batch vs fallback rerun).",
        "- `figures/figure01_af3_pair_metrics_by_variant.png`: interface confidence by variant.",
        "- `figures/figure02_af3_model0_run_stability.png`: run-to-run spread of top AF3 model.",
        "- `figures/figure03_af3_class_separation.png`: deleterious vs neutral separation.",
        "- `figures/figure04_af3_md_bridge.png`: AF3 score vs MD drift.",
        "- `figures/figure05_md_qc_overview.png`: MD runtime QC summary.",
        "",
        "## Interpretation Limits",
        "",
        f"- AF3 classification check uses small N at the variant level (neutral n={n_neutral}, deleterious n={n_del}).",
        "- Current MD is short-timescale minimal runs and not a full 1000 ns protocol replication.",
        "",
    ]

    report_path.write_text("\n".join(lines) + "\n")


def write_af3_expansion_plan(path: Path) -> None:
    lines = [
        "# AF3 Static Expansion Plan (No New MD Required)",
        "",
        "## Priority 1: High-value analyses on existing runs",
        "",
        "1. Bootstrap stability profiles per variant",
        "Run bootstrap CIs for pair ipTM and pair PAE across existing model pools to quantify uncertainty, not just means.",
        "2. Model-rank sensitivity",
        "Compare conclusions using model_0 only versus all-model pooled metrics to test rank robustness.",
        "3. Run-level variance decomposition",
        "Split variance into between-run and within-run components to identify variants with unstable AF3 behavior.",
        "4. Paired metric discordance map",
        "Flag variants where ipTM improves but PAE worsens (or vice versa), suggesting fragile interfaces.",
        "",
        "## Priority 2: New AF3 static runs to request",
        "",
        "1. Add replicate runs where any variant has lower run depth than the rest",
        "Balancing run counts improves fair comparisons and statistical power.",
        "2. Sequence-context perturbation panel",
        "For each variant, run at least one altered RNA context (mutant/non-binder control) while keeping Zn entities constant.",
        "3. Zn-entity perturbation sanity checks",
        "Run a limited set with modified ion setup to test whether deleterious/neutral trends are Zn-model dependent.",
        "4. Single-mutation neighborhood scans",
        "Run nearby substitutions around top signals (e.g., R160, C206-210 region) to map local static sensitivity gradients.",
        "",
        "## Priority 3: AF3-derived features for stronger classification",
        "",
        "1. Contact-map overlap to WT",
        "Compute residue-by-residue protein-RNA contact retention from static structures.",
        "2. Per-residue confidence drop signatures",
        "Track whether confidence loss localizes to zinc fingers, linker, or RNA-contact residues.",
        "3. Interface electrostatics proxy panel",
        "Add charge-change and aromaticity-change descriptors alongside AF3 confidence metrics.",
        "",
        "## Suggested execution order (24-hour window)",
        "",
        "1. Finish Priority 1 analyses first (no new runs required).",
        "2. Launch only C208A balancing runs next (fastest improvement in dataset quality).",
        "3. Use updated balanced dataset for final presentation plots and conclusions.",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Build AF3+MD presentation bundle")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output bundle directory (default: MD_GCP_Pipeline_Consolidated/presentation/compiled_YYYYMMDD)",
    )
    args = parser.parse_args()

    script_path = Path(__file__).resolve()
    consolidated_root = script_path.parents[2]
    workspace_root = consolidated_root.parent

    if args.output_dir is None:
        out_root = consolidated_root / "presentation" / f"compiled_{date.today().strftime('%Y%m%d')}"
    else:
        out_root = args.output_dir

    figs_dir = ensure_dir(out_root / "figures")
    tables_dir = ensure_dir(out_root / "tables")
    copied_figs_dir = ensure_dir(out_root / "figures" / "from_existing")

    # Load data
    df_models = load_per_model_pair_metrics(workspace_root)
    df_af3 = load_af3_variant_summary(workspace_root)
    md_qc, md_ca, md_rna = load_md_tables(consolidated_root)
    md_provenance = load_md_run_provenance(consolidated_root, md_qc)

    # Canonical AF3 table with forced class map
    df_af3 = df_af3[df_af3["variant"].isin(SELECTED_VARIANTS)].copy()
    df_af3["classification"] = df_af3["variant"].map(CLASS_MAP)

    # Composite score + class-separation stats
    df_score = build_af3_composite_score(df_af3)
    df_score_mut = df_score[df_score["classification"].isin(["neutral", "deleterious"])].copy()
    y = (df_score_mut["classification"] == "deleterious").astype(int).to_numpy()
    s = df_score_mut["af3_composite_badness"].to_numpy()
    auc = roc_auc_from_scores(y, s)
    auc_lo, auc_hi = auc_bootstrap_ci(y, s)
    auc_p = auc_permutation_pvalue(y, s)

    # Per-model summary tables
    agg = (
        df_models.groupby(["variant", "classification"], dropna=False)[
            ["pair_ipTM_protein_rna", "pair_pae_min_cross_mean", "iptm", "ptm", "ranking_score"]
        ]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    agg.columns = [
        "_".join([c for c in col if c]).rstrip("_") if isinstance(col, tuple) else col
        for col in agg.columns
    ]

    agg_m0 = (
        df_models[df_models["model_local_idx"] == 0]
        .groupby(["variant", "classification"], dropna=False)[
            ["pair_ipTM_protein_rna", "pair_pae_min_cross_mean", "iptm", "ptm", "ranking_score"]
        ]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    agg_m0.columns = [
        "_".join([c for c in col if c]).rstrip("_") if isinstance(col, tuple) else col
        for col in agg_m0.columns
    ]

    # Merge AF3 + MD
    md_qc_core = md_qc[[
        "variant",
        "classification",
        "prod_final_step",
        "prod_temp_mean_k",
        "prod_speed_mean_nsday",
        "qc_status",
    ]].copy()
    md_qc_core = md_qc_core.rename(columns={"classification": "classification_md"})

    bridge = (
        df_score[["variant", "classification", "af3_composite_badness", "pair_ipTM_protein_rna", "pair_pae_min_cross_mean", "pLDDT_min"]]
        .merge(md_ca[["variant", "ca_rmsd_a_final_vs_start"]], on="variant", how="left")
        .merge(md_rna[["variant", "p_rmsd_b_final_vs_start"]], on="variant", how="left")
        .merge(md_qc_core, on="variant", how="left")
        .merge(md_provenance[["variant", "run_source_group", "production_profile"]], on="variant", how="left")
    )

    # Save tables
    df_models.to_csv(tables_dir / "af3_per_model_pair_metrics_selected.csv", index=False)
    agg.to_csv(tables_dir / "af3_per_model_pair_summary_selected.csv", index=False)
    agg_m0.to_csv(tables_dir / "af3_model0_pair_summary_selected.csv", index=False)
    df_af3.to_csv(tables_dir / "af3_variant_summary_selected.csv", index=False)
    df_score.to_csv(tables_dir / "af3_variant_scores_selected.csv", index=False)
    df_score_mut.to_csv(tables_dir / "af3_deleterious_vs_neutral_selected.csv", index=False)
    bridge.to_csv(tables_dir / "combined_af3_md_selected.csv", index=False)
    md_provenance.to_csv(tables_dir / "md_run_provenance_selected.csv", index=False)

    # Save compact stats table
    stats_df = pd.DataFrame(
        [
            {
                "metric": "auc_deleterious_vs_neutral",
                "value": auc,
                "ci95_lo": auc_lo,
                "ci95_hi": auc_hi,
                "permutation_pvalue": auc_p,
                "n_variants": len(df_score_mut),
            }
        ]
    )
    stats_df.to_csv(tables_dir / "af3_classification_stats_selected.csv", index=False)

    # New figures
    plot_af3_variant_metric_bars(df_af3.rename(columns={"variant": "variant"}), figs_dir / "figure01_af3_pair_metrics_by_variant.png")
    plot_af3_model0_run_stability(df_models, figs_dir / "figure02_af3_model0_run_stability.png")
    plot_af3_class_separation(df_score_mut, figs_dir / "figure03_af3_class_separation.png")
    plot_af3_md_bridge(bridge, figs_dir / "figure04_af3_md_bridge.png")
    plot_md_qc_overview(md_qc.rename(columns={"classification": "classification_md"}), figs_dir / "figure05_md_qc_overview.png")

    # Copy key existing figures for continuity
    legacy = [
        workspace_root / "results" / "plots" / "composite_deviation.png",
        workspace_root / "results" / "plots" / "composite_deviation_all.png",
        workspace_root / "results" / "plots" / "composite_deviation_with_plddt.png",
        workspace_root / "results" / "plots" / "composite_deviation_with_plddt_all.png",
        workspace_root / "results" / "plots" / "pair_iptm_protein_rna_barplot_all.png",
        workspace_root / "results" / "plots" / "pae_cross_mean_barplot_all.png",
        workspace_root / "analysis" / "plots" / "rmsd_vs_contacts.png",
        workspace_root / "analysis" / "plots" / "interface_metric_by_protein_class.png",
    ]
    for src in legacy:
        if src.exists():
            shutil.copy2(src, copied_figs_dir / src.name)

    # Reports
    bundle_rel = str(out_root.relative_to(consolidated_root)) if out_root.is_relative_to(consolidated_root) else str(out_root)
    write_report(
        out_root / "REPORT.md",
        bundle_rel,
        auc,
        auc_lo,
        auc_hi,
        auc_p,
        df_score,
        md_qc,
        md_provenance,
    )
    write_af3_expansion_plan(out_root / "AF3_STATIC_EXPANSION_PLAN.md")

    # Small manifest
    manifest = pd.DataFrame(
        {
            "key": [
                "workspace_root",
                "consolidated_root",
                "output_dir",
                "selected_variants",
                "excluded_variants",
            ],
            "value": [
                str(workspace_root),
                str(consolidated_root),
                str(out_root),
                ",".join(SELECTED_VARIANTS),
                ",".join(EXCLUDED_VARIANTS),
            ],
        }
    )
    manifest.to_csv(out_root / "bundle_manifest.csv", index=False)

    print(f"Wrote bundle: {out_root}")
    print(f"  tables:  {tables_dir}")
    print(f"  figures: {figs_dir}")
    print(f"  report:  {out_root / 'REPORT.md'}")


if __name__ == "__main__":
    main()
