#!/usr/bin/env python3
"""
Parse md.log files from all expanse_md_results variants and produce:
  - per-variant time series CSV (production phase only)
  - summary table: ns_completed, mean T, mean E, drift
  - sanity check (T ~310K, no NaN, energy stable in last 25%)

The md.log files are written by OpenMM StateDataReporter with 3 different
"flavors" because the script appends to the same file across 3 stages:
  Stage 1 (NVT):   step, Epot, Ekin, T, speed   (5 cols, no header row mid-file)
  Stage 2 (NPT-eq):step, Epot, T, V, density, speed   (6 cols)
  Stage 3 (prod):  progress%, step, time_ps, Epot, Ekin, Etot, T, V, density,
                   speed, remaining           (11 cols including progress + remaining)

We detect production rows by the leading "%" sign that OpenMM writes when
`progress=True` is enabled.
"""
from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path
from statistics import mean, pstdev

import pandas as pd

EXPANSE_DIR = Path(__file__).resolve().parent.parent / "expanse_md_results"
OUT_DIR = Path(__file__).resolve().parent
OUT_DIR.mkdir(parents=True, exist_ok=True)

VARIANTS = ["WT", "C206F", "C212A", "D210E", "G204D",
            "R160W", "R184C", "S185N", "T187S", "Y154H"]

CLASS_MAP = {
    "WT": "WT",
    "C206F": "deleterious",
    "C212A": "zn_knockout_control",  # was mislabeled C208A; mutation is at C212
    "G204D": "deleterious",
    "R160W": "deleterious",
    "R184C": "deleterious",
    "Y154H": "deleterious",
    "D210E": "neutral",  # NOTE: actual mutation in this run is D219E
    "S185N": "neutral",
    "T187S": "neutral",
}

PROGRESS_RE = re.compile(r"^\s*\d+(?:\.\d+)?%\s*\t")


def parse_production_rows(log_path: Path):
    """Yield production rows. Each row is a list of floats:
    [progress, step, time_ps, Epot, Ekin, Etot, T, V, density, speed]
    The wall remaining (last col) is HH:MM:SS or D:HH:MM:SS so we drop it.
    """
    with log_path.open() as f:
        for line in f:
            if not PROGRESS_RE.match(line):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            try:
                progress = float(parts[0].rstrip("%"))
                step = int(float(parts[1]))
                time_ps = float(parts[2])
                Epot = float(parts[3])
                Ekin = float(parts[4])
                Etot = float(parts[5])
                T = float(parts[6])
                V = float(parts[7])
                density = float(parts[8])
                speed = float(parts[9])
            except (ValueError, IndexError):
                continue
            yield (progress, step, time_ps, Epot, Ekin, Etot, T, V, density, speed)


def read_progress(variant: str) -> dict:
    p = EXPANSE_DIR / variant / "progress.txt"
    info = {}
    if not p.exists():
        return info
    for line in p.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            info[k.strip()] = v.strip()
    return info


def quick_stats(values):
    if not values:
        return float("nan"), float("nan")
    return mean(values), (pstdev(values) if len(values) > 1 else 0.0)


def main():
    summary_rows = []
    for variant in VARIANTS:
        log = EXPANSE_DIR / variant / "md.log"
        prog = read_progress(variant)
        if not log.exists():
            print(f"[{variant}] no log; skipping")
            continue

        print(f"[{variant}] parsing {log} ...")
        rows = list(parse_production_rows(log))
        if not rows:
            print(f"  no production rows found")
            continue

        # Build dataframe
        df = pd.DataFrame(
            rows,
            columns=["progress_pct", "step", "time_ps", "Epot_kJmol",
                     "Ekin_kJmol", "Etot_kJmol", "T_K", "V_nm3",
                     "density_g_per_mL", "speed_ns_per_day"],
        )
        # ns of production
        df["time_ns"] = df["time_ps"] / 1000.0
        # Save per-variant CSV
        out_csv = OUT_DIR / f"{variant}_production.csv"
        df.to_csv(out_csv, index=False)

        # Sanity: T mean, T std (overall and last 25%)
        T_mean, T_std = quick_stats(df["T_K"].tolist())
        last_q = df.iloc[int(0.75 * len(df)):]
        T_last_mean, T_last_std = quick_stats(last_q["T_K"].tolist())
        E_first = df.iloc[: max(1, len(df) // 20)]["Epot_kJmol"].mean()
        E_last = last_q["Epot_kJmol"].mean()
        E_drift_per_ns = (
            (E_last - E_first) / max(1e-9, df["time_ns"].iloc[-1] - df["time_ns"].iloc[0])
        )
        nan_count = int(df.isna().sum().sum())

        # Time span
        ns_completed = float(prog.get("ns_completed", df["time_ns"].iloc[-1] - df["time_ns"].iloc[0]))

        summary_rows.append({
            "variant": variant,
            "class": CLASS_MAP.get(variant, "?"),
            "status": prog.get("status", "?"),
            "target_ns": prog.get("target_ns", "?"),
            "ns_completed_from_progress": prog.get("ns_completed", ""),
            "rows_in_log": len(df),
            "time_ns_first_logged": float(df["time_ns"].iloc[0]),
            "time_ns_last_logged": float(df["time_ns"].iloc[-1]),
            "production_window_ns": float(df["time_ns"].iloc[-1] - df["time_ns"].iloc[0]),
            "T_mean_K": round(T_mean, 3),
            "T_std_K": round(T_std, 3),
            "T_last25pct_mean_K": round(T_last_mean, 3),
            "T_last25pct_std_K": round(T_last_std, 3),
            "Epot_first1pct_mean_kJmol": round(E_first, 1),
            "Epot_last25pct_mean_kJmol": round(E_last, 1),
            "Epot_drift_kJmol_per_ns": round(E_drift_per_ns, 2),
            "density_last25pct_mean_g_per_mL": round(last_q["density_g_per_mL"].mean(), 4),
            "speed_mean_ns_per_day": round(df["speed_ns_per_day"].mean(), 1),
            "nan_count": nan_count,
            "production_hours": prog.get("production_hours", ""),
        })

    summary = pd.DataFrame(summary_rows)
    summary = summary.sort_values(
        ["class", "variant"],
        key=lambda c: c if c.name != "class" else c.map(
            {"WT": 0, "deleterious": 1, "deleterious_in_set": 1, "neutral": 2}
        ),
    )

    out_summary = OUT_DIR / "md_run_summary.csv"
    summary.to_csv(out_summary, index=False)
    # also a pretty markdown copy (avoid tabulate dep)
    out_md = OUT_DIR / "md_run_summary.md"
    cols = list(summary.columns)
    with out_md.open("w") as f:
        f.write("| " + " | ".join(cols) + " |\n")
        f.write("| " + " | ".join("---" for _ in cols) + " |\n")
        for _, row in summary.iterrows():
            f.write("| " + " | ".join(str(v) for v in row) + " |\n")

    print("\n=== Summary ===")
    print(summary.to_string(index=False))
    print(f"\nWrote: {out_summary}")
    print(f"Wrote: {out_md}")


if __name__ == "__main__":
    main()
