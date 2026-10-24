#!/usr/bin/env python3
"""
Cache the per-frame T156<->RNA H-bond series used by the convergence figure.

The trajectories (production.dcd) are too large for the repository, so this
script runs once locally and writes a small CSV that make_figures.py reads.
It reuses load_traj() and f156_hbond_series() from
expanse_md_analysis/diagnostic_f156_convergence.py, so the numbers are the
same ones behind diagnostic_f156_convergence.csv (1 ns sampling, <=200 ns).

Note: "f156" in that script is a historical misnomer. resSeq 6 of the
construct is Thr156 (the T of the RYKTEL motif), not a phenylalanine.

Output:
    figures/data/t156_rna_hbond_series_1ns.csv
        variant, class, source, t_ns, hbonds

Env: /Users/eric/micromamba/envs/omm/bin/python  (needs mdtraj)
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "expanse_md_analysis"))
import diagnostic_f156_convergence as diag  # noqa: E402

OUT = ROOT / "figures" / "data" / "t156_rna_hbond_series_1ns.csv"


def main() -> None:
    rows = []
    for src_name, base_dir in diag.SOURCES:
        for variant in diag.ALL_VARIANTS:
            traj = diag.load_traj(base_dir, variant)
            if traj is None:
                print(f"  {src_name:13s} {variant:6s} missing")
                continue
            series = diag.f156_hbond_series(traj)
            for i, value in enumerate(series):
                rows.append({
                    "variant": variant,
                    "class": diag.CLASS[variant],
                    "source": src_name,
                    "t_ns": round(i * diag.NS_PER_FRAME, 3),
                    "hbonds": value,
                })
            print(f"  {src_name:13s} {variant:6s} {len(series):4d} frames, "
                  f"mean {series.mean():.4f}", flush=True)
    pd.DataFrame(rows).to_csv(OUT, index=False)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
