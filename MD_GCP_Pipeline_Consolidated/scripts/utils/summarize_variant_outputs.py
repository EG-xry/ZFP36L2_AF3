#!/usr/bin/env python3
"""Summarize per-variant output completeness into a CSV table."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "variant_run_manifest.csv",
        help="Variant manifest CSV path.",
    )
    parser.add_argument(
        "--outputs-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "runs" / "outputs",
        help="Directory containing one output folder per variant.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "variant_output_qc.csv",
        help="Output CSV summary path.",
    )
    args = parser.parse_args()

    if not args.manifest.exists():
        raise SystemExit(f"Manifest not found: {args.manifest}")

    with args.manifest.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    summary: list[dict[str, str]] = []
    for row in rows:
        variant = row["variant"]
        vdir = args.outputs_root / variant
        stats_files = list(vdir.glob("stats*.csv"))
        dcd_files = list(vdir.glob("*.dcd"))
        colvar_files = list(vdir.glob("COLVAR*"))
        hills_files = list(vdir.glob("HILLS*"))
        log_files = (
            list(vdir.glob("run.log"))
            + list(vdir.glob("simulation*.log"))
            + list(vdir.glob("*.out"))
            + list(vdir.glob("*.err"))
        )

        complete = all(
            [
                len(stats_files) > 0,
                len(dcd_files) > 0,
                len(colvar_files) > 0,
            ]
        )

        summary.append(
            {
                "variant": variant,
                "output_dir": str(vdir),
                "has_output_dir": str(vdir.exists()).lower(),
                "n_stats_csv": str(len(stats_files)),
                "n_dcd": str(len(dcd_files)),
                "n_colvar": str(len(colvar_files)),
                "n_hills": str(len(hills_files)),
                "n_logs": str(len(log_files)),
                "qc_status": "complete_minimum" if complete else "incomplete",
            }
        )

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "variant",
        "output_dir",
        "has_output_dir",
        "n_stats_csv",
        "n_dcd",
        "n_colvar",
        "n_hills",
        "n_logs",
        "qc_status",
    ]
    with args.out_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(summary)

    print(f"Wrote output QC summary for {len(summary)} variants: {args.out_csv}")


if __name__ == "__main__":
    main()
