#!/usr/bin/env python3
"""Basic QC/summary analysis for cloud minimal-MD variant outputs."""
from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path


def _to_float(value: str) -> float | None:
    value = value.strip()
    if not value or value == "--":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _find_md_subdir(variant_dir: Path) -> Path | None:
    if not variant_dir.exists():
        return None
    for child in sorted(variant_dir.iterdir()):
        if child.is_dir() and (child / "md.log").exists():
            return child
    return None


def _parse_md_log(md_log: Path) -> dict[str, float | int | str | None]:
    lines = [ln for ln in md_log.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if len(lines) < 2:
        return {
            "n_rows_total": 0,
            "n_stages_detected": 0,
            "prod_rows": 0,
            "prod_final_step": None,
            "prod_temp_mean_k": None,
            "prod_temp_sd_k": None,
            "prod_potential_mean_kjmol": None,
            "prod_potential_final_kjmol": None,
            "prod_speed_mean_nsday": None,
        }

    header = [h.strip().strip('"').lstrip("#") for h in lines[0].split("\t")]
    raw_rows = [dict(zip(header, row.split("\t"))) for row in lines[1:]]

    steps: list[int] = []
    for row in raw_rows:
        step_val = row.get("Step", "").strip()
        steps.append(int(float(step_val)) if step_val else 0)

    stage_indices: list[int] = []
    stage = 0
    last_step = -1
    for step in steps:
        if last_step >= 0 and step < last_step:
            stage += 1
        stage_indices.append(stage)
        last_step = step

    n_stages = stage + 1 if stage_indices else 0
    prod_stage = stage if stage_indices else -1
    prod_rows = [row for row, sidx in zip(raw_rows, stage_indices) if sidx == prod_stage]

    prod_steps = [int(float(r["Step"])) for r in prod_rows] if prod_rows else []
    prod_temps = [_to_float(r.get("Temperature (K)", "")) for r in prod_rows]
    prod_temps = [v for v in prod_temps if v is not None]
    prod_pots = [_to_float(r.get("Potential Energy (kJ/mole)", "")) for r in prod_rows]
    prod_pots = [v for v in prod_pots if v is not None]
    prod_speeds = [_to_float(r.get("Speed (ns/day)", "")) for r in prod_rows]
    prod_speeds = [v for v in prod_speeds if v is not None]

    return {
        "n_rows_total": len(raw_rows),
        "n_stages_detected": n_stages,
        "prod_rows": len(prod_rows),
        "prod_final_step": prod_steps[-1] if prod_steps else None,
        "prod_temp_mean_k": statistics.mean(prod_temps) if prod_temps else None,
        "prod_temp_sd_k": statistics.pstdev(prod_temps) if len(prod_temps) > 1 else None,
        "prod_potential_mean_kjmol": statistics.mean(prod_pots) if prod_pots else None,
        "prod_potential_final_kjmol": prod_pots[-1] if prod_pots else None,
        "prod_speed_mean_nsday": statistics.mean(prod_speeds) if prod_speeds else None,
    }


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
        default=Path(__file__).resolve().parents[2] / "runs" / "cloud_outputs" / "variants",
        help="Local sync root of gs://.../outputs/variants.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "cloud_variant_analysis.csv",
        help="Output CSV path.",
    )
    args = parser.parse_args()

    if not args.manifest.exists():
        raise SystemExit(f"Manifest not found: {args.manifest}")

    with args.manifest.open("r", encoding="utf-8") as fh:
        manifest_rows = list(csv.DictReader(fh))

    out_rows: list[dict[str, str]] = []
    for mrow in manifest_rows:
        variant = mrow["variant"]
        vdir = args.outputs_root / variant
        md_subdir = _find_md_subdir(vdir)
        md_log = md_subdir / "md.log" if md_subdir else None
        dcd = md_subdir / "production.dcd" if md_subdir else None
        final_pdb = md_subdir / "final_frame.pdb" if md_subdir else None
        meta = vdir / "task_metadata.env"

        stats = (
            _parse_md_log(md_log)
            if md_log is not None and md_log.exists()
            else {
                "n_rows_total": 0,
                "n_stages_detected": 0,
                "prod_rows": 0,
                "prod_final_step": None,
                "prod_temp_mean_k": None,
                "prod_temp_sd_k": None,
                "prod_potential_mean_kjmol": None,
                "prod_potential_final_kjmol": None,
                "prod_speed_mean_nsday": None,
            }
        )

        complete = bool(
            md_log
            and md_log.exists()
            and dcd
            and dcd.exists()
            and final_pdb
            and final_pdb.exists()
            and meta.exists()
            and (stats.get("prod_final_step") or 0) >= 25000
        )

        out_rows.append(
            {
                "variant": variant,
                "classification": mrow["classification"],
                "has_output_dir": str(vdir.exists()).lower(),
                "has_md_log": str(bool(md_log and md_log.exists())).lower(),
                "has_production_dcd": str(bool(dcd and dcd.exists())).lower(),
                "has_final_frame_pdb": str(bool(final_pdb and final_pdb.exists())).lower(),
                "has_task_metadata": str(meta.exists()).lower(),
                "production_dcd_bytes": str(dcd.stat().st_size if dcd and dcd.exists() else 0),
                "n_rows_total": str(stats["n_rows_total"]),
                "n_stages_detected": str(stats["n_stages_detected"]),
                "prod_rows": str(stats["prod_rows"]),
                "prod_final_step": "" if stats["prod_final_step"] is None else str(stats["prod_final_step"]),
                "prod_temp_mean_k": "" if stats["prod_temp_mean_k"] is None else f"{stats['prod_temp_mean_k']:.3f}",
                "prod_temp_sd_k": "" if stats["prod_temp_sd_k"] is None else f"{stats['prod_temp_sd_k']:.3f}",
                "prod_potential_mean_kjmol": ""
                if stats["prod_potential_mean_kjmol"] is None
                else f"{stats['prod_potential_mean_kjmol']:.3f}",
                "prod_potential_final_kjmol": ""
                if stats["prod_potential_final_kjmol"] is None
                else f"{stats['prod_potential_final_kjmol']:.3f}",
                "prod_speed_mean_nsday": ""
                if stats["prod_speed_mean_nsday"] is None
                else f"{stats['prod_speed_mean_nsday']:.3f}",
                "qc_status": "complete_minimal_md" if complete else "incomplete",
            }
        )

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "variant",
        "classification",
        "has_output_dir",
        "has_md_log",
        "has_production_dcd",
        "has_final_frame_pdb",
        "has_task_metadata",
        "production_dcd_bytes",
        "n_rows_total",
        "n_stages_detected",
        "prod_rows",
        "prod_final_step",
        "prod_temp_mean_k",
        "prod_temp_sd_k",
        "prod_potential_mean_kjmol",
        "prod_potential_final_kjmol",
        "prod_speed_mean_nsday",
        "qc_status",
    ]
    with args.out_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    n_complete = sum(1 for row in out_rows if row["qc_status"] == "complete_minimal_md")
    print(f"Wrote analysis summary: {args.out_csv}")
    print(f"Completed variants: {n_complete}/{len(out_rows)}")


if __name__ == "__main__":
    main()

