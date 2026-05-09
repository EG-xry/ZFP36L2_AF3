#!/usr/bin/env python3
"""Stage per-variant input folders from a manifest for cloud execution."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from pathlib import Path


def stage_one(row: dict[str, str], out_root: Path, mode: str) -> Path:
    variant = row["variant"]
    source_pdb = Path(row["source_pdb"])
    variant_dir = out_root / variant
    variant_dir.mkdir(parents=True, exist_ok=True)

    target_pdb = variant_dir / "input_af.pdb"
    if target_pdb.exists() or target_pdb.is_symlink():
        target_pdb.unlink()

    if mode == "symlink":
        target_pdb.symlink_to(source_pdb)
    else:
        shutil.copy2(source_pdb, target_pdb)

    metadata = {
        "variant": variant,
        "classification": row["classification"],
        "source_pdb": str(source_pdb),
        "source_tag": row["source_tag"],
        "temperature_k": float(row["temperature_k"]),
        "ph": float(row["ph"]),
        "ionic_strength_m": float(row["ionic_strength_m"]),
        "status": row["status"],
    }
    (variant_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    return variant_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "variant_run_manifest.csv",
        help="CSV manifest path.",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "runs" / "variants",
        help="Output directory for staged per-variant folders.",
    )
    parser.add_argument(
        "--mode",
        choices=["copy", "symlink"],
        default="copy",
        help="Whether to copy or symlink source PDB files.",
    )
    args = parser.parse_args()

    if not args.manifest.exists():
        raise SystemExit(f"Manifest not found: {args.manifest}")

    with args.manifest.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    if not rows:
        raise SystemExit("Manifest has no rows.")

    staged: list[Path] = []
    for row in rows:
        src = Path(row["source_pdb"])
        if not src.exists():
            raise SystemExit(f"Missing source PDB: {src}")
        staged.append(stage_one(row, args.out_root, args.mode))

    print(f"Staged {len(staged)} variants under: {args.out_root}")
    for d in staged:
        print(f"- {d}")


if __name__ == "__main__":
    main()
