#!/usr/bin/env python3
"""Build Google Cloud Batch taskEnvironments from the variant manifest."""
from __future__ import annotations

import argparse
import csv
import json
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
        "--out-json",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "batch_task_environments.json",
        help="Output JSON file containing taskEnvironments array.",
    )
    parser.add_argument(
        "--out-variant-order",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "variant_order.txt",
        help="Output text file listing variant order by task index.",
    )
    parser.add_argument(
        "--input-root-prefix",
        default="variants",
        help="Path prefix used by the job runtime to locate staged variant inputs.",
    )
    args = parser.parse_args()

    if not args.manifest.exists():
        raise SystemExit(f"Manifest not found: {args.manifest}")

    with args.manifest.open("r", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))

    if not rows:
        raise SystemExit("Manifest has no rows.")

    task_envs: list[dict[str, dict[str, str]]] = []
    order_lines: list[str] = []
    for idx, row in enumerate(rows):
        variant = row["variant"]
        rel_input = f"{args.input_root_prefix}/{variant}/input_af.pdb"
        task_envs.append(
            {
                "variables": {
                    "VARIANT": variant,
                    "CLASSIFICATION": row["classification"],
                    "REL_INPUT_PDB": rel_input,
                    "TEMPERATURE_K": row["temperature_k"],
                    "PH": row["ph"],
                    "IONIC_STRENGTH_M": row["ionic_strength_m"],
                }
            }
        )
        order_lines.append(f"{idx}\t{variant}\t{row['classification']}\t{rel_input}")

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(task_envs, indent=2), encoding="utf-8")

    args.out_variant_order.parent.mkdir(parents=True, exist_ok=True)
    args.out_variant_order.write_text(
        "task_index\tvariant\tclassification\trel_input_pdb\n"
        + "\n".join(order_lines)
        + "\n",
        encoding="utf-8",
    )

    print(f"Wrote {len(task_envs)} task environment entries to: {args.out_json}")
    print(f"Wrote variant order table to: {args.out_variant_order}")


if __name__ == "__main__":
    main()
