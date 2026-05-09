#!/usr/bin/env python3
"""Build a one-run-per-variant manifest from run1/model0 PDB files."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

CORE_TO_VARIANT = {
    "wt_nonamer_wt_2zn": "WT",
    "y154h_nonamer_wt_2zn": "Y154H",
    "r160w_nonamer_wt_2zn": "R160W",
    "r160e_deleterious": "R160E",
    "r184c_nonamer_wt_2zn": "R184C",
    "s185n_nonamer_wt_2zn": "S185N",
    "t187s_nonamer_wt_2zn": "T187S",
    "g204d_nonamer_wt_2zn": "G204D",
    "c206f_nonamer_wt_2zn": "C206F",
    "deleterious_c208a": "C208A",
    "d210e_nonamer_wt_2zn": "D210E",
    "deleterious_h212a": "H212A",
}

VARIANT_CLASS = {
    "WT": "WT",
    "Y154H": "deleterious",
    "R160W": "deleterious",
    "R160E": "deleterious",
    "R184C": "deleterious",
    "S185N": "neutral",
    "T187S": "neutral",
    "G204D": "deleterious",
    "C206F": "deleterious",
    "C208A": "deleterious",
    "D210E": "neutral",
    "H212A": "deleterious",
}

VARIANT_ORDER = [
    "WT",
    "Y154H",
    "R160W",
    "R160E",
    "R184C",
    "S185N",
    "T187S",
    "G204D",
    "C206F",
    "C208A",
    "D210E",
    "H212A",
]

DEFAULT_EXCLUDED_VARIANTS = {"R160E", "H212A"}


def core_from_name(file_name: str) -> str:
    prefix = "fold_zfp36l2_"
    suffix = "_model_0_run1.pdb"
    if not file_name.startswith(prefix) or not file_name.endswith(suffix):
        raise ValueError(f"Unexpected file name format: {file_name}")
    return file_name[len(prefix) : -len(suffix)]


def build_rows(input_dir: Path, excluded_variants: set[str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for pdb in sorted(input_dir.glob("fold_zfp36l2*_model_0_run1.pdb")):
        core = core_from_name(pdb.name)
        variant = CORE_TO_VARIANT.get(core)
        if not variant:
            continue
        if variant in excluded_variants:
            continue
        rows.append(
            {
                "variant": variant,
                "classification": VARIANT_CLASS[variant],
                "source_pdb": str(pdb.resolve()),
                "source_tag": "run1_model0",
                "temperature_k": "298",
                "ph": "7.4",
                "ionic_strength_m": "0.2",
                "status": "ready_for_single_md_test",
            }
        )

    order_index = {variant: i for i, variant in enumerate(VARIANT_ORDER)}
    rows.sort(key=lambda r: order_index.get(r["variant"], 999))
    return rows


def write_csv(rows: list[dict[str, str]], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "variant",
        "classification",
        "source_pdb",
        "source_tag",
        "temperature_k",
        "ph",
        "ionic_strength_m",
        "status",
    ]
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "inputs" / "variant_models_run1_model0",
        help="Directory containing fold_zfp36l2*_model_0_run1.pdb files.",
    )
    parser.add_argument(
        "--out-csv",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "variant_run_manifest.csv",
        help="Output manifest CSV path.",
    )
    parser.add_argument(
        "--exclude",
        action="append",
        default=None,
        help=(
            "Variant to exclude from the manifest. "
            "Default excludes are R160E and H212A."
        ),
    )
    args = parser.parse_args()

    if not args.input_dir.exists():
        raise SystemExit(f"Input directory not found: {args.input_dir}")

    excluded = set(args.exclude) if args.exclude else set(DEFAULT_EXCLUDED_VARIANTS)
    rows = build_rows(args.input_dir, excluded)
    if not rows:
        raise SystemExit(f"No run1/model0 variant PDB files found in {args.input_dir}")

    write_csv(rows, args.out_csv)
    if excluded:
        print(f"Excluded variants: {', '.join(sorted(excluded))}")
    print(f"Wrote manifest with {len(rows)} entries: {args.out_csv}")


if __name__ == "__main__":
    main()
