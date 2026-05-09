#!/usr/bin/env python3
"""Compute per-variant structural drift from cloud minimal-MD outputs.

Outputs:
- manifests/cloud_variant_structural_drift.csv      (protein CA RMSD)
- manifests/cloud_variant_structural_drift_rna.csv  (RNA P-atom RMSD)
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np


AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
RNA1 = {"A", "U", "G", "C", "I"}


def parse_atoms(pdb_path: Path) -> List[Dict[str, object]]:
    atoms: List[Dict[str, object]] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if len(line) < 54:
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip() or "_"
        resid = line[22:26].strip()
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        atoms.append(
            {
                "name": name,
                "resname": resname,
                "chain": chain,
                "resid": resid,
                "xyz": np.array([x, y, z], dtype=float),
            }
        )
    return atoms


def kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    if p.shape != q.shape:
        raise ValueError("Coordinate arrays must have same shape")
    if p.shape[0] < 3:
        return float("nan")

    p_cent = p.mean(axis=0)
    q_cent = q.mean(axis=0)
    p0 = p - p_cent
    q0 = q - q_cent

    c = np.dot(q0.T, p0)
    v, s, w_t = np.linalg.svd(c)
    d = np.sign(np.linalg.det(np.dot(v, w_t)))
    d_mat = np.diag([1.0, 1.0, d])
    rot = np.dot(np.dot(v, d_mat), w_t)

    q_aligned = np.dot(q0, rot)
    diff = p0 - q_aligned
    return float(np.sqrt(np.mean(np.sum(diff * diff, axis=1))))


def _key(atom: Dict[str, object]) -> Tuple[str, str, str, str]:
    return (
        str(atom["chain"]),
        str(atom["resid"]),
        str(atom["resname"]),
        str(atom["name"]),
    )


def _select_atoms(atoms: Iterable[Dict[str, object]], mode: str) -> Dict[Tuple[str, str, str, str], np.ndarray]:
    sel: Dict[Tuple[str, str, str, str], np.ndarray] = {}
    for a in atoms:
        name = str(a["name"])
        resname = str(a["resname"]).upper()
        if mode == "protein_ca":
            if name != "CA" or resname not in AA3:
                continue
        else:
            raise ValueError(f"Unknown mode: {mode}")
        sel[_key(a)] = np.asarray(a["xyz"], dtype=float)
    return sel


def _select_rna_p_ordered(atoms: Iterable[Dict[str, object]]) -> List[Tuple[str, np.ndarray]]:
    out: List[Tuple[str, np.ndarray]] = []
    for a in atoms:
        name = str(a["name"])
        resname = str(a["resname"]).upper()
        if name == "P" and resname in RNA1:
            out.append((resname, np.asarray(a["xyz"], dtype=float)))
    return out


def _lcs_pairs(a: List[str], b: List[str]) -> List[Tuple[int, int]]:
    """Longest-common-subsequence index pairs for short RNA base strings."""
    n, m = len(a), len(b)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(n - 1, -1, -1):
        for j in range(m - 1, -1, -1):
            if a[i] == b[j]:
                dp[i][j] = 1 + dp[i + 1][j + 1]
            else:
                dp[i][j] = max(dp[i + 1][j], dp[i][j + 1])

    pairs: List[Tuple[int, int]] = []
    i = j = 0
    while i < n and j < m:
        if a[i] == b[j]:
            pairs.append((i, j))
            i += 1
            j += 1
        elif dp[i + 1][j] >= dp[i][j + 1]:
            i += 1
        else:
            j += 1
    return pairs


def find_final_pdb(variant_dir: Path) -> Path | None:
    if not variant_dir.exists():
        return None
    for child in sorted(variant_dir.iterdir()):
        if child.is_dir():
            f = child / "final_frame.pdb"
            if f.exists():
                return f
    return None


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "variant_run_manifest.csv",
        help="Variant manifest CSV path.",
    )
    parser.add_argument(
        "--inputs-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "runs" / "variants",
        help="Directory containing input_af.pdb under variant subfolders.",
    )
    parser.add_argument(
        "--outputs-root",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "runs" / "cloud_outputs" / "variants",
        help="Directory containing cloud output variant subfolders.",
    )
    parser.add_argument(
        "--out-ca-csv",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "cloud_variant_structural_drift.csv",
        help="Output CSV path for protein CA RMSD.",
    )
    parser.add_argument(
        "--out-rna-csv",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "cloud_variant_structural_drift_rna.csv",
        help="Output CSV path for RNA P RMSD.",
    )
    args = parser.parse_args()

    if not args.manifest.exists():
        raise SystemExit(f"Manifest not found: {args.manifest}")

    rows_ca: List[Dict[str, object]] = []
    rows_rna: List[Dict[str, object]] = []

    with args.manifest.open("r", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            variant = row["variant"]
            classification = row["classification"]
            start_pdb = args.inputs_root / variant / "input_af.pdb"
            final_pdb = find_final_pdb(args.outputs_root / variant)

            if not start_pdb.exists() or final_pdb is None or not final_pdb.exists():
                rows_ca.append(
                    {
                        "variant": variant,
                        "classification": classification,
                        "n_ca": 0,
                        "ca_rmsd_a_final_vs_start": "",
                        "status": "missing_files",
                    }
                )
                rows_rna.append(
                    {
                        "variant": variant,
                        "classification": classification,
                        "n_p": 0,
                        "p_rmsd_b_final_vs_start": "",
                        "status": "missing_files",
                    }
                )
                continue

            try:
                start_atoms = parse_atoms(start_pdb)
                final_atoms = parse_atoms(final_pdb)

                s_ca = _select_atoms(start_atoms, "protein_ca")
                f_ca = _select_atoms(final_atoms, "protein_ca")
                common_ca = sorted(set(s_ca).intersection(f_ca))
                if len(common_ca) >= 3:
                    p = np.vstack([s_ca[k] for k in common_ca])
                    q = np.vstack([f_ca[k] for k in common_ca])
                    rmsd_ca = kabsch_rmsd(p, q)
                    ca_status = "ok"
                else:
                    rmsd_ca = float("nan")
                    ca_status = "insufficient_common_atoms"

                s_rna = _select_rna_p_ordered(start_atoms)
                f_rna = _select_rna_p_ordered(final_atoms)
                seq_s = [x[0] for x in s_rna]
                seq_f = [x[0] for x in f_rna]
                pairs = _lcs_pairs(seq_s, seq_f)
                if len(pairs) >= 3:
                    p = np.vstack([s_rna[i][1] for i, _ in pairs])
                    q = np.vstack([f_rna[j][1] for _, j in pairs])
                    rmsd_rna = kabsch_rmsd(p, q)
                    rna_status = "ok"
                else:
                    rmsd_rna = float("nan")
                    rna_status = "insufficient_common_atoms"

                rows_ca.append(
                    {
                        "variant": variant,
                        "classification": classification,
                        "n_ca": len(common_ca),
                        "ca_rmsd_a_final_vs_start": "" if math.isnan(rmsd_ca) else f"{rmsd_ca:.4f}",
                        "status": ca_status,
                    }
                )
                rows_rna.append(
                    {
                        "variant": variant,
                        "classification": classification,
                        "n_p": len(pairs),
                        "p_rmsd_b_final_vs_start": "" if math.isnan(rmsd_rna) else f"{rmsd_rna:.4f}",
                        "status": rna_status,
                    }
                )
            except Exception as exc:  # pragma: no cover
                rows_ca.append(
                    {
                        "variant": variant,
                        "classification": classification,
                        "n_ca": 0,
                        "ca_rmsd_a_final_vs_start": "",
                        "status": f"error:{type(exc).__name__}",
                    }
                )
                rows_rna.append(
                    {
                        "variant": variant,
                        "classification": classification,
                        "n_p": 0,
                        "p_rmsd_b_final_vs_start": "",
                        "status": f"error:{type(exc).__name__}",
                    }
                )

    args.out_ca_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_ca_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["variant", "classification", "n_ca", "ca_rmsd_a_final_vs_start", "status"],
        )
        w.writeheader()
        w.writerows(rows_ca)

    with args.out_rna_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["variant", "classification", "n_p", "p_rmsd_b_final_vs_start", "status"],
        )
        w.writeheader()
        w.writerows(rows_rna)

    print(f"Wrote protein drift CSV: {args.out_ca_csv}")
    print(f"Wrote RNA drift CSV: {args.out_rna_csv}")


if __name__ == "__main__":
    main()
