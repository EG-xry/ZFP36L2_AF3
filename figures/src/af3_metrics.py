#!/usr/bin/env python3
"""
Per-model AlphaFold 3 confidence metrics for every construct in Data/.

Each Data/<variant>/<run>/ folder is one AlphaFold Server job (one seed,
five ranked models). For every model this records the same metrics that
analyze_best_models.py reports, read straight from the server outputs:

    pair_ipTM_protein_rna    chain_pair_iptm[protein][rna]
    pae_cross_mean           mean of chain_pair_pae_min[protein][rna] and [rna][protein]
    pLDDT                    mean Ca pLDDT of the protein chain (CIF B-factor column)
    pLDDT_min                lowest per-residue Ca pLDDT (pLDDT_min_res = which residue)
    ranking_score, rna_seq, mutation (vs WT, construct + ZFP36L2 numbering)

Output: figures/data/af3_per_model.csv
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "Data"
OUT = ROOT / "figures" / "data" / "af3_per_model.csv"

CONSTRUCT_START = 151  # construct residue 1 = ZFP36L2 residue 151 (NP_008818.3)
WT_SEQ = "STRYKTELCRPFEESGTCKYGEKCQFAHGFHELRSLTRHPKYKTELCRTFHTIGFCPYGPRCHFIHNADE"
RNA_COMPS = {"A", "U", "G", "C"}
AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def read_atom_site(cif_path: Path) -> pd.DataFrame:
    """Minimal mmCIF _atom_site reader (AF3 server files are one loop, no quoting)."""
    cols, rows = [], []
    for line in cif_path.read_text().splitlines():
        if line.startswith("_atom_site."):
            cols.append(line.split(".", 1)[1].strip())
        elif cols and line.startswith(("ATOM", "HETATM")):
            rows.append(line.split())
    return pd.DataFrame(rows, columns=cols)


def model_metrics(cif_path: Path, conf_path: Path) -> dict:
    atoms = read_atom_site(cif_path)
    chains = list(dict.fromkeys(atoms["label_asym_id"]))
    comp_by_chain = atoms.groupby("label_asym_id")["label_comp_id"].agg(set)
    p = next(i for i, c in enumerate(chains) if comp_by_chain[c] & set(AA3))
    r = next(i for i, c in enumerate(chains) if comp_by_chain[c] <= RNA_COMPS)

    prot = atoms[atoms["label_asym_id"] == chains[p]]
    ca = prot[prot["label_atom_id"] == "CA"]
    ca_plddt = ca["B_iso_or_equiv"].astype(float).to_numpy()
    seq = "".join(AA3[x] for x in ca["label_comp_id"])
    rna = atoms[atoms["label_asym_id"] == chains[r]]
    rna_seq = "".join(rna.drop_duplicates("label_seq_id")["label_comp_id"])

    conf = json.loads(conf_path.read_text())
    iptm = conf["chain_pair_iptm"][p][r]
    pae = conf["chain_pair_pae_min"]
    return {
        "pair_ipTM_protein_rna": iptm,
        "pae_cross_mean": 0.5 * (pae[p][r] + pae[r][p]),
        "pLDDT": float(np.mean(ca_plddt)),
        "pLDDT_min": float(np.min(ca_plddt)),
        "pLDDT_min_res": int(np.argmin(ca_plddt)) + CONSTRUCT_START,
        "ranking_score": conf["ranking_score"],
        "protein_seq": seq,
        "rna_seq": rna_seq,
    }


def describe_mutation(seq: str) -> str:
    diffs = [f"{w}{i + 1}{m} ({w}{i + CONSTRUCT_START}{m})"
             for i, (w, m) in enumerate(zip(WT_SEQ, seq)) if w != m]
    return "; ".join(diffs) if diffs else "none"


def main() -> None:
    rows = []
    for variant_dir in sorted(p for p in DATA.iterdir() if p.is_dir()):
        runs = sorted(p for p in variant_dir.iterdir() if p.is_dir())
        for run_idx, run_dir in enumerate(runs, start=1):
            for model_idx in range(5):
                cif = next(run_dir.glob(f"*_model_{model_idx}.cif"))
                conf = next(run_dir.glob(f"*_summary_confidences_{model_idx}.json"))
                m = model_metrics(cif, conf)
                m.update(variant=variant_dir.name, run=run_idx, model=model_idx,
                         mutation=describe_mutation(m.pop("protein_seq")))
                rows.append(m)
    df = pd.DataFrame(rows)[["variant", "run", "model", "mutation", "rna_seq",
                             "pair_ipTM_protein_rna", "pae_cross_mean", "pLDDT",
                             "pLDDT_min", "pLDDT_min_res", "ranking_score"]]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False, float_format="%.4f")
    print(f"wrote {OUT} ({len(df)} models, {df.groupby(['variant', 'run']).ngroups} jobs)")


if __name__ == "__main__":
    main()
