#!/usr/bin/env python3
"""
Compute structural + AF3 confidence interface metrics for ZFP36L2 variants.

Canonical outputs (always produced):
  - analysis/interface_metrics_models.csv / .md
      Per-model table (all models, zinc + non-zinc)
  - analysis/interface_metrics_avg.csv / .md
      Averaged over models per complex (zinc + non-zinc)
  - analysis/interface_metrics_best.csv / .md
      Best model per complex (zinc + non-zinc)
  - analysis/interface_metrics_zinc_avg.csv / .md
      Averaged over models per complex (zinc only)
  - analysis/interface_metrics_zinc_best.csv / .md
      Best model per complex (zinc only)

Additional QC outputs:
  - analysis/zinc_coordination_models.csv / .md
  - analysis/zinc_finger_split_summary.csv / .md
  - analysis/zinc_entity_consistency_models.csv / .md
  - analysis/interface_metrics_zinc_avg_by_znfinger.csv / .md
"""
from __future__ import annotations

import argparse
import math
import tempfile
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align, rms
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.lib import distances
from MDAnalysis.topology import guessers
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.PDB import MMCIFParser, PDBIO

try:
    import freesasa  # type: ignore
except Exception:  # pragma: no cover
    freesasa = None

# Project paths
BASE_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = BASE_DIR / "Data"
DEFAULT_MODELS_DIR = BASE_DIR / "Models"

# Output paths (can be overridden by CLI flags; initialized in main())
OUTPUT_DIR = BASE_DIR / "results"
PLOT_DIR = OUTPUT_DIR / "plots"
CONVERTED_PDB_DIR = OUTPUT_DIR / "converted_pdb"

def init_output_dirs(output_dir: Path) -> None:
    global OUTPUT_DIR, PLOT_DIR, CONVERTED_PDB_DIR
    OUTPUT_DIR = output_dir
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR = OUTPUT_DIR / "plots"
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    CONVERTED_PDB_DIR = OUTPUT_DIR / "converted_pdb"
    CONVERTED_PDB_DIR.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="talk")

ALIGN_SEL = "protein and backbone and resid 1:70"
DOMAIN_SEL = "protein and resid 1:70"
DOMAIN_RMS_SEL = "protein and name CA and resid 1:70"
ZF1_RMS_SEL = "protein and name CA and resid 9:34"
ZF2_RMS_SEL = "protein and name CA and resid 47:66"
LINKER_CA_SEL = "protein and name CA and resid 35:46"
ZF1_CA_SEL = ZF1_RMS_SEL
ZF2_CA_SEL = ZF2_RMS_SEL
PROTEIN_HEAVY_SEL = "protein and resid 1:70 and not name H*"
RNA_HEAVY_SEL = "nucleic and not name H*"
RNA_BASE_SEL = (
    "nucleic and not name P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' O2' C1' H*"
)

DONOR_OXYGEN_NAMES = {"OG", "OG1", "OG2", "OH"}
AA_CLASS = {
    "positive": {"ARG", "LYS", "HIS"},
    "negative": {"ASP", "GLU"},
    "aromatic": {"PHE", "TYR", "TRP", "HIS"},
    "polar": {"SER", "THR", "ASN", "GLN", "TYR", "CYS"},
    "hydrophobic": {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR"},
}

RNA_INFO = {
    "gm_csf_wt": {
        "name": "GM-CSF 17nt (WT)",
        "sequence": "UGAAUUUAUUUAUUUAU",
        "status": "wt",
        "reference_label": "gm_csf_wt",
    },
    "gm_csf_gc_mutant": {
        "name": "GM-CSF 17nt (GC mutant)",
        "sequence": "GGAAUUUAUUUAUUUCU",
        "status": "mutant",
        "reference_label": "gm_csf_wt",
    },
    "nanomer_9_wt": {
        "name": "Nanomer 9nt (WT)",
        "sequence": "UAUAAAUAU",
        "status": "wt",
        "reference_label": "nanomer_9_wt",
    },
    "nanomer_9_mutant": {
        "name": "Nanomer 9nt (mutant / non-binder)",
        "sequence": "UACAAACAU",
        "status": "mutant",
        "reference_label": "nanomer_9_wt",
    },
    "human_lhr_20_mer_wt": {
        "name": "Human LHR 20nt (WT)",
        "sequence": "UAGGAAAUUAUUUAUUUUUA",
        "status": "wt",
        "reference_label": "human_lhr_20_mer_wt",
    },
    "human_lhr_20_mer_cc_mutant": {
        "name": "Human LHR 20nt (CC mutant)",
        "sequence": "UAGGAAAUUAUUUAUUCCUA",
        "status": "mutant",
        "reference_label": "human_lhr_20_mer_wt",
    },
    "mouse_apol11b_21_mer_wt": {
        "name": "Mouse Apol11b 21nt (WT)",
        "sequence": "UCGAAAUAUUUAUUUUCUUGA",
        "status": "wt",
        "reference_label": "mouse_apol11b_21_mer_wt",
    },
    "mouse_apol11b_21_mer_gc_mutant": {
        "name": "Mouse Apol11b 21nt (GC mutant)",
        "sequence": "UCGGAAUAUUUAUUUCCUUGA",
        "status": "mutant",
        "reference_label": "mouse_apol11b_21_mer_wt",
    },
    "nonamer_wt_2zn": {
        "name": "Nonamer 9nt (WT) + 2Zn",
        "sequence": "UAUAAAUAU",
        "status": "wt",
        "reference_label": "nonamer_wt_2zn",
    },
}

PROTEIN_LABEL_MAP = {
    "wt": "WT",
    "y145h": "Y145H",
    "y154h": "Y154H",
    "r160w": "R160W",
    "r160e": "R160E",
    "r184c": "R184C",
    "s185n": "S185N",
    "t187s": "T187S",
    "g204d": "G204D",
    "c206f": "C206F",
    "c208a": "C208A",
    "d210e": "D210E",
    "h212a": "H212A",
}

PROTEIN_EXP_CLASS_MAP = {
    "wt": "WT",
    "r160e": "deleterious",
    "r160w": "deleterious",
    "r184c": "deleterious",
    "g204d": "deleterious",
    "c206f": "deleterious",
    "c208a": "deleterious",
    "y145h": "neutral",
    "y154h": "deleterious",
    "s185n": "neutral",
    "t187s": "neutral",
    "d210e": "neutral",
    "h212a": "deleterious",
}
CLASS_ORDER = {
    "ProteinExpClass": ["WT", "neutral", "deleterious"],
    "RNAExpClass": ["WT", "mutant_bad"],
}
CLASS_COLOR_MAP = {
    "WT": "#4CAF50",
    "neutral": "#1E88E5",
    "deleterious": "#E53935",
    "mutant_bad": "#FB8C00",
}

IPTM_THRESHOLD = 0.5


def pdb_chain_order(pdb_path: Path) -> List[str]:
    """
    Extract chain IDs from a PDB file in first-appearance order.
    Works for our converted AF3 PDBs where each entity is a separate chain.
    """
    seen: List[str] = []
    try:
        for line in pdb_path.read_text().splitlines():
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            if len(line) < 22:
                continue
            ch = line[21].strip() or " "
            if ch not in seen:
                seen.append(ch)
    except Exception:
        return []
    return seen

# Zinc finger regions (domain numbering resid 1..70 corresponds to 151..220 in full-length)
ZF1_RANGE = (9, 34)
ZF2_RANGE = (47, 66)

# Multi-run support: AF3 server jobs can be repeated N times, each producing 5 models.
MODELS_PER_RUN = 5


def extract_run_index(folder_name: str) -> int:
    """Extract run number from a folder name like 'fold_... (3)'.  Returns 1 if absent."""
    m = re.search(r"\((\d+)\)\s*$", folder_name)
    return int(m.group(1)) if m else 1


def strip_run_suffix(folder_name: str) -> str:
    """Remove trailing ' (N)' from a folder name."""
    return re.sub(r"\s*\(\d+\)\s*$", "", folder_name)


def class_order(column: str, values: List[str]) -> List[str]:
    predefined = CLASS_ORDER.get(column)
    if predefined:
        return [cls for cls in predefined if cls in values]
    return sorted(values)


def class_palette(order: List[str]) -> Dict[str, str]:
    return {cls: CLASS_COLOR_MAP.get(cls, "#6a6a6a") for cls in order}


@dataclass
class ComplexMetadata:
    pdb_path: Path
    model_index: int
    protein_label: str
    rna_label: str
    protein_name: str
    rna_name: str
    rna_status: str
    protein_exp_class: str
    rna_exp_class: str
    reference_label: str
    complex_base: str


def parse_metadata(pdb_path: Path) -> ComplexMetadata:
    stem = pdb_path.stem  # fold_zfp36l2_<protein>_<rna>_model_<idx>[_runN]
    # Strip run suffix added during CIF->PDB conversion (e.g., _run3)
    stem = re.sub(r"_run\d+$", "", stem)
    if not stem.startswith("fold_zfp36l2_"):
        raise ValueError(f"Unexpected file name: {pdb_path.name}")
    core = stem[len("fold_zfp36l2_") :]
    m = re.match(r"(.+)_model_(\d+)$", core)
    if not m:
        raise ValueError(f"Unexpected model naming: {pdb_path.name}")
    core = m.group(1)
    model_index = int(m.group(2))
    complex_base = "fold_zfp36l2_" + core
    inferred_class = None
    core_tokens = core.lower().split("_")
    class_tokens = {"deleterious", "neutral", "wt"}
    is_mut = lambda tok: re.fullmatch(r"[a-z]\d+[a-z]", tok) is not None
    if len(core_tokens) == 2 and (
        (core_tokens[0] in class_tokens and is_mut(core_tokens[1]))
        or (is_mut(core_tokens[0]) and core_tokens[1] in class_tokens)
    ):
        # Some newer AF3 jobs are named like either:
        #   "deleterious_c208a" or "r160e_deleterious" (no RNA suffix).
        # Keep the variant identity and map RNA context to 2Zn nonamer.
        if core_tokens[0] in class_tokens:
            class_tok, mut_tok = core_tokens[0], core_tokens[1]
        else:
            class_tok, mut_tok = core_tokens[1], core_tokens[0]
        protein_label = mut_tok
        rna_label = "nonamer_wt_2zn"
        inferred_class = "WT" if class_tok == "wt" else class_tok
    else:
        protein_label, rna_part = core.split("_", 1)
        protein_label = protein_label.lower()
        rna_label = rna_part.lower()
    protein_name = PROTEIN_LABEL_MAP.get(protein_label, protein_label.upper())
    rna_info = RNA_INFO.get(rna_label, None)
    if rna_info is None:
        pretty = rna_label.replace("_", " ").title()
        rna_info = {
            "name": pretty,
            "sequence": "",
            "status": "wt",
            "reference_label": rna_label,
        }
        RNA_INFO[rna_label] = rna_info
    rna_name = rna_info["name"]
    rna_status = rna_info["status"]
    reference_label = rna_info.get("reference_label", rna_label)
    protein_exp_class = PROTEIN_EXP_CLASS_MAP.get(protein_label, inferred_class or "neutral")
    rna_exp_class = "WT" if rna_status == "wt" else "mutant_bad"
    return ComplexMetadata(
        pdb_path=pdb_path,
        model_index=model_index,
        protein_label=protein_label,
        rna_label=rna_label,
        protein_name=protein_name,
        rna_name=rna_name,
        rna_status=rna_status,
        protein_exp_class=protein_exp_class,
        rna_exp_class=rna_exp_class,
        reference_label=reference_label,
        complex_base=complex_base,
    )


def list_complexes(input_root: Path) -> List[ComplexMetadata]:
    structure_files = sorted(input_root.glob("**/*_model_*.pdb")) + sorted(
        input_root.glob("**/*_model_*.cif")
    )
    if not structure_files:
        raise RuntimeError(
            f"No structure files found under {input_root} (expected *_model_*.pdb or *_model_*.cif)"
        )
    results: List[ComplexMetadata] = []
    for orig_path in structure_files:
        run_idx = extract_run_index(orig_path.parent.name)
        pdb_path = ensure_structure_as_pdb(orig_path, run_idx)
        meta = parse_metadata(pdb_path)
        # Override model_index with a global index across all runs
        # (run 1 → 0-4, run 2 → 5-9, … , run 8 → 35-39)
        meta.model_index = (run_idx - 1) * MODELS_PER_RUN + meta.model_index
        results.append(meta)
    return results


def ensure_structure_as_pdb(structure_path: Path, run_index: int = 1) -> Path:
    """
    MDAnalysis can't read mmCIF here, so we convert CIF -> PDB once and cache it.
    A ``run_index`` suffix is appended to avoid name collisions when multiple
    AF3 runs share the same base name.
    """
    if structure_path.suffix.lower() == ".pdb":
        return structure_path
    if structure_path.suffix.lower() != ".cif":
        return structure_path
    out_path = CONVERTED_PDB_DIR / (structure_path.stem + f"_run{run_index}.pdb")
    if out_path.exists() and out_path.stat().st_size > 0:
        return out_path
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(structure_path.stem, str(structure_path))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path))
    return out_path


def pae_min_cross_from_matrix(chain_pair_pae_min: object) -> Optional[float]:
    if not isinstance(chain_pair_pae_min, list):
        return None
    try:
        n = len(chain_pair_pae_min)
        vals = []
        for i in range(n):
            row = chain_pair_pae_min[i]
            for j in range(n):
                if i == j:
                    continue
                vals.append(float(row[j]))
        if not vals:
            return None
        return float(sum(vals) / len(vals))
    except Exception:
        return None


def pae_offdiag_01_10(chain_pair_pae_min: object) -> Tuple[Optional[float], Optional[float]]:
    """
    For 2-chain complexes this corresponds to:
      [0][1] = chain0 -> chain1
      [1][0] = chain1 -> chain0
    """
    try:
        if not isinstance(chain_pair_pae_min, list) or len(chain_pair_pae_min) < 2:
            return None, None
        v01 = float(chain_pair_pae_min[0][1])
        v10 = float(chain_pair_pae_min[1][0])
        return v01, v10
    except Exception:
        return None, None


def load_job_chain_indices(input_root: Path) -> Dict[str, Dict[str, Optional[int]]]:
    """
    Map each AF3 run folder/base name to chain indices in confidence matrices by
    parsing the alphafoldserver job_request.json sequences order.

    Returns:
        base -> {'protein_idx': int, 'rna_idx': int, 'ion_idx': int|None, 'n_entities': int}
    """
    out: Dict[str, Dict[str, Optional[int]]] = {}
    for req_path in input_root.glob("**/*_job_request.json"):
        base = req_path.name.replace("_job_request.json", "")
        try:
            payload = json.loads(req_path.read_text())
            # request is a list with one dict
            req = payload[0] if isinstance(payload, list) and payload else payload
            seqs = req.get("sequences", [])
        except Exception:
            continue
        protein_idx = None
        rna_idx = None
        ion_idx = None
        for i, entry in enumerate(seqs):
            if isinstance(entry, dict) and "proteinChain" in entry:
                protein_idx = i
            elif isinstance(entry, dict) and "rnaSequence" in entry:
                rna_idx = i
            elif isinstance(entry, dict) and "ion" in entry:
                ion_idx = i
        out[base] = {
            "protein_idx": protein_idx,
            "rna_idx": rna_idx,
            "ion_idx": ion_idx,
            "n_entities": len(seqs) if isinstance(seqs, list) else None,
        }
    return out


def load_confidence_metrics(input_root: Path) -> Dict[str, Dict[int, Dict[str, Optional[float]]]]:
    """
    Parse AF3/AlphaFold Server confidence JSONs:
      - iptm (float)
      - ptm (float)
      - chain_pair_pae_min (matrix) -> we store cross-chain proxy as average off-diagonal

    Returns a mapping:
        complex_base -> model_index -> metrics
    """
    metrics: Dict[str, Dict[int, Dict[str, Optional[float]]]] = {}
    for json_path in input_root.glob("**/*_summary_confidences_*.json"):
        m = re.match(r"(.+)_summary_confidences_(\d+)\.json$", json_path.name)
        if not m:
            continue
        base = m.group(1)
        local_model_index = int(m.group(2))
        run_idx = extract_run_index(json_path.parent.name)
        model_index = (run_idx - 1) * MODELS_PER_RUN + local_model_index
        try:
            data = json.loads(json_path.read_text())
        except Exception:
            continue

        iptm = data.get("iptm", None)
        ptm = data.get("ptm", None)
        ranking_score = data.get("ranking_score", None)
        chain_ptm = data.get("chain_ptm", None)
        chain_pair_iptm = data.get("chain_pair_iptm", None)
        pae_matrix = data.get("chain_pair_pae_min", None)
        pae_01, pae_10 = pae_offdiag_01_10(pae_matrix)
        pae_mean = pae_min_cross_from_matrix(pae_matrix)

        # Robust inference (fixes the protein/RNA/Zn mis-indexing you pointed out):
        # - ion-like entity: chain_ptm[i] is null
        # - protein: entity with highest non-null chain_ptm
        # - RNA: entity j where PAE(protein->j) and PAE(j->protein) are both non-null
        n_entities = len(chain_ptm) if isinstance(chain_ptm, list) else None
        valid: List[int] = []
        ion_like: List[int] = []
        if isinstance(chain_ptm, list):
            for i, v in enumerate(chain_ptm):
                if v is None:
                    ion_like.append(i)
                else:
                    valid.append(i)

        pidx: int = 0
        if isinstance(chain_ptm, list) and valid:
            pidx = max(
                valid,
                key=lambda i: float(chain_ptm[i] if chain_ptm[i] is not None else -1.0),
            )

        candidates = [i for i in valid if i != pidx]
        ridx: Optional[int] = None

        def is_good_partner(j: int) -> bool:
            try:
                if not isinstance(pae_matrix, list):
                    return False
                a = pae_matrix[pidx][j]
                b = pae_matrix[j][pidx]
                return (a is not None) and (b is not None)
            except Exception:
                return False

        rna_candidates = [j for j in candidates if is_good_partner(j)]
        if len(rna_candidates) == 1:
            ridx = rna_candidates[0]
        elif candidates:
            # fallback: choose max ipTM(protein->j) among candidates with non-null PAE(protein->j)
            viable: List[int] = []
            if isinstance(pae_matrix, list):
                for j in candidates:
                    try:
                        if pae_matrix[pidx][j] is not None:
                            viable.append(j)
                    except Exception:
                        continue

            def iptm_score(j: int) -> float:
                try:
                    if not isinstance(chain_pair_iptm, list):
                        return -1.0
                    v = chain_pair_iptm[pidx][j]
                    return float(v) if v is not None else -1.0
                except Exception:
                    return -1.0

            ridx = max(viable, key=iptm_score) if viable else candidates[0]

        pair_iptm_pr = None
        pair_pae_pr = None
        pair_pae_rp = None

        if isinstance(chain_pair_iptm, list) and ridx is not None:
            try:
                v = chain_pair_iptm[pidx][ridx]
                pair_iptm_pr = float(v) if v is not None else None
            except Exception:
                pair_iptm_pr = None

        if isinstance(pae_matrix, list) and ridx is not None:
            try:
                v = pae_matrix[pidx][ridx]
                pair_pae_pr = float(v) if v is not None else None
            except Exception:
                pair_pae_pr = None
            try:
                v = pae_matrix[ridx][pidx]
                pair_pae_rp = float(v) if v is not None else None
            except Exception:
                pair_pae_rp = None

        metrics.setdefault(base, {})[model_index] = {
            "n_entities": float(n_entities) if n_entities is not None else None,
            "protein_idx": float(pidx) if pidx is not None else None,
            "rna_idx": float(ridx) if ridx is not None else None,
            "n_ion_like": float(len(ion_like)) if ion_like else 0.0,
            "ion_idx": float(ion_like[0]) if ion_like else None,
            "ion_idx_2": float(ion_like[1]) if len(ion_like) >= 2 else None,
            "iptm": float(iptm) if iptm is not None else None,
            "ptm": float(ptm) if ptm is not None else None,
            "ranking_score": float(ranking_score) if ranking_score is not None else None,
            # Correct protein↔RNA pair metrics (2×2 and 3×3 cases)
            "pair_iptm_protein_rna": pair_iptm_pr,
            "pair_pae_min_protein_to_rna": pair_pae_pr,
            "pair_pae_min_rna_to_protein": pair_pae_rp,
            # Legacy 2-chain convenience fields (kept for debugging)
            "chain_pair_pae_min_01": pae_01,
            "chain_pair_pae_min_10": pae_10,
            "chain_pair_pae_min": pae_01 if pae_01 is not None else pae_mean,
            "chain_pair_pae_min_cross_mean": pae_mean,
        }
    return metrics


def load_full_data_plddts(input_root: Path) -> Dict[str, Dict[int, List[float]]]:
    """
    Load per-atom pLDDT arrays from AlphaFold 3 ``full_data_*.json`` files.

    These JSON files are the **authoritative** source for per-atom pLDDT.
    The CIF ``B_iso_or_equiv`` column stores similar values but with minor
    rounding differences.

    Returns
    -------
    complex_base -> model_index -> list of per-atom pLDDT values
    """
    plddts: Dict[str, Dict[int, List[float]]] = {}
    for json_path in input_root.glob("**/*_full_data_*.json"):
        m = re.match(r"(.+)_full_data_(\d+)\.json$", json_path.name)
        if not m:
            continue
        base = m.group(1)
        local_model_index = int(m.group(2))
        run_idx = extract_run_index(json_path.parent.name)
        model_index = (run_idx - 1) * MODELS_PER_RUN + local_model_index
        try:
            data = json.loads(json_path.read_text())
        except Exception:
            continue
        atom_plddts = data.get("atom_plddts", None)
        if atom_plddts is not None and isinstance(atom_plddts, list):
            plddts.setdefault(base, {})[model_index] = atom_plddts
    return plddts


def build_reference_map(
    complexes: List[ComplexMetadata],
) -> Dict[str, Dict[int, Path]]:
    """
    Map each reference_label (e.g. nanomer_9_wt) to WT structure paths per model_index.
    """
    desired = {meta.reference_label for meta in complexes}
    refs: Dict[str, Dict[int, Path]] = {label: {} for label in desired}
    for meta in complexes:
        if meta.protein_label != "wt":
            continue
        # only use WT + reference RNA as the reference set
        if meta.rna_label != meta.reference_label:
            continue
        if meta.rna_label in refs:
            refs[meta.rna_label][meta.model_index] = meta.pdb_path

    missing = [label for label, by_model in refs.items() if not by_model]
    if missing:
        raise RuntimeError(
            "Missing WT reference structures for reference labels: "
            + ", ".join(sorted(missing))
        )
    return refs


def load_canonical_sequence(ref_path: Path) -> Dict[int, str]:
    u_ref = mda.Universe(ref_path)
    seq_map = {}
    protein = u_ref.select_atoms(DOMAIN_SEL)
    for residue in protein.residues:
        seq_map[residue.resid] = residue.resname
    return seq_map


def ensure_bonds(u: mda.Universe) -> None:
    try:
        _ = u.bonds
        return
    except NoDataError:
        pass
    guessed = guessers.guess_bonds(u.atoms, u.atoms.positions)
    if guessed:
        u.add_TopologyAttr("bonds", guessed)


def atom_is_donor(atom: mda.core.groups.Atom) -> bool:
    element = (atom.element or "").upper()
    name = atom.name.strip()
    if element == "N":
        return True
    if element == "O" and name in DONOR_OXYGEN_NAMES:
        return True
    return False


def atom_is_acceptor(atom: mda.core.groups.Atom) -> bool:
    element = (atom.element or "").upper()
    if element == "O":
        return True
    if element == "N" and atom.resname.strip() in {"A", "U", "G", "C"}:
        return True
    return False


def donor_pairs(atoms: mda.core.groups.AtomGroup) -> List[Tuple[mda.core.groups.Atom, np.ndarray]]:
    pairs = []
    for atom in atoms:
        if not atom_is_donor(atom):
            continue
        bonded = [b for b in atom.bonded_atoms if (b.element or "").upper() != "H"]
        if not bonded:
            continue
        base = bonded[0]
        base_vec = atom.position - base.position
        norm = np.linalg.norm(base_vec)
        if norm < 1e-3:
            continue
        pairs.append((atom, base_vec / norm))
    return pairs


def acceptor_atoms(atoms: mda.core.groups.AtomGroup) -> List[mda.core.groups.Atom]:
    return [atom for atom in atoms if atom_is_acceptor(atom)]


def mean_plddt_from_atoms(atoms: mda.core.groups.AtomGroup) -> float:
    """
    AlphaFold commonly stores pLDDT in the PDB/mmCIF temperature factor field.
    We report mean pLDDT over a selection (typically protein CA atoms).

    Returns NaN if tempfactors are unavailable.
    """
    try:
        if atoms is None or len(atoms) == 0:
            return float("nan")
        # MDAnalysis exposes PDB B-factors as tempfactors
        vals = np.array(atoms.tempfactors, dtype=float)
        vals = vals[np.isfinite(vals)]
        return float(vals.mean()) if len(vals) else float("nan")
    except Exception:
        return float("nan")


def plddt_subregion_stats(
    protein_ca: mda.core.groups.AtomGroup,
) -> Dict[str, float]:
    """
    Compute per-subregion and per-residue pLDDT statistics.

    All metrics use **Cα atoms** (one per residue), which is the standard
    per-residue pLDDT convention.  This matches the AlphaFold 3 server's
    colour scheme: the server colours the structure by per-residue (Cα)
    pLDDT, so ``pLDDT_min`` here corresponds to what you actually see on
    the server.

    Returns dict with keys:
        pLDDT_ZF1, pLDDT_ZF2, pLDDT_linker,
        pLDDT_min (lowest per-residue Cα pLDDT in the protein domain),
        pLDDT_min_resid (which resid has the minimum).
    """
    nan = float("nan")
    out: Dict[str, float] = {
        "pLDDT_ZF1": nan,
        "pLDDT_ZF2": nan,
        "pLDDT_linker": nan,
        "pLDDT_min": nan,
        "pLDDT_min_resid": nan,
    }
    try:
        if protein_ca is None or len(protein_ca) == 0:
            return out
        # Build per-residue pLDDT (one CA per residue)
        resid_plddt: Dict[int, float] = {}
        for atom in protein_ca:
            b = float(atom.tempfactor)
            if math.isfinite(b):
                resid_plddt[atom.resid] = b

        if not resid_plddt:
            return out

        def _mean_range(start: int, end: int) -> float:
            vals = [resid_plddt[r] for r in range(start, end + 1) if r in resid_plddt]
            return float(np.mean(vals)) if vals else nan

        out["pLDDT_ZF1"] = _mean_range(*ZF1_RANGE)
        out["pLDDT_ZF2"] = _mean_range(*ZF2_RANGE)
        out["pLDDT_linker"] = _mean_range(35, 46)

        # Minimum per-residue (Cα) pLDDT across the entire domain
        min_resid = min(resid_plddt, key=resid_plddt.get)  # type: ignore[arg-type]
        out["pLDDT_min"] = resid_plddt[min_resid]
        out["pLDDT_min_resid"] = float(min_resid)
    except Exception:
        pass
    return out


def mean_plddt_from_pdb_text(pdb_path: Path) -> float:
    """
    Fallback parser: compute mean B-factor over protein CA atoms directly from PDB text.
    This is intentionally conservative and returns NaN if nothing is found.
    """
    try:
        vals: List[float] = []
        for line in pdb_path.read_text().splitlines():
            if not line.startswith("ATOM"):
                continue
            # Atom name (cols 13-16), residue name (cols 18-20), B-factor (cols 61-66)
            if len(line) < 66:
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            resname = line[17:20].strip().upper()
            # skip nucleic acids / unknowns
            if resname in {"A", "U", "G", "C"}:
                continue
            try:
                b = float(line[60:66].strip())
            except Exception:
                continue
            if math.isfinite(b):
                vals.append(b)
        return float(sum(vals) / len(vals)) if vals else float("nan")
    except Exception:
        return float("nan")


def compute_hbond_count(
    protein_atoms: mda.core.groups.AtomGroup,
    rna_atoms: mda.core.groups.AtomGroup,
    distance_cutoff: float = 3.5,
    angle_cutoff: float = 120.0,
) -> int:
    protein_donors = donor_pairs(protein_atoms)
    rna_donors = donor_pairs(rna_atoms)
    protein_acceptors = acceptor_atoms(protein_atoms)
    rna_acceptors = acceptor_atoms(rna_atoms)

    def count_pairs(
        donor_list: List[Tuple[mda.core.groups.Atom, np.ndarray]],
        acceptor_list: List[mda.core.groups.Atom],
    ) -> int:
        if not donor_list or not acceptor_list:
            return 0
        donor_positions = np.array([atom.position for atom, _ in donor_list])
        acceptor_positions = np.array([atom.position for atom in acceptor_list])
        dist_matrix = distances.distance_array(donor_positions, acceptor_positions)
        hits = np.argwhere(dist_matrix <= distance_cutoff)
        pairs = set()
        for donor_idx, acc_idx in hits:
            donor_atom, donor_vec = donor_list[donor_idx]
            acceptor_atom = acceptor_list[acc_idx]
            vec = acceptor_atom.position - donor_atom.position
            norm = np.linalg.norm(vec)
            if norm < 1e-3:
                continue
            vec /= norm
            cos_angle = np.clip(np.dot(donor_vec, vec), -1.0, 1.0)
            angle = math.degrees(math.acos(cos_angle))
            if angle >= angle_cutoff:
                pairs.add(
                    (
                        donor_atom.resid,
                        acceptor_atom.resid,
                        donor_atom.segid,
                        acceptor_atom.segid,
                    )
                )
        return len(pairs)

    total = 0
    total += count_pairs(protein_donors, rna_acceptors)
    total += count_pairs(rna_donors, protein_acceptors)
    return total


def count_contacts(
    protein_atoms: mda.core.groups.AtomGroup,
    rna_atoms: mda.core.groups.AtomGroup,
    cutoff: float = 4.0,
) -> int:
    if not protein_atoms or not rna_atoms:
        return 0
    dist_matrix = distances.distance_array(
        protein_atoms.positions, rna_atoms.positions
    )
    contact_indices = np.argwhere(dist_matrix <= cutoff)
    residue_pairs = {
        (protein_atoms[i].resid, rna_atoms[j].resid)
        for i, j in contact_indices
    }
    return len(residue_pairs)


def write_temp_pdb(atom_group: mda.core.groups.AtomGroup) -> Optional[str]:
    if not atom_group:
        return None
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    tmp.close()
    atom_group.write(tmp.name)
    return tmp.name


def compute_bsa(
    pdb_path: Path,
    protein_atoms: mda.core.groups.AtomGroup,
    rna_atoms: mda.core.groups.AtomGroup,
) -> float:
    # Optional dependency: if FreeSASA isn't installed, skip BSA but keep the rest of the metrics.
    if freesasa is None:
        return float("nan")
    # Zinc runs include ZN ions; FreeSASA might not know their radii.
    # We skip unknown atoms to avoid warnings and keep BSA comparable.
    opts = {"hetatm": True, "hydrogen": False, "skip-unknown": True}
    complex_structure = freesasa.Structure(str(pdb_path), options=opts)
    complex_result = freesasa.calc(complex_structure)
    complex_area = complex_result.totalArea()

    protein_file = write_temp_pdb(protein_atoms)
    rna_file = write_temp_pdb(rna_atoms)
    try:
        protein_area = 0.0
        rna_area = 0.0
        if protein_file:
            protein_result = freesasa.calc(freesasa.Structure(protein_file, options=opts))
            protein_area = protein_result.totalArea()
        if rna_file:
            rna_result = freesasa.calc(freesasa.Structure(rna_file, options=opts))
            rna_area = rna_result.totalArea()
    finally:
        for tmp in (protein_file, rna_file):
            if tmp and Path(tmp).exists():
                Path(tmp).unlink()

    return protein_area + rna_area - complex_area


def rmsd_between(
    ref_u: mda.Universe, mobile_u: mda.Universe, selection: str
) -> float:
    ref_atoms = ref_u.select_atoms(selection)
    mob_atoms = mobile_u.select_atoms(selection)
    if len(ref_atoms) != len(mob_atoms):
        raise ValueError(f"Selection mismatch for {selection}")
    return float(rms.rmsd(mob_atoms.positions, ref_atoms.positions))


def detect_mutations(
    mobile_u: mda.Universe, canonical_seq: Dict[int, str]
) -> List[Tuple[int, str, str]]:
    mutations = []
    residues = mobile_u.select_atoms(DOMAIN_SEL).residues
    for res in residues:
        resid = res.resid
        wt_res = canonical_seq.get(resid)
        if wt_res and res.resname != wt_res:
            mutations.append((resid, wt_res, res.resname))
    return mutations


def classify_change(wt_res: str, mut_res: str) -> str:
    notes = []
    if wt_res in AA_CLASS["positive"] and mut_res not in AA_CLASS["positive"]:
        notes.append("lost positive charge")
    if wt_res in AA_CLASS["negative"] and mut_res not in AA_CLASS["negative"]:
        notes.append("lost negative charge")
    if wt_res in AA_CLASS["aromatic"] and mut_res not in AA_CLASS["aromatic"]:
        notes.append("lost aromatic stacking")
    if mut_res in AA_CLASS["positive"] and wt_res not in AA_CLASS["positive"]:
        notes.append("gained positive charge")
    if mut_res in AA_CLASS["negative"] and wt_res not in AA_CLASS["negative"]:
        notes.append("gained negative charge")
    if not notes:
        notes.append("side-chain chemistry conserved")
    return "; ".join(notes)


def mutation_local_effects(
    mobile_u: mda.Universe,
    mutations: List[Tuple[int, str, str]],
    rna_base_atoms: mda.core.groups.AtomGroup,
) -> List[str]:
    notes = []
    for resid, wt_res, mut_res in mutations:
        residue = mobile_u.select_atoms(f"protein and resid {resid}").residues[0]
        sidechain = residue.atoms.select_atoms("not name N C O CA OXT")
        if not sidechain:
            continue
        centroid = sidechain.positions.mean(axis=0)
        ca = residue.atoms.select_atoms("name CA")
        if not ca:
            continue
        ca_pos = ca.positions[0]
        vec_sidechain = centroid - ca_pos
        vec_sidechain_norm = np.linalg.norm(vec_sidechain)
        if vec_sidechain_norm < 1e-3:
            continue
        vec_sidechain /= vec_sidechain_norm
        distances_to_bases = distances.distance_array(
            sidechain.positions, rna_base_atoms.positions
        )
        min_idx = np.unravel_index(np.argmin(distances_to_bases), distances_to_bases.shape)
        min_dist = float(distances_to_bases[min_idx])
        closest_atom = rna_base_atoms[min_idx[1]]
        vec_to_base = closest_atom.position - ca_pos
        vec_to_base_norm = np.linalg.norm(vec_to_base)
        orientation = "undetermined"
        if vec_to_base_norm >= 1e-3:
            vec_to_base /= vec_to_base_norm
            angle = math.degrees(
                math.acos(np.clip(np.dot(vec_sidechain, vec_to_base), -1.0, 1.0))
            )
            orientation = "toward RNA" if angle <= 90 else "solvated"
        change_note = classify_change(wt_res, mut_res)
        abs_resid = resid + 150
        notes.append(
            f"{wt_res}{abs_resid}{mut_res}: {min_dist:.2f} Å to {closest_atom.resname}{closest_atom.resid} base "
            f"({orientation}); {change_note}"
        )
    return notes


def residues_to_rna_sequence(rna_residues: mda.core.groups.ResidueGroup) -> str:
    # mdanalysis uses resnames like 'A','U','G','C' for RNA bases
    return "".join([res.resname.strip() for res in rna_residues])


def model_entity_signature(u: mda.Universe) -> Dict[str, object]:
    protein = u.select_atoms("protein")
    rna = u.select_atoms("nucleic")
    zn = u.select_atoms("resname ZN or name ZN")
    rna_seq = residues_to_rna_sequence(rna.residues)
    return {
        "n_protein_res": int(len(protein.residues)),
        "n_rna_res": int(len(rna.residues)),
        "n_zn_atoms": int(len(zn)),
        "rna_seq": rna_seq,
        "rna_seq_rev": rna_seq[::-1],
        "protein_chain_ids": ",".join(sorted(set(protein.chainIDs))) if hasattr(protein, "chainIDs") else "",
        "rna_chain_ids": ",".join(sorted(set(rna.chainIDs))) if hasattr(rna, "chainIDs") else "",
        "zn_chain_ids": ",".join(sorted(set(zn.chainIDs))) if hasattr(zn, "chainIDs") else "",
    }


def zinc_coordination_table(
    u: mda.Universe, zn_atom: mda.core.groups.Atom, finger_range: Tuple[int, int]
) -> Optional[Dict[str, object]]:
    """
    Choose coordinating residues inside a finger range.
    - Cys: first 3 Cys residues (SG)
    - His: choose the His residue whose min(Zn-NE2, Zn-ND1) is smallest (closest His)

    Returns dict with residue IDs and distances.
    """
    start, end = finger_range
    cys_atoms = u.select_atoms(
        f"protein and resid {start}:{end} and resname CYS and name SG"
    )
    cys_resids = sorted(set(cys_atoms.resids))
    if len(cys_resids) < 3:
        return None
    c1, c2, c3 = cys_resids[:3]

    his_resids = sorted(
        set(u.select_atoms(f"protein and resid {start}:{end} and resname HIS").resids)
    )
    if not his_resids:
        return None

    def atom_or_none(sel: str):
        ag = u.select_atoms(sel)
        return ag[0] if len(ag) else None

    best_h = None
    best_h_score = None
    best_ne2 = None
    best_nd1 = None
    for h in his_resids:
        ne2 = atom_or_none(f"protein and resid {h} and name NE2")
        nd1 = atom_or_none(f"protein and resid {h} and name ND1")
        d_ne2 = float(np.linalg.norm(zn_atom.position - ne2.position)) if ne2 else None
        d_nd1 = float(np.linalg.norm(zn_atom.position - nd1.position)) if nd1 else None
        d_min = min([d for d in [d_ne2, d_nd1] if d is not None], default=None)
        if d_min is None:
            continue
        if best_h_score is None or d_min < best_h_score:
            best_h_score = d_min
            best_h = h
            best_ne2 = d_ne2
            best_nd1 = d_nd1

    if best_h is None:
        return None

    # Cys distances
    d_c1 = float(
        np.linalg.norm(zn_atom.position - u.select_atoms(f"protein and resid {c1} and name SG")[0].position)
    )
    d_c2 = float(
        np.linalg.norm(zn_atom.position - u.select_atoms(f"protein and resid {c2} and name SG")[0].position)
    )
    d_c3 = float(
        np.linalg.norm(zn_atom.position - u.select_atoms(f"protein and resid {c3} and name SG")[0].position)
    )

    d_h_min = min([d for d in [best_ne2, best_nd1] if d is not None], default=None)
    return {
        "Cys1_resid": c1,
        "Cys2_resid": c2,
        "Cys3_resid": c3,
        "His_resid": best_h,
        "Zn_Cys1_SG": d_c1,
        "Zn_Cys2_SG": d_c2,
        "Zn_Cys3_SG": d_c3,
        "Zn_His_NE2": best_ne2,
        "Zn_His_ND1": best_nd1,
        "Zn_His_min": d_h_min,
    }


def compute_zinc_coordination_rows(
    meta: ComplexMetadata,
    pdb_path: Path,
) -> List[Dict[str, object]]:
    """
    For zinc runs, compute Zn->(Cys SG / His NE2) distances to sanity-check coordination.
    If only one Zn is present, we report which finger it best matches (ZF1 vs ZF2).
    """
    u = mda.Universe(pdb_path)
    zn_atoms = u.select_atoms("resname ZN or name ZN")
    if len(zn_atoms) == 0:
        return []

    rows: List[Dict[str, object]] = []
    for i, zn in enumerate(zn_atoms):
        best = None
        for label, frange in (("ZF1", ZF1_RANGE), ("ZF2", ZF2_RANGE)):
            cand = zinc_coordination_table(u, zn, frange)
            if cand is None:
                continue
            # conservative score: max of 3 Cys distances and min(His) distance
            score = max(
                cand["Zn_Cys1_SG"],
                cand["Zn_Cys2_SG"],
                cand["Zn_Cys3_SG"],
                cand["Zn_His_min"] if cand["Zn_His_min"] is not None else 999.0,
            )
            if best is None or score < best[0]:
                best = (score, label, cand)

        comment = ""
        if best is None:
            comment = "missing expected Cys/His atoms in both ZF regions"
            rows.append(
                {
                    "ComplexBase": meta.complex_base,
                    "Model": meta.model_index,
                    "ZN_index": i,
                    "ZN_chain": zn.chainID if hasattr(zn, "chainID") else "",
                    "ZN_resid": int(zn.resid),
                    "AssignedFinger": None,
                    "Cys1_resid": None,
                    "Cys2_resid": None,
                    "Cys3_resid": None,
                    "His_resid": None,
                    "Zn_Cys1_SG": None,
                    "Zn_Cys2_SG": None,
                    "Zn_Cys3_SG": None,
                    "Zn_His_NE2": None,
                    "Zn_His_ND1": None,
                    "coord_atoms_within_2p6": None,
                    "comment": comment,
                }
            )
            continue

        _, assigned, cand = best
        ligand_atom = None
        if cand["Zn_His_NE2"] is not None and cand["Zn_His_NE2"] <= 2.6:
            if cand["Zn_His_ND1"] is None or cand["Zn_His_NE2"] <= cand["Zn_His_ND1"]:
                ligand_atom = "NE2"
        if ligand_atom is None and cand["Zn_His_ND1"] is not None and cand["Zn_His_ND1"] <= 2.6:
            ligand_atom = "ND1"
        cys = [cand["Zn_Cys1_SG"], cand["Zn_Cys2_SG"], cand["Zn_Cys3_SG"]]
        his_min = cand["Zn_His_min"]
        # coordination count among the 4 expected ligands (3 Cys + 1 His N)
        coord_atoms_within = 0
        coord_atoms_within += sum(1 for d in cys if d is not None and d <= 2.6)
        coord_atoms_within += 1 if (his_min is not None and his_min <= 2.6) else 0

        if his_min is None or any(d is None for d in cys):
            comment = "coord atoms missing"
        elif coord_atoms_within == 4:
            comment = "coordination plausible (4 atoms within 2.6 Å)"
        elif coord_atoms_within == 3:
            comment = "3-coordinate only (His not within 2.6 Å)"
        else:
            comment = "NOT coordinating (insufficient atoms within 2.6 Å)"
        rows.append(
            {
                "ComplexBase": meta.complex_base,
                "Model": meta.model_index,
                "ZN_index": i,
                "ZN_chain": zn.chainID if hasattr(zn, "chainID") else "",
                "ZN_resid": int(zn.resid),
                "AssignedFinger": assigned,
                "ZnLigandAtom": ligand_atom,
                "Cys1_resid": cand["Cys1_resid"],
                "Cys2_resid": cand["Cys2_resid"],
                "Cys3_resid": cand["Cys3_resid"],
                "His_resid": cand["His_resid"],
                "Zn_Cys1_SG": cys[0],
                "Zn_Cys2_SG": cys[1],
                "Zn_Cys3_SG": cys[2],
                "Zn_His_NE2": cand["Zn_His_NE2"],
                "Zn_His_ND1": cand["Zn_His_ND1"],
                "coord_atoms_within_2p6": coord_atoms_within,
                "comment": comment,
            }
        )
    return rows


def plot_rmsd_box(
    df: pd.DataFrame, group_col: str, filename: str, title: str
) -> Path:
    order = class_order(group_col, df[group_col].unique().tolist())
    palette = class_palette(order)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    sns.boxplot(
        data=df,
        x=group_col,
        y="RMSD_whole",
        hue=group_col,
        order=order,
        hue_order=order,
        palette=palette,
        dodge=False,
        ax=ax,
        legend=False,
    )
    sns.stripplot(
        data=df,
        x=group_col,
        y="RMSD_whole",
        order=order,
        dodge=True,
        color="#333333",
        size=6,
        ax=ax,
    )
    ax.set_ylabel("Backbone RMSD vs WT (Å)")
    ax.set_xlabel("")
    ax.set_title(title)
    ax.set_ylim(bottom=-0.05)
    fig.tight_layout()
    out_path = PLOT_DIR / filename
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


def plot_rmsd_contacts(df: pd.DataFrame) -> Path:
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    palette = class_palette(
        class_order("ProteinExpClass", df["ProteinExpClass"].unique().tolist())
    )
    sns.scatterplot(
        data=df,
        x="RMSD_whole",
        y="Contacts",
        hue="ProteinExpClass",
        hue_order=class_order(
            "ProteinExpClass", df["ProteinExpClass"].unique().tolist()
        ),
        style="RNAExpClass",
        style_order=class_order(
            "RNAExpClass", df["RNAExpClass"].unique().tolist()
        ),
        size="BSA",
        sizes=(50, 300),
        palette=palette,
        ax=ax,
    )
    ax.set_xlabel("Backbone RMSD vs WT (Å)")
    ax.set_ylabel("Protein–RNA contacts (<4 Å)")
    ax.set_title("RMSD vs interface packing")
    ax.set_xlim(left=-0.05)
    fig.tight_layout()
    out_path = PLOT_DIR / "rmsd_vs_contacts.png"
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


def plot_interface_summary(
    df: pd.DataFrame, group_col: str, filename: str, title: str
) -> Path:
    agg = (
        df.groupby(group_col)[["HBonds", "Contacts", "BSA"]]
        .mean()
        .reset_index()
    )
    tidy = agg.melt(id_vars=group_col, var_name="Metric", value_name="Value")
    order = class_order(group_col, agg[group_col].tolist())
    palette = class_palette(order)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sns.barplot(
        data=tidy,
        x="Metric",
        y="Value",
        hue=group_col,
        hue_order=order,
        palette=palette,
        ax=ax,
    )
    ax.set_ylabel("Mean value per class")
    ax.set_title(title)
    fig.tight_layout()
    out_path = PLOT_DIR / filename
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


def make_plots(df: pd.DataFrame) -> List[Path]:
    return [
        plot_rmsd_box(
            df,
            group_col="ProteinExpClass",
            filename="rmsd_by_protein_class.png",
            title="Backbone deviations by protein class",
        ),
        plot_rmsd_contacts(df),
        plot_interface_summary(
            df,
            group_col="ProteinExpClass",
            filename="interface_metric_by_protein_class.png",
            title="Interface metrics by protein class",
        ),
        plot_interface_summary(
            df,
            group_col="RNAExpClass",
            filename="interface_metric_by_rna_class.png",
            title="Interface metrics by RNA class",
        ),
    ]


def analyze_complex(
    meta: ComplexMetadata,
    ref_universe: mda.Universe,
    canonical_seq: Dict[int, str],
    confidence_metrics: Dict[str, Dict[int, Dict[str, Optional[float]]]],
    full_data_plddts: Optional[Dict[str, Dict[int, List[float]]]] = None,
) -> Dict[str, object]:
    mobile_u = mda.Universe(meta.pdb_path)
    align.alignto(mobile_u, ref_universe, select=ALIGN_SEL)
    ensure_bonds(mobile_u)
    protein_atoms = mobile_u.select_atoms(PROTEIN_HEAVY_SEL)
    rna_atoms = mobile_u.select_atoms(RNA_HEAVY_SEL)
    rna_base_atoms = mobile_u.select_atoms(RNA_BASE_SEL)
    protein_ca = mobile_u.select_atoms("protein and name CA")

    # Inject authoritative per-atom pLDDT from full_data JSON into the
    # universe tempfactors.  The CIF B_iso_or_equiv column carries similar
    # values but with minor rounding differences; the JSON is the source of
    # truth for AlphaFold 3.
    if full_data_plddts is not None:
        json_plddts = full_data_plddts.get(meta.complex_base, {}).get(
            meta.model_index, None
        )
        if json_plddts is not None and len(json_plddts) == len(mobile_u.atoms):
            mobile_u.atoms.tempfactors = np.array(json_plddts, dtype=float)

    plddt = mean_plddt_from_atoms(protein_ca)
    # Fallback in case MDAnalysis couldn't populate tempfactors for a given file.
    if not math.isfinite(plddt):
        plddt = mean_plddt_from_pdb_text(meta.pdb_path)
    # Per-subregion and per-residue pLDDT metrics (all Cα-based)
    plddt_stats = plddt_subregion_stats(protein_ca)

    whole_rmsd = rmsd_between(ref_universe, mobile_u, DOMAIN_RMS_SEL)
    zf1_rmsd = rmsd_between(ref_universe, mobile_u, ZF1_RMS_SEL)
    zf2_rmsd = rmsd_between(ref_universe, mobile_u, ZF2_RMS_SEL)

    hbonds = compute_hbond_count(protein_atoms, rna_atoms)
    contacts = count_contacts(protein_atoms, rna_atoms)
    bsa = compute_bsa(meta.pdb_path, protein_atoms, rna_atoms)

    mutations = detect_mutations(mobile_u, canonical_seq)
    mutation_notes = mutation_local_effects(mobile_u, mutations, rna_base_atoms)

    notes = []
    if mutation_notes:
        notes.extend(mutation_notes)
    if meta.rna_status == "mutant":
        notes.append("RNA substitution expected to weaken binding")

    conf = confidence_metrics.get(meta.complex_base, {}).get(meta.model_index, {})
    # Convenience: cross-mean PAE for protein↔RNA when both directions exist
    pae_pr = conf.get("pair_pae_min_protein_to_rna", None)
    pae_rp = conf.get("pair_pae_min_rna_to_protein", None)
    pae_cross_mean = None
    try:
        vals = [v for v in [pae_pr, pae_rp] if v is not None]
        pae_cross_mean = float(sum(vals) / len(vals)) if vals else None
    except Exception:
        pae_cross_mean = None
    return {
        "ComplexBase": meta.complex_base,
        "Model": meta.model_index,
        "pdb": meta.pdb_path.name,
        "Protein": meta.protein_name,
        "RNA": meta.rna_name,
        "ProteinExpClass": meta.protein_exp_class,
        "RNAExpClass": meta.rna_exp_class,
        "n_entities": conf.get("n_entities", None),
        "protein_idx": conf.get("protein_idx", None),
        "rna_idx": conf.get("rna_idx", None),
        "n_ion_like": conf.get("n_ion_like", None),
        "ion_idx": conf.get("ion_idx", None),
        "ion_idx_2": conf.get("ion_idx_2", None),
        "ranking_score": conf.get("ranking_score", None),
        "ipTM": conf.get("iptm", None),
        "pTM": conf.get("ptm", None),
        "pair_ipTM_protein_rna": conf.get("pair_iptm_protein_rna", None),
        "pair_pae_min_protein_to_rna": conf.get(
            "pair_pae_min_protein_to_rna", None
        ),
        "pair_pae_min_rna_to_protein": conf.get(
            "pair_pae_min_rna_to_protein", None
        ),
        "pair_pae_min_cross_mean": pae_cross_mean,
        "chain_pair_pae_min": conf.get("chain_pair_pae_min", None),
        "chain_pair_pae_min_01": conf.get("chain_pair_pae_min_01", None),
        "chain_pair_pae_min_10": conf.get("chain_pair_pae_min_10", None),
        "chain_pair_pae_min_cross_mean": conf.get(
            "chain_pair_pae_min_cross_mean", None
        ),
        "pLDDT": plddt,
        "pLDDT_ZF1": plddt_stats["pLDDT_ZF1"],
        "pLDDT_ZF2": plddt_stats["pLDDT_ZF2"],
        "pLDDT_linker": plddt_stats["pLDDT_linker"],
        "pLDDT_min": plddt_stats["pLDDT_min"],
        "pLDDT_min_resid": plddt_stats["pLDDT_min_resid"],
        "RMSD_whole": whole_rmsd,
        "RMSD_ZF1": zf1_rmsd,
        "RMSD_ZF2": zf2_rmsd,
        "HBonds": hbonds,
        "Contacts": contacts,
        "BSA": bsa,
        "Notes": " | ".join(notes) if notes else "",
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute structural + AF3 confidence interface metrics for ZFP36L2 variants."
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=None,
        help=(
            "Root directory containing AF3 runs. Defaults to ./Models if present, else ./Data."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(BASE_DIR / "results"),
        help="Output directory for CSV/MD tables and plots (default: ./results).",
    )
    args = parser.parse_args()

    input_root = Path(args.input_dir) if args.input_dir else (
        DEFAULT_MODELS_DIR if DEFAULT_MODELS_DIR.exists() else DEFAULT_DATA_DIR
    )
    output_root = Path(args.output_dir)
    init_output_dirs(output_root)

    complexes = list_complexes(input_root)
    refs = build_reference_map(complexes)
    # canonical residue names (for mutation detection) from any reference model
    any_ref_path = next(iter(next(iter(refs.values())).values()))
    canonical_seq = load_canonical_sequence(any_ref_path)
    ref_cache: Dict[str, mda.Universe] = {}
    confidence_metrics = load_confidence_metrics(input_root)
    full_data_plddts = load_full_data_plddts(input_root)

    records = []
    zinc_coord_rows: List[Dict[str, object]] = []
    entity_rows: List[Dict[str, object]] = []
    for meta in complexes:
        ref_label = meta.reference_label
        ref_by_model = refs[ref_label]
        # Prefer matching model index reference; fall back to model_0 (or any available)
        ref_path = ref_by_model.get(meta.model_index) or ref_by_model.get(0) or next(
            iter(ref_by_model.values())
        )
        ref_cache_key = f"{ref_label}:{Path(ref_path).name}"
        if ref_cache_key not in ref_cache:
            ref_cache[ref_cache_key] = mda.Universe(ref_path)
        record = analyze_complex(
            meta,
            ref_cache[ref_cache_key],
            canonical_seq,
            confidence_metrics,
            full_data_plddts,
        )
        records.append(record)

        # Entity consistency + Zn coordination checks (metal runs only: _zinc, _2zn)
        if meta.complex_base.endswith("_zinc") or meta.complex_base.endswith("_2zn"):
            u_tmp = mda.Universe(meta.pdb_path)
            sig = model_entity_signature(u_tmp)
            entity_rows.append(
                {
                    "ComplexBase": meta.complex_base,
                    "Model": meta.model_index,
                    **sig,
                }
            )
            zinc_coord_rows.extend(compute_zinc_coordination_rows(meta, meta.pdb_path))

    df_models = pd.DataFrame(records)
    df_models.sort_values(
        ["ComplexBase", "Model"], inplace=True
    )

    # Per-model outputs (all complexes)
    csv_models = OUTPUT_DIR / "interface_metrics_models.csv"
    md_models = OUTPUT_DIR / "interface_metrics_models.md"
    df_models.to_csv(csv_models, index=False)
    md_models.write_text(
        "\n".join(
            [
                "# ZFP36L2 AlphaFold3 interface metrics (per-model)",
                "",
                f"- Rows: {len(df_models)} (one per AF3 model)",
                "- RMSD columns are in Å; BSA in Å²",
                "",
                df_models.to_markdown(index=False),
            ]
        )
    )

    # Zinc sanity check reports
    if entity_rows:
        df_entity = pd.DataFrame(entity_rows).sort_values(["ComplexBase", "Model"])
        ent_csv = OUTPUT_DIR / "zinc_entity_consistency_models.csv"
        ent_md = OUTPUT_DIR / "zinc_entity_consistency_models.md"
        df_entity.to_csv(ent_csv, index=False)
        ent_md.write_text(
            "\n".join(
                [
                    "# Zinc runs: entity consistency per model",
                    "",
                    "- `rna_seq` is extracted from the model residue order; compare across models for 5'→3' consistency.",
                    "- `n_zn_atoms` should be stable across models if Zn is defined consistently.",
                    "",
                    df_entity.to_markdown(index=False),
                ]
            )
        )

    if zinc_coord_rows:
        df_zn = pd.DataFrame(zinc_coord_rows).sort_values(
            ["ComplexBase", "Model", "ZN_index"]
        )
        zn_csv = OUTPUT_DIR / "zinc_coordination_models.csv"
        zn_md = OUTPUT_DIR / "zinc_coordination_models.md"
        df_zn.to_csv(zn_csv, index=False)
        zn_md.write_text(
            "\n".join(
                [
                    "# Zinc runs: Zn coordination sanity check (per model)",
                    "",
                    "- Distances are in Å.",
                    "- Expected ballpark for true CCCH coordination:",
                    "  - Zn–S(Cys SG) ~ 2.2–2.4 Å",
                    "  - Zn–N(His NE2/ND1) ~ 2.0–2.2 Å",
                    "",
                    df_zn.to_markdown(index=False),
                ]
            )
        )

        # Zn finger split summary (% models in ZF1 vs ZF2) per complex
        split = (
            df_zn.groupby(["ComplexBase", "AssignedFinger"])["Model"]
            .nunique()
            .rename("n_models")
            .reset_index()
        )
        total = df_zn.groupby("ComplexBase")["Model"].nunique().rename("n_models_total")
        split = split.merge(total.reset_index(), on="ComplexBase", how="left")
        split["pct_models"] = 100.0 * split["n_models"] / split["n_models_total"]
        split.sort_values(["ComplexBase", "AssignedFinger"], inplace=True)
        split_csv = OUTPUT_DIR / "zinc_finger_split_summary.csv"
        split_md = OUTPUT_DIR / "zinc_finger_split_summary.md"
        split.to_csv(split_csv, index=False)
        split_md.write_text(
            "\n".join(
                [
                    "# Zinc runs: ZnFinger split summary",
                    "",
                    "- If a complex has substantial mass in both ZF1 and ZF2, do not average across them.",
                    "",
                    split.to_markdown(index=False),
                ]
            )
        )

    def make_avg_table(df_models_in: pd.DataFrame) -> pd.DataFrame:
        group_cols = ["ComplexBase", "Protein", "RNA", "ProteinExpClass", "RNAExpClass"]
        numeric_cols = [
            c
            for c in df_models_in.columns
            if c not in group_cols + ["pdb", "Notes"]
            and pd.api.types.is_numeric_dtype(df_models_in[c])
        ]
        numeric_cols = [c for c in numeric_cols if c != "Model"]
        grouped = df_models_in.groupby(group_cols, dropna=False)
        df_mean = grouped[numeric_cols].mean().reset_index()
        df_sd = grouped[numeric_cols].std(ddof=1).add_suffix("_sd").reset_index()
        df_out = df_mean.merge(df_sd, on=group_cols, how="left")
        if "ipTM" in df_models_in.columns:
            iptm_frac = grouped["ipTM"].apply(
                lambda s: float((s > IPTM_THRESHOLD).mean())
                if len(s)
                else float("nan")
            ).rename(f"ipTM_frac_gt_{IPTM_THRESHOLD}")
            df_out = df_out.merge(iptm_frac.reset_index(), on=group_cols, how="left")
            df_out[f"ipTM_pct_gt_{IPTM_THRESHOLD}"] = 100.0 * df_out[
                f"ipTM_frac_gt_{IPTM_THRESHOLD}"
            ]
            df_out.drop(columns=[f"ipTM_frac_gt_{IPTM_THRESHOLD}"], inplace=True)
        # Also compute % of models with protein↔RNA pair ipTM above threshold (pairwise interface)
        if "pair_ipTM_protein_rna" in df_models_in.columns:
            pair_frac = grouped["pair_ipTM_protein_rna"].apply(
                lambda s: float((s > IPTM_THRESHOLD).mean())
                if len(s)
                else float("nan")
            ).rename(f"pair_ipTM_protein_rna_frac_gt_{IPTM_THRESHOLD}")
            df_out = df_out.merge(pair_frac.reset_index(), on=group_cols, how="left")
            df_out[f"pair_ipTM_protein_rna_pct_gt_{IPTM_THRESHOLD}"] = 100.0 * df_out[
                f"pair_ipTM_protein_rna_frac_gt_{IPTM_THRESHOLD}"
            ]
            df_out.drop(
                columns=[f"pair_ipTM_protein_rna_frac_gt_{IPTM_THRESHOLD}"],
                inplace=True,
            )
        counts = (
            df_models_in.groupby("ComplexBase")["Model"]
            .count()
            .rename("n_models")
        )
        df_out = df_out.merge(counts.reset_index(), on="ComplexBase", how="left")
        df_out.sort_values(["ProteinExpClass", "Protein", "RNA"], inplace=True)
        return df_out

    def make_best_table(df_models_in: pd.DataFrame) -> pd.DataFrame:
        # pick best model per ComplexBase by ranking_score; fallback to ipTM; then model index
        df = df_models_in.copy()
        df["_rank"] = df["ranking_score"]
        if df["_rank"].isna().all():
            df["_rank"] = df["ipTM"]
        df.sort_values(
            ["ComplexBase", "_rank", "Model"],
            ascending=[True, False, True],
            inplace=True,
        )
        best = df.groupby("ComplexBase", as_index=False).head(1).copy()
        counts = (
            df_models_in.groupby("ComplexBase")["Model"]
            .count()
            .rename("n_models")
        )
        best = best.merge(counts.reset_index(), on="ComplexBase", how="left")
        best.rename(columns={"Model": "best_model"}, inplace=True)
        best.drop(columns=["_rank"], inplace=True)
        best.sort_values(["ProteinExpClass", "Protein", "RNA"], inplace=True)
        return best

    # Best-model outputs (zinc + non-zinc)
    df_best = make_best_table(df_models)
    csv_best = OUTPUT_DIR / "interface_metrics_best.csv"
    md_best = OUTPUT_DIR / "interface_metrics_best.md"
    df_best.to_csv(csv_best, index=False)
    md_best.write_text(
        "\n".join(
            [
                "# ZFP36L2 AlphaFold3 interface metrics (best model, zinc + non-zinc)",
                "",
                "- Best model chosen by max `ranking_score` (fallback: max ipTM)",
                "",
                df_best.to_markdown(index=False),
            ]
        )
    )

    # Averaged outputs (mean across models per complex) + SD + % above thresholds
    group_cols = ["ComplexBase", "Protein", "RNA", "ProteinExpClass", "RNAExpClass"]
    numeric_cols = [
        c
        for c in df_models.columns
        if c not in group_cols + ["pdb", "Notes"]
        and pd.api.types.is_numeric_dtype(df_models[c])
    ]
    # Don't average over model index; keep it out of the averaged table
    numeric_cols = [c for c in numeric_cols if c != "Model"]
    df_avg = make_avg_table(df_models)

    csv_avg = OUTPUT_DIR / "interface_metrics_avg.csv"
    md_avg = OUTPUT_DIR / "interface_metrics_avg.md"
    df_avg.to_csv(csv_avg, index=False)
    md_avg.write_text(
        "\n".join(
            [
                "# ZFP36L2 AlphaFold3 interface metrics (averaged over models)",
                "",
                f"- Complexes analyzed: {len(df_avg)}",
                "- Means are computed over available AF3 models per complex",
                f"- SDs are computed across models (ddof=1); ipTM threshold = {IPTM_THRESHOLD}",
                "- RMSD columns are in Å; BSA in Å²",
                "",
                df_avg.to_markdown(index=False),
            ]
        )
    )

    # Plots from averaged table (less cluttered than per-model)
    plot_paths = make_plots(df_avg)

    # Zinc-only outputs (folders/complexes ending with "_zinc")
    df_models_zinc = df_models[df_models["ComplexBase"].str.endswith("_zinc")].copy()
    if len(df_models_zinc):
        df_avg_zinc = make_avg_table(df_models_zinc)
        df_best_zinc = make_best_table(df_models_zinc)

        # If ZnFinger splits, compute stratified averages (ComplexBase + ZnFinger)
        if zinc_coord_rows:
            df_zn_map = pd.DataFrame(zinc_coord_rows)[
                ["ComplexBase", "Model", "AssignedFinger", "ZnLigandAtom"]
            ].drop_duplicates()
            df_models_zinc = df_models_zinc.merge(
                df_zn_map, on=["ComplexBase", "Model"], how="left"
            )
            df_models_zinc.rename(
                columns={"AssignedFinger": "ZnFinger"}, inplace=True
            )

            df_avg_zinc_by_finger = make_avg_table(
                df_models_zinc.drop(columns=["ZnFinger", "ZnLigandAtom"], errors="ignore")
            )
            # also write a stratified table (so you can avoid mixing when split)
            group_cols_zf = [
                "ComplexBase",
                "ZnFinger",
                "Protein",
                "RNA",
                "ProteinExpClass",
                "RNAExpClass",
            ]
            numeric_cols_zf = [
                c
                for c in df_models_zinc.columns
                if c not in group_cols_zf + ["pdb", "Notes", "ZnLigandAtom"]
                and pd.api.types.is_numeric_dtype(df_models_zinc[c])
            ]
            numeric_cols_zf = [c for c in numeric_cols_zf if c != "Model"]
            grouped_zf = df_models_zinc.groupby(group_cols_zf, dropna=False)
            df_zf_mean = grouped_zf[numeric_cols_zf].mean().reset_index()
            df_zf_sd = grouped_zf[numeric_cols_zf].std(ddof=1).add_suffix("_sd").reset_index()
            df_zf = df_zf_mean.merge(df_zf_sd, on=group_cols_zf, how="left")
            counts_zf = grouped_zf["Model"].count().rename("n_models")
            df_zf = df_zf.merge(counts_zf.reset_index(), on=group_cols_zf, how="left")
            df_zf.sort_values(["ComplexBase", "ZnFinger"], inplace=True)

            csv_zf = OUTPUT_DIR / "interface_metrics_zinc_avg_by_znfinger.csv"
            md_zf = OUTPUT_DIR / "interface_metrics_zinc_avg_by_znfinger.md"
            df_zf.to_csv(csv_zf, index=False)
            md_zf.write_text(
                "\n".join(
                    [
                        "# ZFP36L2 AlphaFold3 interface metrics (zinc runs, averaged by ZnFinger)",
                        "",
                        "- Use this when a complex has ZnFinger split across models (ZF1 vs ZF2).",
                        "",
                        df_zf.to_markdown(index=False),
                    ]
                )
            )

        csv_avg_zinc = OUTPUT_DIR / "interface_metrics_zinc_avg.csv"
        md_avg_zinc = OUTPUT_DIR / "interface_metrics_zinc_avg.md"
        df_avg_zinc.to_csv(csv_avg_zinc, index=False)
        md_avg_zinc.write_text(
            "\n".join(
                [
                    "# ZFP36L2 AlphaFold3 interface metrics (zinc runs, averaged over models)",
                    "",
                    f"- Complexes analyzed: {len(df_avg_zinc)}",
                    "- Means are computed over available AF3 models per complex",
                    f"- SDs are computed across models (ddof=1); ipTM threshold = {IPTM_THRESHOLD}",
                    "",
                    df_avg_zinc.to_markdown(index=False),
                ]
            )
        )

        csv_best_zinc = OUTPUT_DIR / "interface_metrics_zinc_best.csv"
        md_best_zinc = OUTPUT_DIR / "interface_metrics_zinc_best.md"
        df_best_zinc.to_csv(csv_best_zinc, index=False)
        md_best_zinc.write_text(
            "\n".join(
                [
                    "# ZFP36L2 AlphaFold3 interface metrics (zinc runs, best model)",
                    "",
                    "- Best model chosen by max `ranking_score` (fallback: max ipTM)",
                    "",
                    df_best_zinc.to_markdown(index=False),
                ]
            )
        )

    # 2Zn-only outputs (folders/complexes ending with "_2zn")
    df_models_2zn = df_models[df_models["ComplexBase"].str.endswith("_2zn")].copy()
    if len(df_models_2zn):
        df_avg_2zn = make_avg_table(df_models_2zn)
        df_best_2zn = make_best_table(df_models_2zn)

        csv_avg_2zn = OUTPUT_DIR / "interface_metrics_2zn_avg.csv"
        md_avg_2zn = OUTPUT_DIR / "interface_metrics_2zn_avg.md"
        df_avg_2zn.to_csv(csv_avg_2zn, index=False)
        md_avg_2zn.write_text(
            "\n".join(
                [
                    "# ZFP36L2 AlphaFold3 interface metrics (2Zn runs, averaged over models)",
                    "",
                    f"- Complexes analyzed: {len(df_avg_2zn)}",
                    "- Means are computed over available AF3 models per complex",
                    f"- SDs are computed across models (ddof=1); ipTM threshold = {IPTM_THRESHOLD}",
                    "",
                    df_avg_2zn.to_markdown(index=False),
                ]
            )
        )

        csv_best_2zn = OUTPUT_DIR / "interface_metrics_2zn_best.csv"
        md_best_2zn = OUTPUT_DIR / "interface_metrics_2zn_best.md"
        df_best_2zn.to_csv(csv_best_2zn, index=False)
        md_best_2zn.write_text(
            "\n".join(
                [
                    "# ZFP36L2 AlphaFold3 interface metrics (2Zn runs, best model)",
                    "",
                    "- Best model chosen by max `ranking_score` (fallback: max ipTM)",
                    "",
                    df_best_2zn.to_markdown(index=False),
                ]
            )
        )

        # 2Zn mapping + per-model partner metrics/QC tables
        # - Entity mapping: entity_idx -> type (protein/RNA/Zn1/Zn2) and PDB chain ID
        map_rows: List[Dict[str, object]] = []
        qc_rows: List[Dict[str, object]] = []
        for _, r in df_models_2zn.sort_values(["ComplexBase", "Model"]).iterrows():
            base = str(r["ComplexBase"])
            model = int(r["Model"])
            pdb_name = str(r["pdb"])
            pdb_path = OUTPUT_DIR / "converted_pdb" / pdb_name
            chains = pdb_chain_order(pdb_path)
            n_ent = int(r["n_entities"]) if not pd.isna(r["n_entities"]) else None
            pidx = int(r["protein_idx"]) if not pd.isna(r["protein_idx"]) else None
            ridx = int(r["rna_idx"]) if not pd.isna(r["rna_idx"]) else None
            ion1 = int(r["ion_idx"]) if not pd.isna(r.get("ion_idx", np.nan)) else None
            ion2 = int(r["ion_idx_2"]) if not pd.isna(r.get("ion_idx_2", np.nan)) else None

            idx_to_chain: Dict[int, Optional[str]] = {}
            if n_ent is not None and chains and len(chains) >= n_ent:
                for i in range(n_ent):
                    idx_to_chain[i] = chains[i]
            else:
                for i in range(n_ent or 0):
                    idx_to_chain[i] = None

            def etype(i: int) -> str:
                if pidx is not None and i == pidx:
                    return "protein"
                if ridx is not None and i == ridx:
                    return "RNA"
                if ion1 is not None and i == ion1:
                    return "Zn1"
                if ion2 is not None and i == ion2:
                    return "Zn2"
                return "other"

            for i in range(n_ent or 0):
                map_rows.append(
                    {
                        "ComplexBase": base,
                        "Model": model,
                        "entity_idx": i,
                        "entity_type": etype(i),
                        "pdb_chain_id": idx_to_chain.get(i, None),
                    }
                )

            # partner-metric QC: ensure cross mean is exactly (PAE0,3 + PAE3,0)/2 when both present
            pae_pr = r.get("pair_pae_min_protein_to_rna", np.nan)
            pae_rp = r.get("pair_pae_min_rna_to_protein", np.nan)
            expected = None
            if not pd.isna(pae_pr) and not pd.isna(pae_rp):
                expected = 0.5 * (float(pae_pr) + float(pae_rp))
            reported = r.get("pair_pae_min_cross_mean", np.nan)
            diff = None
            if expected is not None and not pd.isna(reported):
                diff = float(reported) - float(expected)
            qc_rows.append(
                {
                    "ComplexBase": base,
                    "Model": model,
                    "protein_idx": pidx,
                    "rna_idx": ridx,
                    "ion_idx": ion1,
                    "ion_idx_2": ion2,
                    "pair_ipTM_protein_rna": r.get("pair_ipTM_protein_rna", np.nan),
                    "pair_pae_min_protein_to_rna": pae_pr,
                    "pair_pae_min_rna_to_protein": pae_rp,
                    "pair_pae_min_cross_mean_reported": reported,
                    "pair_pae_min_cross_mean_expected": expected,
                    "pair_pae_min_cross_mean_diff": diff,
                }
            )

        df_map = pd.DataFrame(map_rows).sort_values(
            ["ComplexBase", "Model", "entity_idx"]
        )
        map_csv = OUTPUT_DIR / "entity_mapping_2zn_models.csv"
        map_md = OUTPUT_DIR / "entity_mapping_2zn_models.md"
        df_map.to_csv(map_csv, index=False)
        map_md.write_text(
            "\n".join(
                [
                    "# 2Zn runs: entity index → type and PDB chain ID (per model)",
                    "",
                    "- Expect: protein=0, Zn1=1, Zn2=2, RNA=3, and chain IDs A/B/C/D respectively.",
                    "",
                    df_map.to_markdown(index=False),
                ]
            )
        )

        df_qc = pd.DataFrame(qc_rows).sort_values(["ComplexBase", "Model"])
        qc_csv = OUTPUT_DIR / "protein_rna_pair_metrics_2zn_models.csv"
        qc_md = OUTPUT_DIR / "protein_rna_pair_metrics_2zn_models.md"
        df_qc.to_csv(qc_csv, index=False)
        qc_md.write_text(
            "\n".join(
                [
                    "# 2Zn runs: per-model protein↔RNA pair metrics (explicit [protein_idx][rna_idx])",
                    "",
                    f"- ipTM threshold for % summaries: {IPTM_THRESHOLD}",
                    "- `pair_pae_min_cross_mean_expected` is exactly (PAEprot→RNA + PAERNA→prot)/2.",
                    "",
                    df_qc.to_markdown(index=False),
                ]
            )
        )

        # Simple summary per complex: distributions + % models with pair ipTM > threshold
        summ_rows: List[Dict[str, object]] = []
        for base, sub in df_qc.groupby("ComplexBase"):
            s = pd.to_numeric(sub["pair_ipTM_protein_rna"], errors="coerce")
            summ_rows.append(
                {
                    "ComplexBase": base,
                    "n_models": int(len(sub)),
                    "rna_idx_mode": int(pd.Series(sub["rna_idx"]).mode().iloc[0])
                    if len(sub["rna_idx"]) and not pd.isna(pd.Series(sub["rna_idx"]).mode().iloc[0])
                    else None,
                    "pair_ipTM_mean": float(s.mean()) if not s.isna().all() else None,
                    "pair_ipTM_sd": float(s.std(ddof=1)) if s.notna().sum() > 1 else None,
                    "pair_ipTM_min": float(s.min()) if s.notna().sum() else None,
                    "pair_ipTM_max": float(s.max()) if s.notna().sum() else None,
                    f"pair_ipTM_pct_gt_{IPTM_THRESHOLD}": float((s > IPTM_THRESHOLD).mean() * 100.0)
                    if len(s)
                    else None,
                }
            )
        df_summ = pd.DataFrame(summ_rows).sort_values("ComplexBase")
        summ_csv = OUTPUT_DIR / "protein_rna_pair_metrics_2zn_summary.csv"
        summ_md = OUTPUT_DIR / "protein_rna_pair_metrics_2zn_summary.md"
        df_summ.to_csv(summ_csv, index=False)
        summ_md.write_text(
            "\n".join(
                [
                    "# 2Zn runs: protein↔RNA pair ipTM distribution summary",
                    "",
                    df_summ.to_markdown(index=False),
                ]
            )
        )

        # Zn1/Zn2 finger assignment stability (per model + summary)
        try:
            df_zn_all = pd.read_csv(OUTPUT_DIR / "zinc_coordination_models.csv")
            df_zn_2zn = df_zn_all[df_zn_all["ComplexBase"].str.endswith("_2zn")].copy()
        except Exception:
            df_zn_2zn = pd.DataFrame()

        if len(df_zn_2zn) and len(df_map):
            # build chain->entity_idx per (ComplexBase, Model)
            chain_to_idx = (
                df_map.pivot_table(
                    index=["ComplexBase", "Model"],
                    columns="pdb_chain_id",
                    values="entity_idx",
                    aggfunc="first",
                )
                .reset_index()
            )
            long = []
            for _, rr in chain_to_idx.iterrows():
                base = rr["ComplexBase"]
                model = rr["Model"]
                for ch in [c for c in chain_to_idx.columns if c not in ["ComplexBase", "Model"]]:
                    v = rr.get(ch, np.nan)
                    if pd.isna(v):
                        continue
                    long.append({"ComplexBase": base, "Model": model, "ZN_chain": ch, "ZnEntityIdx": int(v)})
            df_chain_idx = pd.DataFrame(long)
            df_zn_2zn = df_zn_2zn.merge(
                df_chain_idx, on=["ComplexBase", "Model", "ZN_chain"], how="left"
            )

            def zn_label(row):
                if row.get("ZnEntityIdx", None) == 1:
                    return "Zn1"
                if row.get("ZnEntityIdx", None) == 2:
                    return "Zn2"
                return "Zn?"

            df_zn_2zn["ZnLabel"] = df_zn_2zn.apply(zn_label, axis=1)
            df_zn_2zn.sort_values(["ComplexBase", "Model", "ZnLabel"], inplace=True)

            zn_csv = OUTPUT_DIR / "zinc_finger_assignment_2zn_models.csv"
            zn_md = OUTPUT_DIR / "zinc_finger_assignment_2zn_models.md"
            df_zn_2zn.to_csv(zn_csv, index=False)
            zn_md.write_text(
                "\n".join(
                    [
                        "# 2Zn runs: per-model Zn1/Zn2 finger assignment",
                        "",
                        "- Expect: Zn1→ZF1 and Zn2→ZF2 consistently across models.",
                        "",
                        df_zn_2zn.to_markdown(index=False),
                    ]
                )
            )

            stab = (
                df_zn_2zn.groupby(["ComplexBase", "ZnLabel", "AssignedFinger"])["Model"]
                .nunique()
                .rename("n_models")
                .reset_index()
            )
            total = (
                df_zn_2zn.groupby(["ComplexBase", "ZnLabel"])["Model"]
                .nunique()
                .rename("n_models_total")
                .reset_index()
            )
            stab = stab.merge(total, on=["ComplexBase", "ZnLabel"], how="left")
            stab["pct_models"] = 100.0 * stab["n_models"] / stab["n_models_total"]
            stab.sort_values(["ComplexBase", "ZnLabel", "AssignedFinger"], inplace=True)
            stab_csv = OUTPUT_DIR / "zinc_finger_assignment_2zn_summary.csv"
            stab_md = OUTPUT_DIR / "zinc_finger_assignment_2zn_summary.md"
            stab.to_csv(stab_csv, index=False)
            stab_md.write_text(
                "\n".join(
                    [
                        "# 2Zn runs: Zn1/Zn2 finger assignment stability summary",
                        "",
                        stab.to_markdown(index=False),
                    ]
                )
            )

    print(f"Wrote {csv_models}")
    print(f"Wrote {md_models}")
    print(f"Wrote {csv_best}")
    print(f"Wrote {md_best}")
    print(f"Wrote {csv_avg}")
    print(f"Wrote {md_avg}")
    if len(df_models_zinc):
        print(f"Wrote {csv_avg_zinc}")
        print(f"Wrote {md_avg_zinc}")
        print(f"Wrote {csv_best_zinc}")
        print(f"Wrote {md_best_zinc}")
    if len(df_models_2zn):
        print(f"Wrote {csv_avg_2zn}")
        print(f"Wrote {md_avg_2zn}")
        print(f"Wrote {csv_best_2zn}")
        print(f"Wrote {md_best_2zn}")
        # additional 2Zn QC outputs (if created)
        for extra in [
            OUTPUT_DIR / "entity_mapping_2zn_models.csv",
            OUTPUT_DIR / "entity_mapping_2zn_models.md",
            OUTPUT_DIR / "protein_rna_pair_metrics_2zn_models.csv",
            OUTPUT_DIR / "protein_rna_pair_metrics_2zn_models.md",
            OUTPUT_DIR / "protein_rna_pair_metrics_2zn_summary.csv",
            OUTPUT_DIR / "protein_rna_pair_metrics_2zn_summary.md",
            OUTPUT_DIR / "zinc_finger_assignment_2zn_models.csv",
            OUTPUT_DIR / "zinc_finger_assignment_2zn_models.md",
            OUTPUT_DIR / "zinc_finger_assignment_2zn_summary.csv",
            OUTPUT_DIR / "zinc_finger_assignment_2zn_summary.md",
        ]:
            if extra.exists():
                print(f"Wrote {extra}")
        if zinc_coord_rows:
            split_csv = OUTPUT_DIR / "zinc_finger_split_summary.csv"
            split_md = OUTPUT_DIR / "zinc_finger_split_summary.md"
            if split_csv.exists():
                print(f"Wrote {split_csv}")
            if split_md.exists():
                print(f"Wrote {split_md}")
            csv_zf = OUTPUT_DIR / "interface_metrics_zinc_avg_by_znfinger.csv"
            md_zf = OUTPUT_DIR / "interface_metrics_zinc_avg_by_znfinger.md"
            if csv_zf.exists():
                print(f"Wrote {csv_zf}")
            if md_zf.exists():
                print(f"Wrote {md_zf}")
        if entity_rows:
            print(f"Wrote {ent_csv}")
            print(f"Wrote {ent_md}")
        if zinc_coord_rows:
            print(f"Wrote {zn_csv}")
            print(f"Wrote {zn_md}")
    for path in plot_paths:
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
