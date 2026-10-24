#!/usr/bin/env python3
"""
Utility to normalize AF3 protein-RNA models for OpenMM.

Actions performed per input PDB:
- keep exactly one protein chain (renamed to A) and one RNA chain (renamed to B)
- drop waters/ions/unknown cofactors
- ensure (and optionally add) Zn2+ sites for CCCH motifs
- strip stray TER/END records and rebuild a clean record order
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Sequence, Tuple

PROTEIN_RESNAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ASX",
    "GLX",
}

RNA_RESNAMES = {
    "A",
    "C",
    "G",
    "U",
    "RA",
    "RC",
    "RG",
    "RU",
    "DA",
    "DC",
    "DG",
    "DT",
    "ADE",
    "CYT",
    "GUA",
    "URI",
    "URA",
}

WATER_RESNAMES = {"HOH", "H2O", "WAT", "TIP", "SOL"}
ION_RESNAMES = {"NA", "CL", "K", "MG", "CA", "MN", "CO", "CU", "FE", "IOD", "I"}
METAL_RESNAMES = {"ZN"}

RNA_TERMINI = {
    "A": {"5": "A5", "3": "A3"},
    "C": {"5": "C5", "3": "C3"},
    "G": {"5": "G5", "3": "G3"},
    "U": {"5": "U5", "3": "U3"},
}


@dataclass
class AtomRecord:
    raw: str
    resname: str
    atom: str
    resseq_field: str
    resseq_int: int | None
    resseq_label: str
    coords: Tuple[float, float, float]

    def format_line(
        self,
        serial: int,
        chain: str,
        resname_override: str | None = None,
        resseq_override: int | None = None,
    ) -> str:
        base = self.raw.rstrip("\n").ljust(80)
        base = f"{base[:6]}{serial:5d}{base[11:]}"
        base = f"{base[:21]}{chain}{base[22:]}"
        if resname_override is not None:
            base = f"{base[:17]}{resname_override:>3}{base[20:]}"
        if resseq_override is not None:
            base = f"{base[:22]}{resseq_override:4d}{base[26:]}"
        return base[:80] + "\n"


def parse_atom_line(line: str) -> AtomRecord:
    resname = line[17:20].strip().upper()
    atom = line[12:16].strip().upper()
    resseq_field = line[22:26]
    resseq_int = None
    stripped = resseq_field.strip()
    if stripped:
        try:
            resseq_int = int(stripped)
        except ValueError:
            pass
    coords = (
        float(line[30:38]),
        float(line[38:46]),
        float(line[46:54]),
    )
    return AtomRecord(
        raw=line.rstrip("\n"),
        resname=resname,
        atom=atom,
        resseq_field=resseq_field,
        resseq_int=resseq_int,
        resseq_label=resseq_field.strip(),
        coords=coords,
    )


def _canonical_rna_letter(resname: str) -> str | None:
    upper = resname.upper()
    if not upper:
        return None
    if upper in RNA_TERMINI:
        return upper
    if len(upper) > 1 and upper[-1] in RNA_TERMINI:
        return upper[-1]
    return None


def _terminal_resname(resname: str, terminus: str) -> str | None:
    letter = _canonical_rna_letter(resname)
    if letter is None:
        return None
    return RNA_TERMINI[letter][terminus]


def categorize_resname(resname: str) -> str | None:
    if resname in PROTEIN_RESNAMES:
        return "protein"
    if resname in RNA_RESNAMES:
        return "rna"
    if resname in METAL_RESNAMES:
        return "metal"
    if resname in WATER_RESNAMES or resname in ION_RESNAMES:
        return "skip"
    return None


def distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))


def detect_zinc_sites(records: Sequence[AtomRecord]) -> List[Tuple[Tuple[float, float, float], List[int]]]:
    cys_sg: dict[int, Tuple[float, float, float]] = {}
    his_coords: dict[int, Tuple[float, float, float]] = {}
    for rec in records:
        if rec.resseq_int is None:
            continue
        if rec.resname == "CYS" and rec.atom == "SG":
            cys_sg[rec.resseq_int] = rec.coords
        elif rec.resname == "HIS":
            if rec.atom == "NE2":
                his_coords[rec.resseq_int] = rec.coords
            elif rec.atom == "ND1" and rec.resseq_int not in his_coords:
                his_coords[rec.resseq_int] = rec.coords

    sites: List[Tuple[Tuple[float, float, float], List[int]]] = []
    used_cys: set[int] = set()
    for his_resseq in sorted(his_coords):
        hcoord = his_coords[his_resseq]
        nearby = [
            (distance(hcoord, ccoord), cseq, ccoord)
            for cseq, ccoord in cys_sg.items()
            if cseq not in used_cys and distance(hcoord, ccoord) <= 4.0
        ]
        if len(nearby) < 3:
            continue
        nearby.sort()
        selected = nearby[:3]
        used_cys.update(seq for _, seq, _ in selected)
        lig_coords = [hcoord] + [coord for _, _, coord in selected]
        centroid = tuple(sum(vals) / len(lig_coords) for vals in zip(*lig_coords))
        sites.append((centroid, [seq for _, seq, _ in selected]))
    return sites


def format_new_zn_line(serial: int, chain: str, resseq: int, coord: Tuple[float, float, float]) -> str:
    x, y, z = coord
    record = "HETATM"
    atom_name = "ZN"
    resname = "ZN"
    altloc = " "
    icode = " "
    element = "ZN"
    charge = "  "
    return (
        f"{record:<6}{serial:5d} {atom_name:>4}{altloc}"
        f"{resname:>3} {chain}{resseq:4d}{icode}   "
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
        f"{1.00:6.2f}{0.00:6.2f}          {element:>2}{charge}\n"
    )


def format_ter_line(serial: int, chain: str) -> str:
    return f"TER   {serial:5d}      {chain}\n"


def clean_pdb(input_path: Path, output_path: Path) -> None:
    header_lines: List[str] = []
    protein_records: List[AtomRecord] = []
    rna_records: List[AtomRecord] = []
    metal_records: List[AtomRecord] = []
    skipped: List[str] = []

    with input_path.open() as handle:
        for raw_line in handle:
            if raw_line.startswith(("ATOM", "HETATM")):
                rec = parse_atom_line(raw_line)
                category = categorize_resname(rec.resname)
                if category == "protein":
                    protein_records.append(rec)
                elif category == "rna":
                    rna_records.append(rec)
                elif category == "metal":
                    metal_records.append(rec)
                elif category == "skip":
                    continue
                else:
                    skipped.append(rec.resname)
                continue
            if raw_line.startswith(("TER", "END", "MODEL", "ENDMDL", "CONECT")):
                continue
            if raw_line.startswith(("HEADER", "TITLE", "REMARK", "COMPND", "SOURCE", "CRYST1")):
                header_lines.append(raw_line.rstrip("\n"))

    if not protein_records:
        raise ValueError(f"{input_path.name}: no protein atoms detected")
    if not rna_records:
        raise ValueError(f"{input_path.name}: no RNA atoms detected")

    rna_residue_order: List[Tuple[int | None, str, str]] = []
    seen_rna: set[Tuple[int | None, str]] = set()
    for rec in rna_records:
        key = (rec.resseq_int, rec.resseq_label)
        if key in seen_rna:
            continue
        seen_rna.add(key)
        rna_residue_order.append((rec.resseq_int, rec.resseq_label, rec.resname))

    rna_residue_overrides: dict[Tuple[int | None, str], str] = {}
    rna_first_key: Tuple[int | None, str] | None = None
    rna_last_key: Tuple[int | None, str] | None = None
    if rna_residue_order:
        rna_first_key = (rna_residue_order[0][0], rna_residue_order[0][1])
        rna_last_key = (rna_residue_order[-1][0], rna_residue_order[-1][1])
        first = rna_residue_order[0]
        first_override = _terminal_resname(first[2], "5")
        if first_override:
            rna_residue_overrides[(first[0], first[1])] = first_override
        last = rna_residue_order[-1]
        last_override = _terminal_resname(last[2], "3")
        if last_override:
            rna_residue_overrides[(last[0], last[1])] = last_override

    existing_zn_coords: List[Tuple[float, float, float]] = [rec.coords for rec in metal_records]
    inferred_zn_sites: List[Tuple[Tuple[float, float, float], List[int]]] = []
    if not existing_zn_coords:
        inferred_zn_sites = detect_zinc_sites(protein_records)

    # Normalize Zn records regardless of source so element/residue columns are
    # strictly PDB-aligned (some AF3 exports shift the element field).
    zn_coords: List[Tuple[float, float, float]] = []
    if existing_zn_coords:
        zn_coords.extend(existing_zn_coords)
    else:
        zn_coords.extend(coord for coord, _ in inferred_zn_sites)

    protein_resseqs = [rec.resseq_int for rec in protein_records if rec.resseq_int is not None]
    next_resseq = (max(protein_resseqs) if protein_resseqs else 0) + 1
    zn_resseqs = list(range(next_resseq, next_resseq + len(zn_coords)))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    serial = 1
    with output_path.open("w") as out:
        for line in header_lines:
            out.write(line + "\n")
        out.write(
            "REMARK 999 Cleaned with clean_af3_structures.py "
            "(protein=A, RNA=B, Zn enforced)\n"
        )
        for rec in protein_records:
            out.write(rec.format_line(serial, "A"))
            serial += 1
        for idx, coord in enumerate(zn_coords):
            out.write(format_new_zn_line(serial, "A", zn_resseqs[idx], coord))
            serial += 1
        out.write(format_ter_line(serial, "A"))
        serial += 1
        first_skip_atoms = {"P", "OP1", "OP2", "OP3"}
        last_skip_atoms = {"OP3"}
        for rec in rna_records:
            override = rna_residue_overrides.get((rec.resseq_int, rec.resseq_label))
            key = (rec.resseq_int, rec.resseq_label)
            if key == rna_first_key and rec.atom in first_skip_atoms:
                continue
            if key == rna_last_key and rec.atom in last_skip_atoms:
                continue
            out.write(rec.format_line(serial, "B", resname_override=override))
            serial += 1
        out.write(format_ter_line(serial, "B"))
        serial += 1
        out.write("END\n")

    info = (
        f"{input_path.name} -> {output_path.name}: "
        f"{len(protein_records)} protein atoms, {len(rna_records)} RNA atoms"
    )
    if inferred_zn_sites:
        zn_detail = ", ".join(
            f"Cys{','.join(str(res) for res in sorted(res_list))}"
            for _, res_list in inferred_zn_sites
        )
        info += f", added {len(inferred_zn_sites)} Zn (~{zn_detail})"
    if metal_records:
        info += f", normalized {len(metal_records)} existing Zn atoms"
    if skipped:
        info += f", skipped residues: {sorted(set(skipped))}"
    print(info)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="directory to place cleaned structures",
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="raw AF3 PDB files",
    )
    args = parser.parse_args(argv)
    for pdb_path in args.inputs:
        output_path = args.output_dir / f"{pdb_path.stem}_cleaned.pdb"
        clean_pdb(pdb_path, output_path)


if __name__ == "__main__":
    main()
