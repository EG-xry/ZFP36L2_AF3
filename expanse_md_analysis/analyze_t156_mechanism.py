#!/usr/bin/env python3
"""Mechanistic capstone for the T156 RYKTEL-motif backbone-reading model (zero GPU).

The discriminating feature is the T156-RNA hydrogen bond, made ~96-100% by the T156
MAIN CHAIN (N-H, C=O) - the canonical Hudson-2004 backbone recognition. This script,
per variant, computes:

  Part A - which ARE base(s) the T156 backbone reads (per-base H-bond occupancy),
           and whether the deleterious loss maps to a specific base.
  Part B - the per-residue protein-RNA H-bond profile and its change vs WT, to test
           whether distal deleterious mutations disrupt the finger-1 RYKTEL region
           (resSeq 3-8: R153-Y154-K155-T156-E157-L158), i.e. long-range propagation.

Env: /Users/eric/micromamba/envs/omm/bin/python
"""
from pathlib import Path
import numpy as np
import pandas as pd
import mdtraj as md
from diagnostic_f156_convergence import (
    load_traj, SOURCES, CLASS, TRUTH_DEL, NEU, MIN_NS, ALL_VARIANTS, RNA_RES, NS_PER_FRAME)

OUT = Path(__file__).resolve().parent
T156_RESSEQ = 6
BACKBONE = {"N", "O"}                 # H-bonding main-chain atoms of T156
RYKTEL = list(range(3, 9))            # resSeq 3-8 = R153 Y154 K155 T156 E157 L158
MUT_SITE = {"Y154H": 4, "R160W": 10, "R184C": 34, "G204D": 54, "C206F": 56,
            "S185N": 35, "T187S": 37, "D210E": 60}   # resSeq of the mutation


def analyze_traj(traj):
    top = traj.topology
    atom_res = {a.index: a.residue for a in top.atoms}
    rna_objs = sorted([r for r in top.residues if r.name.upper() in RNA_RES],
                      key=lambda r: (r.chain.index, r.resSeq))
    rna_ord = {r.index: (i + 1, r.name.upper()[0]) for i, r in enumerate(rna_objs)}   # 1..9, base letter
    rna_atom = {a.index for a in top.atoms if a.residue.index in rna_ord}
    t156_bb = {a.index for a in top.atoms
               if a.residue.is_protein and a.residue.resSeq == T156_RESSEQ and a.name in BACKBONE}
    n = traj.n_frames
    base_hit = {i: 0 for i in range(1, len(rna_objs) + 1)}     # T156-bb occupancy per base (per-frame)
    perres = {}                                                # protein resSeq -> total H-bond count
    for fi in range(n):
        hbs = md.baker_hubbard(traj[fi], freq=0.0, exclude_water=True)
        bh = set()
        for d, h, a in hbs:
            # Part B: any protein residue <-> RNA
            if d in rna_atom and atom_res[a].is_protein:
                perres[atom_res[a].resSeq] = perres.get(atom_res[a].resSeq, 0) + 1
            elif a in rna_atom and atom_res[d].is_protein:
                perres[atom_res[d].resSeq] = perres.get(atom_res[d].resSeq, 0) + 1
            # Part A: T156 backbone <-> which RNA base
            if d in t156_bb and a in rna_atom:
                bh.add(rna_ord[atom_res[a].index][0])
            elif a in t156_bb and d in rna_atom:
                bh.add(rna_ord[atom_res[d].index][0])
        for b in bh:
            base_hit[b] += 1
    base_occ = {b: base_hit[b] / n for b in base_hit}
    perres_occ = {rs: perres.get(rs, 0) / n for rs in range(1, 71)}
    base_letters = {i: rna_ord[r.index][1] for i, r in enumerate(rna_objs, 1)}
    return base_occ, perres_occ, base_letters


def main():
    baseA, perresB, letters = {}, {}, None
    for v in ALL_VARIANTS:
        for src, base in SOURCES:
            traj = load_traj(base, v)
            if traj is None:
                continue
            ns = traj.n_frames * NS_PER_FRAME
            if ns < MIN_NS:
                continue
            bo, pr, lets = analyze_traj(traj)
            letters = letters or lets
            baseA.setdefault(v, []).append(bo)
            perresB.setdefault(v, []).append(pr)
            print(f"  {v:6s} {src:13s} {ns:.0f}ns  done", flush=True)

    # average across replicates
    baseA = {v: pd.DataFrame(rows).mean().to_dict() for v, rows in baseA.items()}
    perresB = {v: pd.DataFrame(rows).mean().to_dict() for v, rows in perresB.items()}

    # ---------- Part A: which base does T156 backbone read? ----------
    print("\n" + "=" * 70)
    print("PART A - T156 BACKBONE H-bond occupancy per ARE base (RNA = 5'-UUAUUUAUU-3')")
    print("=" * 70)
    bases = sorted(letters)
    hdr = "  variant  " + " ".join(f"{letters[b]}{b}" for b in bases)
    print(hdr)
    def classavg(cls):
        vs = [v for v in baseA if CLASS.get(v) == cls]
        return {b: np.mean([baseA[v][b] for v in vs]) for b in bases} if vs else None
    for label, data in [("WT", baseA.get("WT")),
                        ("del(avg)", classavg("deleterious") if False else
                         ({b: np.mean([baseA[v][b] for v in TRUTH_DEL if v in baseA]) for b in bases})),
                        ("neu(avg)", {b: np.mean([baseA[v][b] for v in NEU if v in baseA]) for b in bases})]:
        if data:
            print(f"  {label:8s}  " + " ".join(f"{data[b]:.2f}" for b in bases))
    pd.DataFrame(baseA).T.to_csv(OUT / "t156_base_reading.csv")

    # ---------- Part B: does deleterious disrupt the RYKTEL region? ----------
    print("\n" + "=" * 70)
    print("PART B - RYKTEL-region (resSeq 3-8) RNA H-bond contact vs mutation location")
    print("=" * 70)
    wt = perresB.get("WT", {})
    ryktel_wt = sum(wt.get(rs, 0) for rs in RYKTEL)
    print(f"  WT RYKTEL-region total RNA H-bonds/frame = {ryktel_wt:.2f}")
    print(f"  {'variant':8s} {'class':4s} {'RYKTEL_sum':>10s} {'Δvs WT':>8s} {'mut@resSeq':>11s} {'mut_local_Δ':>12s}")
    rows = []
    for v in sorted(perresB, key=lambda x: sum(perresB[x].get(rs, 0) for rs in RYKTEL)):
        pr = perresB[v]
        ryk = sum(pr.get(rs, 0) for rs in RYKTEL)
        ms = MUT_SITE.get(v)
        mut_local_d = (pr.get(ms, 0) - wt.get(ms, 0)) if ms else np.nan
        rows.append({"variant": v, "class": CLASS.get(v), "RYKTEL_sum": round(ryk, 3),
                     "d_vs_WT": round(ryk - ryktel_wt, 3),
                     "mut_resSeq": ms, "mut_local_d": round(mut_local_d, 3) if ms else None})
        print(f"  {v:8s} {str(CLASS.get(v)):4s} {ryk:10.2f} {ryk-ryktel_wt:+8.2f} "
              f"{str(ms):>11s} {('%+.2f'%mut_local_d) if ms else '   n/a':>12s}")
    dfB = pd.DataFrame(rows); dfB.to_csv(OUT / "t156_ryktel_propagation.csv", index=False)
    pd.DataFrame(perresB).to_csv(OUT / "perresidue_rna_hbond_profile.csv")

    # test: is RYKTEL-region contact lower in deleterious than neutral?
    dels = [r["RYKTEL_sum"] for r in rows if r["class"] == "deleterious"]
    neus = [r["RYKTEL_sum"] for r in rows if r["class"] == "neutral"]
    if dels and neus:
        print(f"\n  RYKTEL-region contact: del mean={np.mean(dels):.2f}  neu mean={np.mean(neus):.2f}  "
              f"({'del < neu (propagation supported)' if np.mean(dels) < np.mean(neus) else 'no clear separation'})")
    print(f"\nwrote t156_base_reading.csv, t156_ryktel_propagation.csv, perresidue_rna_hbond_profile.csv")


if __name__ == "__main__":
    main()
