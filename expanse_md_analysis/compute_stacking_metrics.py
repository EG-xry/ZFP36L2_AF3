#!/usr/bin/env python3
"""Geometric aromatic pi-stacking metric for the REAL aromatic residues vs RNA (zero GPU).

Context: the AUROC-winning "F156" feature is actually T156 (RYKTEL-motif H-bond, NOT
stacking; verified from sequence). This script tests whether GENUINE aromatic stacking
discriminates del/neu by targeting the actual aromatic residues:
  Y154 (RYKTEL Tyr) + F162 F176 F180 F200 F205 F214,
plus an unbiased scan of ALL aromatics in WT to find real RNA-stackers.

Per-frame stacking = ring-plane centroid distance < D_MAX AND inter-plane angle < ANGLE_MAX,
using the TDD-validated stacking_geometry core. Baseline to beat: T156 H-bond AUROC 0.800.

Env: /Users/eric/micromamba/envs/omm/bin/python
"""
from pathlib import Path
import numpy as np
import pandas as pd
from stacking_geometry import stacking_series
from diagnostic_f156_convergence import (
    load_traj, SOURCES, CLASS, TRUTH_DEL, NEU, MIN_NS, ALL_VARIANTS, auroc,
    RNA_RES, NS_PER_FRAME)

OUT = Path(__file__).resolve().parent

AROM = {
    "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"], "HID": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "HIE": ["CG", "ND1", "CD2", "CE1", "NE2"], "HIP": ["CG", "ND1", "CD2", "CE1", "NE2"],
    "TRP": ["CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"],
}
PURINE = ["N9", "C8", "N7", "C5", "C4", "N3", "C2", "N1", "C6"]
PYRIMIDINE = ["N1", "C2", "N3", "C4", "C5", "C6"]

# The genuine aromatics to test (ZFP36L2 label -> resSeq = num-150). Y154 = RYKTEL Tyr.
# Y170,Y208 added after the unbiased WT scan flagged them as real RNA-stackers.
POS = {"Y154": 4, "F162": 12, "Y170": 20, "F176": 26, "F180": 30,
       "F200": 50, "F205": 55, "Y208": 58, "F214": 64}
FINGER1 = [k for k, v in POS.items() if v <= 34]   # Y154,F162,F176,F180
FINGER2 = [k for k, v in POS.items() if v > 34]     # F200,F205,F214
D_MAX, ANGLE_MAX = 5.0, 40.0

_sanity_done = False
_warned = set()


def prot_ring(top, resSeq):
    for r in top.residues:
        if r.is_protein and r.resSeq == resSeq:
            names = AROM.get(r.name.upper())
            if not names:
                return None, r.name
            idx = [a.index for a in r.atoms if a.name in names]
            return (idx if len(idx) >= 3 else None), r.name
    return None, None


def base_ring_idx_lists(top):
    lists = []
    for r in top.residues:
        rn = r.name.upper()
        if rn in RNA_RES:
            names = PURINE if rn[0] in ("A", "G") else PYRIMIDINE
            idx = [a.index for a in r.atoms if a.name in names]
            if len(idx) >= 3:
                lists.append(idx)
    return lists


def occ_for_traj(traj, variant):
    global _sanity_done
    xyz = traj.xyz * 10.0
    top = traj.topology
    bases = [xyz[:, bi, :] for bi in base_ring_idx_lists(top)]
    if not bases:
        return None
    if not _sanity_done:
        print(f"  [sanity] resSeq4={prot_ring(top,4)[1]}  RNA_bases={len(bases)}  "
              f"testing aromatics: {list(POS)}")
        _sanity_done = True

    series = {}
    out = {}
    for lab, rs in POS.items():
        ridx, rname = prot_ring(top, rs)
        if ridx is None:
            series[lab] = np.zeros(traj.n_frames, bool)
            if (variant, lab) not in _warned:
                print(f"    [note] {variant}: {lab} (resSeq{rs})={rname} non-aromatic -> occ=0")
                _warned.add((variant, lab))
        else:
            series[lab] = stacking_series(xyz[:, ridx, :], bases, D_MAX, ANGLE_MAX)
        out[lab] = float(series[lab].mean())
    S = np.array([series[l] for l in POS])                       # (nres, F)
    out["aromatic_any"] = float(np.any(S, axis=0).mean())        # any of the 7 stacks
    out["aromatic_sum"] = float(S.sum(axis=0).mean())            # mean # stacked at once
    out["finger1_any"] = float(np.any([series[l] for l in FINGER1], axis=0).mean())
    out["finger2_any"] = float(np.any([series[l] for l in FINGER2], axis=0).mean())
    return out


def feature_auroc(valid_df, feat):
    avg = valid_df.groupby("variant")[feat].mean()
    if not all(v in avg.index for v in TRUTH_DEL | NEU):
        return None, None
    vs = [v for v in avg.index if v in TRUTH_DEL or v in NEU]
    a = auroc(np.array([-avg[v] for v in vs]),
              np.array([1 if v in TRUTH_DEL else 0 for v in vs]))
    dl = {}
    for v in avg.index:
        m = np.sort(valid_df[valid_df.variant == v][feat].values)
        dl[v] = m[1:].mean() if len(m) > 1 else m.mean()
    vs2 = [v for v in dl if v in TRUTH_DEL or v in NEU]
    a2 = auroc(np.array([-dl[v] for v in vs2]),
               np.array([1 if v in TRUTH_DEL else 0 for v in vs2]))
    return a, a2


def scan_all_aromatics_wt():
    """Unbiased: which aromatic residues actually stack with RNA in WT?"""
    traj = load_traj(SOURCES[0][1], "WT")
    if traj is None:
        return
    xyz = traj.xyz * 10.0
    top = traj.topology
    bases = [xyz[:, bi, :] for bi in base_ring_idx_lists(top)]
    rows = []
    for r in top.residues:
        if r.is_protein and r.name.upper() in AROM:
            idx = [a.index for a in r.atoms if a.name in AROM[r.name.upper()]]
            if len(idx) >= 3:
                occ = float(stacking_series(xyz[:, idx, :], bases, D_MAX, ANGLE_MAX).mean())
                rows.append((r.resSeq + 150, r.name, occ))
    rows.sort(key=lambda x: -x[2])
    print("\n" + "=" * 60)
    print("UNBIASED WT SCAN: stacking occupancy of ALL aromatic residues")
    print("  (which residues actually pi-stack with RNA in WT?)")
    print("=" * 60)
    for zfp, nm, occ in rows:
        star = "  <-- STACKS" if occ > 0.15 else ""
        print(f"    {nm}{zfp:<4d}  occ={occ:.3f}{star}")


def main():
    rows = []
    for v in ALL_VARIANTS:
        for src, base in SOURCES:
            traj = load_traj(base, v)
            if traj is None:
                continue
            ns = traj.n_frames * NS_PER_FRAME
            o = occ_for_traj(traj, v)
            if o is None:
                continue
            o.update(variant=v, cls=CLASS.get(v, "?"), source=src, ns=round(ns, 1))
            rows.append(o)
            print(f"  {v:6s} {src:13s} {ns:5.0f}ns  any={o['aromatic_any']:.2f} "
                  f"Y154={o['Y154']:.2f} f1={o['finger1_any']:.2f} f2={o['finger2_any']:.2f}"
                  + ("  [PARTIAL]" if ns < MIN_NS else ""), flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "stacking_metrics.csv", index=False)
    valid = df[df["ns"] >= MIN_NS]
    wt_occ = valid[valid.variant == "WT"][list(POS)].mean()

    feats = list(POS) + ["aromatic_any", "aromatic_sum", "finger1_any", "finger2_any"]
    print("\n" + "=" * 74)
    print("AROMATIC-STACKING AUROC (lower stacking = deleterious; 5 del vs 3 neu, >=180ns)")
    print("  baseline: T156 RYKTEL H-bond AUROC = 0.800 (drop-lowest 0.800)")
    print("=" * 74)
    for f in feats:
        a, a2 = feature_auroc(valid, f)
        wt = f"  WT_occ={wt_occ[f]:.2f}" if f in wt_occ.index else ""
        if a is None:
            print(f"  {f:14s}  (missing variants)")
        else:
            flag = "  <<< beats/matches T156" if a >= 0.80 else ""
            print(f"  {f:14s}  AUROC={a:.3f}  drop-lowest={a2:.3f}{wt}{flag}")

    scan_all_aromatics_wt()
    print(f"\nwrote {OUT / 'stacking_metrics.csv'}")


if __name__ == "__main__":
    main()
