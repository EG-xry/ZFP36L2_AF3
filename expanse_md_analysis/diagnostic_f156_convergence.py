#!/usr/bin/env python3
"""
Make-or-break diagnostic: is the F156-RNA contact signal REAL or a SAMPLING ARTIFACT?

Interrogates EXISTING trajectories (zero GPU). For each (source, variant) it computes
the per-frame F156(resSeq 6)<->RNA H-bond count (Baker-Hubbard, matching
multi_replicate_auroc_v2.py) at 1 ns sampling, and reports per trajectory:
  - full mean (AUROC input) | 1st vs 2nd half (drift) | frac_zero | blocking SEM | tau_int | n_eff
Per variant it prints replicate spread (valid >=180 ns reps only) and the leverage of the
lowest replicate, then RECOMPUTES the F156 AUROC (5 del vs 3 neu) under this clean sampling
plus a drop-lowest-replicate sensitivity AUROC.

Env: /Users/eric/micromamba/envs/omm/bin/python
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import mdtraj as md

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
PLOTS = OUT / "plots"; PLOTS.mkdir(exist_ok=True)

RNA_RES = {"A", "C", "G", "U", "A5", "A3", "C5", "C3", "G5", "G3", "U5", "U3"}
F156_PDB_RESSEQ = 6
SOURCES = [
    ("expanse_orig", ROOT / "expanse_md_results"),
    ("anvil_rep1",   ROOT / "anvil_md_results" / "rep1"),
    ("anvil_rep2",   ROOT / "anvil_md_results" / "rep2"),
    ("expanse_rep2", ROOT / "expanse_md_results_rep2"),
]
CLASS = {"WT": "WT", "C206F": "del", "C212A": "zn_ctrl", "G204D": "del",
         "R160W": "del", "R184C": "del", "Y154H": "del",
         "D210E": "neu", "S185N": "neu", "T187S": "neu"}
TRUTH_DEL = {"Y154H", "R160W", "R184C", "G204D", "C206F"}
NEU = {"D210E", "S185N", "T187S"}
ALL_VARIANTS = ["WT", "Y154H", "R160W", "R184C", "G204D", "C206F",
                "C212A", "D210E", "S185N", "T187S"]

LOAD_STRIDE = 10
NS_PER_FRAME = 0.1 * LOAD_STRIDE   # 1.0 ns
MAX_NS = 200.0
MIN_NS = 180.0   # runs shorter than this are "partial" -> excluded from variant averages (matches AUROC dataset)


def rna_sel():
    return " or ".join(f"resname {r}" for r in RNA_RES)


def load_traj(base_dir: Path, variant: str):
    for vname in (["C212A", "C208A"] if variant == "C212A" else [variant]):
        d = base_dir / vname
        top, dcd = d / "equilibrated.pdb", d / "production.dcd"
        if top.exists() and dcd.exists():
            break
    else:
        return None
    try:
        full_top = md.load_topology(str(top))
        keep = np.unique(np.concatenate([
            full_top.select("protein"), full_top.select(rna_sel()),
            full_top.select("name ZN or element Zn")]))
        traj = md.load(str(dcd), top=str(top), atom_indices=keep, stride=LOAD_STRIDE)
    except Exception as e:
        print(f"    load error {variant}/{base_dir.name}: {e}")
        return None
    nmax = int(MAX_NS / NS_PER_FRAME) + 1
    return traj[:nmax] if traj.n_frames > nmax else traj


def f156_hbond_series(traj) -> np.ndarray:
    top = traj.topology
    rna_atoms = {a.index for a in top.atoms if a.residue.name.upper() in RNA_RES}
    f156 = {a.index for a in top.atoms
            if a.residue.is_protein and a.residue.resSeq == F156_PDB_RESSEQ}
    if not f156:
        return np.full(traj.n_frames, np.nan)
    out = np.zeros(traj.n_frames)
    for fi in range(traj.n_frames):
        try:
            hbs = md.baker_hubbard(traj[fi], freq=0.0, exclude_water=True)
        except Exception:
            continue
        out[fi] = sum(1 for d, h, a in hbs
                      if (d in rna_atoms and a in f156) or (a in rna_atoms and d in f156))
    return out


def integrated_autocorr_time(x, c=5.0):
    x = np.asarray(x, float); n = len(x); xm = x - x.mean()
    if np.allclose(xm, 0) or n < 4:
        return 0.5
    f = np.fft.fft(xm, n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[:n].real
    if acf[0] == 0:
        return 0.5
    acf /= acf[0]; tau = 0.5
    for t in range(1, n):
        tau += acf[t]
        if t >= c * tau:
            break
    return max(tau, 0.5)


def blocking_sem(x):
    x = np.asarray(x, float).copy()
    if len(x) < 4 or np.allclose(x, x[0]):
        return 0.0
    sems = []
    while len(x) >= 4:
        sems.append(np.std(x, ddof=1) / np.sqrt(len(x)))
        if len(x) % 2:
            x = x[:-1]
        x = 0.5 * (x[0::2] + x[1::2])
    return float(max(sems)) if sems else 0.0


def auroc(scores, labels):
    pos = scores[labels == 1]; neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    nc = sum(1 for p in pos for q in neg if p > q)
    nt = sum(1 for p in pos for q in neg if p == q)
    return (nc + 0.5 * nt) / (len(pos) * len(neg))


def f156_auroc(var_means: dict):
    """AUROC for lower-F156 = deleterious (5 del vs 3 neu, exclude WT/C212A)."""
    vs = [v for v in var_means if v in TRUTH_DEL or v in NEU]
    scores = np.array([-var_means[v] for v in vs])   # sign: lower F156 -> higher del score
    labels = np.array([1 if v in TRUTH_DEL else 0 for v in vs])
    return auroc(scores, labels), len([v for v in vs if v in TRUTH_DEL]), len([v for v in vs if v in NEU])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="+", default=ALL_VARIANTS)
    args = ap.parse_args()

    rows = []
    for v in args.variants:
        for src_name, base in SOURCES:
            traj = load_traj(base, v)
            if traj is None:
                continue
            ns = traj.n_frames * NS_PER_FRAME
            s = f156_hbond_series(traj)
            half = len(s) // 2
            tau = integrated_autocorr_time(s)
            rows.append({"variant": v, "class": CLASS.get(v, "?"), "source": src_name,
                         "ns": round(ns, 1), "n_frames": len(s),
                         "f156_mean": round(float(np.nanmean(s)), 4),
                         "f156_1sthalf": round(float(np.nanmean(s[:half])), 4),
                         "f156_2ndhalf": round(float(np.nanmean(s[half:])), 4),
                         "frac_zero": round(float(np.mean(s == 0)), 3),
                         "blocking_sem": round(blocking_sem(s), 4),
                         "tau_int_ns": round(tau * NS_PER_FRAME, 2),
                         "n_eff": round(len(s) / (2 * tau), 1), "_series": s})
            print(f"  {v:6s} {src_name:13s} {ns:5.0f}ns  mean={rows[-1]['f156_mean']:.3f}  "
                  f"1st={rows[-1]['f156_1sthalf']:.3f} 2nd={rows[-1]['f156_2ndhalf']:.3f}  "
                  f"frac0={rows[-1]['frac_zero']:.2f} n_eff={rows[-1]['n_eff']:.0f}"
                  + ("  [PARTIAL <180ns]" if ns < MIN_NS else ""), flush=True)

    df = pd.DataFrame([{k: r[k] for k in r if k != "_series"} for r in rows])
    df.to_csv(OUT / "diagnostic_f156_convergence.csv", index=False)

    print("\n" + "=" * 80)
    print(f"PER-VARIANT SPREAD (valid >={MIN_NS:.0f}ns reps only; partials shown but excluded)")
    print("=" * 80)
    var_means, var_means_droplow = {}, {}
    for v in args.variants:
        sub = df[df["variant"] == v].sort_values("f156_mean")
        valid = sub[sub["ns"] >= MIN_NS]
        if valid.empty:
            continue
        m = valid["f156_mean"].values
        var_means[v] = float(m.mean())
        var_means_droplow[v] = float(m[1:].mean()) if len(m) > 1 else float(m.mean())
        med_oth = np.median(m[1:]) if len(m) > 1 else np.nan
        outlier = (len(m) > 1 and med_oth > 0 and m[0] < 0.6 * med_oth)
        print(f"\n{v} ({CLASS.get(v,'?')}): valid_reps={len(m)}  mean={var_means[v]:.3f}  "
              f"drop-lowest={var_means_droplow[v]:.3f}" + ("  <-- OUTLIER-DEP" if outlier else ""))
        for _, r in sub.iterrows():
            tag = "  [PARTIAL]" if r["ns"] < MIN_NS else ""
            print(f"    {r['source']:13s} {r['ns']:5.0f}ns  F156={r['f156_mean']:.3f} "
                  f"frac0={r['frac_zero']:.2f} drift={r['f156_2ndhalf']-r['f156_1sthalf']:+.3f} "
                  f"n_eff={r['n_eff']:.0f}{tag}")

    # ---- Recompute F156 AUROC under this clean sampling ----
    print("\n" + "=" * 80)
    print("F156 AUROC UNDER CLEAN 1ns SAMPLING  (lower F156 = deleterious; 5 del vs 3 neu)")
    print("=" * 80)
    if all(v in var_means for v in TRUTH_DEL | NEU):
        a_full, nd, nn = f156_auroc(var_means)
        a_drop, _, _ = f156_auroc(var_means_droplow)
        print(f"  replicate-averaged AUROC = {a_full:.3f}   ({nd} del vs {nn} neu)")
        print(f"  drop-lowest-rep  AUROC   = {a_drop:.3f}   (sensitivity to the lower replicates)")
        print("\n  per-variant means (sorted):")
        for v in sorted(var_means, key=var_means.get):
            if v in TRUTH_DEL or v in NEU:
                print(f"    {v:6s} {CLASS[v]:3s}  {var_means[v]:.3f}")
    else:
        missing = [v for v in TRUTH_DEL | NEU if v not in var_means]
        print(f"  (need all del+neu variants; missing valid data for: {missing})")

    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        vs = [v for v in args.variants if not df[df["variant"] == v].empty]
        ncol = 3; nrow = int(np.ceil(len(vs) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(5 * ncol, 3.2 * nrow), squeeze=False)
        for i, v in enumerate(vs):
            ax = axes[i // ncol][i % ncol]
            for r in rows:
                if r["variant"] != v:
                    continue
                s = r["_series"]; t = np.arange(len(s)) * NS_PER_FRAME
                run = np.cumsum(s) / (np.arange(len(s)) + 1)
                ls = "--" if r["ns"] < MIN_NS else "-"
                ax.plot(t, run, lw=1.3, ls=ls, label=f"{r['source']} ({r['f156_mean']:.2f})")
            ax.set_title(f"{v} ({CLASS.get(v,'?')})"); ax.set_xlabel("ns")
            ax.set_ylabel("cum. mean F156-RNA Hb"); ax.grid(alpha=0.3)
            ax.legend(fontsize=6, frameon=False)
        for j in range(len(vs), nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        fig.tight_layout(); fig.savefig(PLOTS / "diag_f156_running_means.png", dpi=150)
        print(f"\nsaved plot: {PLOTS / 'diag_f156_running_means.png'}")
    except Exception as e:
        print(f"\n(plot skipped: {e})")
    print(f"wrote {OUT / 'diagnostic_f156_convergence.csv'}")


if __name__ == "__main__":
    main()
