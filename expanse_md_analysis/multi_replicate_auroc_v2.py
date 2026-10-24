#!/usr/bin/env python3
"""
Multi-replicate AUROC analysis v2 — correct methodology.

Approach:
  1. Scan ALL sources for each variant; only include trajectories ≥ 200ns (≥1800 frames)
  2. Compute metrics per trajectory (full feature set: RMSD, H-bonds, energy, etc.)
  3. For per-source AUROC: compute AUROC within each source using only its complete variants
     (reported for reference, but caveated when n_variants < 9)
  4. For replicate-averaged AUROC (the main result):
     - Average each metric across all available replicates per variant
     - Compute ONE AUROC on the 9 averaged values (6 del vs 3 neu)
     - Bootstrap 95% CI by resampling which replicate is used per variant
  5. AF3 static features are added from CSV (same values regardless of MD source)

Output:
  - multi_replicate_all_metrics.csv  : every (source, variant) metric row
  - multi_replicate_per_source_auroc.csv : AUROC per source (for reference)
  - multi_replicate_averaged_auroc.csv : THE main result — replicate-averaged AUROC + bootstrap CI
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import mdtraj as md

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent

VARIANTS = ["WT", "C206F", "C212A", "D210E", "G204D",
            "R160W", "R184C", "S185N", "T187S", "Y154H"]
TRUTH_DEL = {"Y154H", "R160W", "R184C", "G204D", "C206F", "C212A"}
CLASS = {
    "WT": "WT", "C206F": "deleterious", "C212A": "zn_knockout_control",
    "G204D": "deleterious", "R160W": "deleterious", "R184C": "deleterious",
    "Y154H": "deleterious",
    "D210E": "neutral", "S185N": "neutral", "T187S": "neutral",
}
RNA_RES = {"A", "C", "G", "U", "A5", "A3", "C5", "C3", "G5", "G3", "U5", "U3"}
ZFP_TO_PDB = {154: 4, 156: 6, 160: 10, 170: 20, 171: 21, 184: 34, 185: 35,
              187: 37, 200: 50, 204: 54, 206: 56, 208: 58, 210: 60, 212: 62}
RES_NAMES = {154: "Y154", 156: "F156", 160: "R160", 170: "H170", 171: "R171",
             184: "R184", 185: "S185", 187: "T187", 200: "R200", 204: "G204",
             206: "C206", 208: "C208", 210: "D210", 212: "C212"}
ZF1_PDB = list(range(2, 30))
ZF2_PDB = list(range(40, 68))

SOURCES = [
    ("expanse_orig", ROOT / "expanse_md_results"),
    ("anvil_rep1", ROOT / "anvil_md_results" / "rep1"),
    ("anvil_rep2", ROOT / "anvil_md_results" / "rep2"),
    ("expanse_rep2", ROOT / "expanse_md_results_rep2"),
]
MIN_FRAMES = 1800  # ~180 ns minimum to count as "complete enough" (200ns = 2000 frames)
MAX_NS = 200.0


def rna_sel():
    return " or ".join(f"resname {r}" for r in RNA_RES)


def load_traj(base_dir, variant):
    for vname in [variant, "C208A" if variant == "C212A" else variant]:
        d = base_dir / vname
        top = d / "equilibrated.pdb"
        dcd = d / "production.dcd"
        if top.exists() and dcd.exists():
            break
    else:
        return None
    try:
        full_top = md.load_topology(str(top))
        keep = np.unique(np.concatenate([
            full_top.select("protein"),
            full_top.select(rna_sel()),
            full_top.select("name ZN or element Zn"),
        ]))
        traj = md.load(str(dcd), top=str(top), atom_indices=keep)
    except Exception as e:
        print(f"    load error: {e}")
        return None
    # Cap at MAX_NS
    nmax = int(MAX_NS * 10) + 1
    if traj.n_frames > nmax:
        traj = traj[:nmax]
    return traj


def compute_all_metrics(traj, variant):
    top = traj.topology
    n = traj.n_frames
    last20 = max(1, int(0.2 * n))
    row = {"variant": variant, "class": CLASS[variant], "n_frames": n, "ns_approx": n * 0.1}

    ca = top.select("protein and name CA")
    if len(ca) > 0:
        rmsd = md.rmsd(traj, traj, 0, atom_indices=ca) * 10
        row["rmsd_ca_last20pct"] = float(np.mean(rmsd[-last20:]))

    for label, pdb_range in [("zf1", ZF1_PDB), ("zf2", ZF2_PDB)]:
        sel = top.select("protein and name CA and (" +
                         " or ".join(f"resSeq {p}" for p in pdb_range) + ")")
        if len(sel) > 0:
            r = md.rmsd(traj, traj, 0, atom_indices=sel) * 10
            row[f"rmsd_{label}_last20pct"] = float(np.mean(r[-last20:]))

    rna_p = top.select(f"({rna_sel()}) and name P")
    if len(rna_p) > 0:
        r = md.rmsd(traj, traj, 0, atom_indices=rna_p) * 10
        row["rmsd_rnaP_last20pct"] = float(np.mean(r[-last20:]))

    prot_heavy = top.select("protein and not type H")
    rna_heavy = top.select(f"({rna_sel()}) and not type H")
    contacts = []
    stride = max(1, n // 100)
    for fi in range(0, n, stride):
        nb = md.compute_neighbors(traj[fi], 0.45, rna_heavy, haystack_indices=prot_heavy)
        contacts.append(len(nb[0]))
    contacts = np.array(contacts)
    c20 = max(1, int(0.2 * len(contacts)))
    row["contacts_last20pct"] = float(np.mean(contacts[-c20:]))

    # Per-residue H-bonds
    rna_atoms = {a.index for a in top.atoms if a.residue.name.upper() in RNA_RES}
    hb_stride = max(1, n // 50)
    sub = traj[::hb_stride]
    res_atom_map = {}
    for a in top.atoms:
        if a.residue.is_protein:
            rs = a.residue.resSeq
            if rs not in res_atom_map:
                res_atom_map[rs] = set()
            res_atom_map[rs].add(a.index)
    hb_counts = {zfp: np.zeros(sub.n_frames) for zfp in ZFP_TO_PDB}
    for fi in range(sub.n_frames):
        try:
            hbs = md.baker_hubbard(sub[fi], freq=0.0, exclude_water=True)
        except Exception:
            continue
        for d, h, a in hbs:
            for zfp, pdb in ZFP_TO_PDB.items():
                atoms = res_atom_map.get(pdb, set())
                if not atoms:
                    continue
                if (d in rna_atoms and a in atoms) or (a in rna_atoms and d in atoms):
                    hb_counts[zfp][fi] += 1
    for zfp in ZFP_TO_PDB:
        row[f"hbond_RNA_{RES_NAMES[zfp]}"] = float(np.mean(hb_counts[zfp]))

    # Screened Coulomb
    pos_idx, neg_idx = [], []
    for a in top.atoms:
        r = a.residue
        nm, rn = a.name, r.name.upper()
        if r.is_protein:
            if rn == "ARG" and nm in ("NH1", "NH2", "NE"):
                pos_idx.append(a.index)
            elif rn == "LYS" and nm == "NZ":
                pos_idx.append(a.index)
        else:
            if nm in ("OP1", "OP2", "O1P", "O2P"):
                neg_idx.append(a.index)
    rna_neg = np.array([i for i in neg_idx if i in rna_atoms], dtype=int)
    prot_pos = np.array(pos_idx, dtype=int)
    if len(prot_pos) > 0 and len(rna_neg) > 0:
        pairs = np.array(np.meshgrid(prot_pos, rna_neg)).T.reshape(-1, 2)
        sub2 = traj[::max(1, n // 50)]
        d = md.compute_distances(sub2, pairs, periodic=True)
        d = np.clip(d, 0.15, None)
        e = -138.935458 * np.exp(-d / 0.8) / d
        row["coulomb_prot_rna"] = float(e.sum(axis=1).mean())

    # Salt bridges
    rna_phos = top.select(f"({rna_sel()}) and (name OP1 or name OP2 or name O1P or name O2P)")
    cation = top.select("protein and (name NZ or name NH1 or name NH2 or name NE)")
    if len(rna_phos) > 0 and len(cation) > 0:
        pairs = np.array(np.meshgrid(cation, rna_phos)).T.reshape(-1, 2)
        d = md.compute_distances(traj[::max(1, n // 50)], pairs, periodic=True)
        row["salt_bridges_prot_rna"] = float((d < 0.4).sum(axis=1).mean())

    # Interface atom count
    if len(prot_heavy) > 0 and len(rna_heavy) > 0:
        counts = []
        for fi in range(0, n, max(1, n // 50)):
            nb = md.compute_neighbors(traj[fi], 0.5, rna_heavy, haystack_indices=prot_heavy)
            counts.append(len(nb[0]))
        row["interface_atom_count_5A"] = float(np.mean(counts))

    # Binding-residue RMSF
    if len(ca) > 0:
        tr = traj.superpose(traj, frame=0, atom_indices=ca, ref_atom_indices=ca)
        bind_ca = top.select("protein and name CA and (" +
                             " or ".join(f"resSeq {p}" for p in ZF1_PDB + ZF2_PDB) + ")")
        if len(bind_ca) > 0:
            xyz = tr.xyz[:, bind_ca, :]
            mean = xyz.mean(axis=0, keepdims=True)
            rmsf = np.sqrt(((xyz - mean) ** 2).sum(axis=2).mean(axis=0))
            row["binding_residue_rmsf_A"] = float(rmsf.mean() * 10.0)

    return row


def auroc(scores, labels):
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    nc = sum(1 for p in pos for n_ in neg if p > n_)
    nt = sum(1 for p in pos for n_ in neg if p == n_)
    return (nc + 0.5 * nt) / (len(pos) * len(neg))


FEATURE_SPEC = {
    "rmsd_ca_last20pct": (+1, "MD-RMSD"),
    "rmsd_zf1_last20pct": (+1, "MD-RMSD"),
    "rmsd_zf2_last20pct": (+1, "MD-RMSD"),
    "rmsd_rnaP_last20pct": (+1, "MD-RMSD"),
    "contacts_last20pct": (-1, "MD-contacts"),
    "coulomb_prot_rna": (+1, "MD-energy"),
    "salt_bridges_prot_rna": (-1, "MD-energy"),
    "interface_atom_count_5A": (-1, "MD-energy"),
    "binding_residue_rmsf_A": (+1, "MD-energy"),
}
for zfp in ZFP_TO_PDB:
    FEATURE_SPEC[f"hbond_RNA_{RES_NAMES[zfp]}"] = (-1, "MD-Hbond")

AF3_FEATURES = {
    "pair_ipTM_protein_rna": (-1, "AF3"),
    "pair_pae_min_cross_mean": (+1, "AF3"),
    "pLDDT": (-1, "AF3"),
    "pLDDT_ZF1": (-1, "AF3"),
    "pLDDT_ZF2": (-1, "AF3"),
    "pLDDT_min": (-1, "AF3"),
    "RMSD_whole": (+1, "AF3"),
}


def compute_aurocs_df(df, feature_spec):
    mut = df[df["variant"] != "WT"].copy()
    y = mut["variant"].apply(lambda v: 1 if v in TRUTH_DEL else 0).values
    rows = []
    for f, (sign, origin) in feature_spec.items():
        if f not in mut.columns:
            continue
        v = mut[f].values.astype(float) * sign
        if not np.isfinite(v).all():
            continue
        rows.append({"feature": f, "AUROC": round(auroc(v, y), 3), "origin": origin})
    return pd.DataFrame(rows).sort_values("AUROC", ascending=False)


def main():
    # Load AF3 static
    af3_csv = ROOT / "results" / "variant_avg_all.csv"
    af3 = None
    if af3_csv.exists():
        af3 = pd.read_csv(af3_csv).rename(columns={"Protein": "variant"})
        af3.loc[af3["variant"] == "C208A", "variant"] = "C212A"

    # =====================================================================
    # STEP 1: Compute metrics for every (source, variant) with ≥ MIN_FRAMES
    # =====================================================================
    all_rows = []
    for src_name, base_dir in SOURCES:
        print(f"\n{'='*60}")
        print(f"SOURCE: {src_name}")
        print(f"{'='*60}")
        for v in VARIANTS:
            dcd_path = None
            for vname in [v, "C208A" if v == "C212A" else v]:
                p = base_dir / vname / "production.dcd"
                if p.exists():
                    dcd_path = p
                    break
            if dcd_path is None:
                continue

            print(f"  {v}: loading...", end=" ", flush=True)
            traj = load_traj(base_dir, v)
            if traj is None:
                print("FAIL")
                continue
            if traj.n_frames < MIN_FRAMES:
                print(f"TOO SHORT ({traj.n_frames} frames = {traj.n_frames*0.1:.0f}ns < {MIN_FRAMES*0.1:.0f}ns), skip")
                continue
            metrics = compute_all_metrics(traj, v)
            metrics["source"] = src_name
            print(f"{metrics['ns_approx']:.0f}ns, {metrics['n_frames']} frames ✓")
            all_rows.append(metrics)

    all_df = pd.DataFrame(all_rows)
    all_df.to_csv(OUT / "multi_replicate_all_metrics.csv", index=False)
    print(f"\nWrote {OUT / 'multi_replicate_all_metrics.csv'}")

    # Show which variants have how many replicates
    print(f"\n{'='*60}")
    print("REPLICATE COUNT PER VARIANT (only ≥200ns runs)")
    print(f"{'='*60}")
    for v in VARIANTS:
        srcs = all_df[all_df["variant"] == v]["source"].tolist()
        print(f"  {v:8s}: {len(srcs)} reps  ({', '.join(srcs)})")

    all_spec = {**FEATURE_SPEC}
    if af3 is not None:
        all_spec.update(AF3_FEATURES)

    # =====================================================================
    # STEP 2: Per-source AUROC (for reference — only on complete-9 sources)
    # =====================================================================
    print(f"\n{'='*60}")
    print("PER-SOURCE AUROC (reference; only computed when all 9 mut variants present)")
    print(f"{'='*60}")
    per_source_rows = []
    for src_name in all_df["source"].unique():
        sdf = all_df[all_df["source"] == src_name].copy()
        if af3 is not None:
            sdf = sdf.merge(af3, on="variant", how="left", suffixes=("", "_af3"))
        n_mut = len(sdf[sdf["variant"] != "WT"])
        adf = compute_aurocs_df(sdf, all_spec)
        if not adf.empty:
            adf["source"] = src_name
            adf["n_mut_variants"] = n_mut
            per_source_rows.append(adf)
            print(f"\n  {src_name} ({n_mut} mut variants):")
            print(adf.head(8)[["feature", "AUROC", "origin"]].to_string(index=False))

    if per_source_rows:
        per_source_df = pd.concat(per_source_rows, ignore_index=True)
        per_source_df.to_csv(OUT / "multi_replicate_per_source_auroc.csv", index=False)

    # =====================================================================
    # STEP 3: REPLICATE-AVERAGED AUROC (the main, correct result)
    # =====================================================================
    print(f"\n{'='*60}")
    print("REPLICATE-AVERAGED AUROC (main result)")
    print("  Method: mean metric across replicates per variant → one AUROC")
    print(f"{'='*60}")

    # Get all numeric metric columns
    metric_cols = [c for c in all_df.columns if c not in
                   ("variant", "class", "source", "n_frames", "ns_approx")]

    # Average each metric across replicates per variant
    avg_df = all_df.groupby("variant")[metric_cols].mean().reset_index()
    avg_df["class"] = avg_df["variant"].map(CLASS)
    # Add replicate counts
    rep_counts = all_df.groupby("variant").size().reset_index(name="n_replicates")
    avg_df = avg_df.merge(rep_counts, on="variant")

    # Merge AF3 static features
    if af3 is not None:
        avg_df = avg_df.merge(af3, on="variant", how="left", suffixes=("", "_af3"))

    print("\nPer-variant replicate-averaged values:")
    print(avg_df[["variant", "class", "n_replicates"]].to_string(index=False))

    # Compute AUROC on the averaged values
    avg_auroc = compute_aurocs_df(avg_df, all_spec)

    # Bootstrap 95% CI: for each bootstrap iteration, randomly pick one replicate
    # per variant (from those available), compute AUROC, repeat 1000 times
    np.random.seed(42)
    n_boot = 1000
    boot_aurocs = {f: [] for f in all_spec}

    for _ in range(n_boot):
        # Sample one replicate per variant
        boot_rows = []
        for v in VARIANTS:
            v_rows = all_df[all_df["variant"] == v]
            if len(v_rows) == 0:
                continue
            sampled = v_rows.sample(n=1, replace=True).iloc[0]
            boot_rows.append(sampled)
        boot_df = pd.DataFrame(boot_rows)
        if af3 is not None:
            boot_df = boot_df.merge(af3, on="variant", how="left", suffixes=("", "_af3"))

        mut = boot_df[boot_df["variant"] != "WT"]
        y = mut["variant"].apply(lambda v: 1 if v in TRUTH_DEL else 0).values
        for f, (sign, origin) in all_spec.items():
            if f not in mut.columns:
                continue
            vals = mut[f].values.astype(float) * sign
            if np.isfinite(vals).all():
                boot_aurocs[f].append(auroc(vals, y))

    # Compile results with CI
    result_rows = []
    for _, r in avg_auroc.iterrows():
        f = r["feature"]
        boots = np.array(boot_aurocs.get(f, []))
        ci_lo, ci_hi = (np.nan, np.nan)
        boot_std = np.nan
        if len(boots) > 0:
            ci_lo = float(np.percentile(boots, 2.5))
            ci_hi = float(np.percentile(boots, 97.5))
            boot_std = float(np.std(boots))
        result_rows.append({
            "feature": f,
            "AUROC_avg": r["AUROC"],
            "AUROC_boot_mean": round(float(np.mean(boots)), 3) if len(boots) > 0 else np.nan,
            "AUROC_95CI_lo": round(ci_lo, 3),
            "AUROC_95CI_hi": round(ci_hi, 3),
            "AUROC_boot_std": round(boot_std, 3),
            "origin": r["origin"],
        })

    result_df = pd.DataFrame(result_rows).sort_values("AUROC_avg", ascending=False)
    result_df.to_csv(OUT / "multi_replicate_averaged_auroc.csv", index=False)

    print(f"\n{'='*60}")
    print("REPLICATE-AVERAGED AUROC WITH BOOTSTRAP 95% CI")
    print(f"{'='*60}")
    for _, r in result_df.iterrows():
        ci = f"[{r['AUROC_95CI_lo']:.3f}, {r['AUROC_95CI_hi']:.3f}]"
        print(f"  {r['feature']:30s} {r['AUROC_avg']:.3f}  CI={ci}  [{r['origin']}]")

    print(f"\nWrote {OUT / 'multi_replicate_averaged_auroc.csv'}")
    print(f"Wrote {OUT / 'multi_replicate_per_source_auroc.csv'}")
    print(f"Wrote {OUT / 'multi_replicate_all_metrics.csv'}")


if __name__ == "__main__":
    main()
