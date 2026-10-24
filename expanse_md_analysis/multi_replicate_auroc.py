#!/usr/bin/env python3
"""
Multi-replicate AUROC analysis across all available MD runs — FULL feature set.

Sources:
  1. Expanse original (400ns runs, capped at first 200ns)
  2. Anvil rep1 (all 10 variants, ~185-200ns)
  3. Anvil rep2 (most variants ~200ns, some incomplete excluded)
  4. Expanse rep2 (only completed variants included)

Features computed per source (matching the full catalog from feature_aurocs.csv):
  MD bulk: rmsd_ca, rmsd_zf1, rmsd_zf2, rmsd_rnaP, contacts (all last-20% mean)
  MD H-bonds: per-residue H-bonds to RNA for 14 key residues
  MD energy: screened Coulomb, salt bridges, interface atom count, binding-residue RMSF
  AF3 static: pair_ipTM, PAE, pLDDT (whole/ZF1/ZF2/min), RMSD_whole (from CSV, same for all sources)
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

# ZFP36L2 residue number → PDB local resSeq (offset = -150)
ZFP_TO_PDB = {154: 4, 156: 6, 160: 10, 170: 20, 171: 21, 184: 34, 185: 35,
              187: 37, 200: 50, 204: 54, 206: 56, 208: 58, 210: 60, 212: 62}
RES_NAMES = {154: "Y154", 156: "F156", 160: "R160", 170: "H170", 171: "R171",
             184: "R184", 185: "S185", 187: "T187", 200: "R200", 204: "G204",
             206: "C206", 208: "C208", 210: "D210", 212: "C212"}
ZF1_PDB = list(range(2, 30))   # 152-179
ZF2_PDB = list(range(40, 68))  # 190-217

SOURCES = [
    ("expanse_orig", ROOT / "expanse_md_results", 500),
    ("anvil_rep1", ROOT / "anvil_md_results" / "rep1", 450),
    ("anvil_rep2", ROOT / "anvil_md_results" / "rep2", 450),
    ("expanse_rep2", ROOT / "expanse_md_results_rep2", 450),
]
MAX_NS = 200.0


def rna_sel():
    return " or ".join(f"resname {r}" for r in RNA_RES)


def load_traj(base_dir, variant, max_ns=MAX_NS):
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
        print(f"  load error {variant}: {e}")
        return None
    nmax = int(max_ns * 10) + 1
    if traj.n_frames > nmax:
        traj = traj[:nmax]
    return traj


def compute_all_metrics(traj, variant):
    """Compute the full feature set for one trajectory."""
    top = traj.topology
    n = traj.n_frames
    last20 = max(1, int(0.2 * n))
    row = {"variant": variant, "class": CLASS[variant], "n_frames": n, "ns_approx": n * 0.1}

    # --- Bulk RMSD metrics ---
    ca = top.select("protein and name CA")
    if len(ca) > 0:
        rmsd = md.rmsd(traj, traj, 0, atom_indices=ca) * 10
        row["rmsd_ca_last20pct"] = float(np.mean(rmsd[-last20:]))

    # ZF1 / ZF2 Cα RMSD
    for label, pdb_range in [("zf1", ZF1_PDB), ("zf2", ZF2_PDB)]:
        sel = top.select("protein and name CA and (" +
                         " or ".join(f"resSeq {p}" for p in pdb_range) + ")")
        if len(sel) > 0:
            r = md.rmsd(traj, traj, 0, atom_indices=sel) * 10
            row[f"rmsd_{label}_last20pct"] = float(np.mean(r[-last20:]))

    # RNA-P RMSD
    rna_p = top.select(f"({rna_sel()}) and name P")
    if len(rna_p) > 0:
        r = md.rmsd(traj, traj, 0, atom_indices=rna_p) * 10
        row["rmsd_rnaP_last20pct"] = float(np.mean(r[-last20:]))

    # --- Protein-RNA contacts ---
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

    # --- Per-residue H-bonds to RNA (baker_hubbard on subsampled frames) ---
    rna_atoms = {a.index for a in top.atoms if a.residue.name.upper() in RNA_RES}
    hb_stride = max(1, n // 50)
    sub = traj[::hb_stride]
    # Build residue→atom map for key residues
    res_atom_map = {}  # pdb_resSeq -> set of atom indices
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

    for zfp, pdb in ZFP_TO_PDB.items():
        row[f"hbond_RNA_{RES_NAMES[zfp]}"] = float(np.mean(hb_counts[zfp]))

    # --- Screened Coulomb ---
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

    # --- Salt bridges (cation-phosphate <4Å) ---
    rna_phos = top.select(f"({rna_sel()}) and (name OP1 or name OP2 or name O1P or name O2P)")
    cation = top.select("protein and (name NZ or name NH1 or name NH2 or name NE)")
    if len(rna_phos) > 0 and len(cation) > 0:
        pairs = np.array(np.meshgrid(cation, rna_phos)).T.reshape(-1, 2)
        d = md.compute_distances(traj[::max(1, n // 50)], pairs, periodic=True)
        row["salt_bridges_prot_rna"] = float((d < 0.4).sum(axis=1).mean())

    # --- Interface atom count (<5Å) ---
    if len(prot_heavy) > 0 and len(rna_heavy) > 0:
        counts = []
        for fi in range(0, n, max(1, n // 50)):
            nb = md.compute_neighbors(traj[fi], 0.5, rna_heavy, haystack_indices=prot_heavy)
            counts.append(len(nb[0]))
        row["interface_atom_count_5A"] = float(np.mean(counts))

    # --- Binding-residue RMSF ---
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


# Feature → (sign, origin). sign: +1 = higher value → more deleterious
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
# Add per-residue H-bonds (fewer H-bonds → deleterious → sign = -1)
for zfp in ZFP_TO_PDB:
    FEATURE_SPEC[f"hbond_RNA_{RES_NAMES[zfp]}"] = (-1, "MD-Hbond")

# AF3 static features (loaded from CSV, same for all sources)
AF3_FEATURES = {
    "pair_ipTM_protein_rna": (-1, "AF3"),
    "pair_pae_min_cross_mean": (+1, "AF3"),
    "pLDDT": (-1, "AF3"),
    "pLDDT_ZF1": (-1, "AF3"),
    "pLDDT_ZF2": (-1, "AF3"),
    "pLDDT_min": (-1, "AF3"),
    "RMSD_whole": (+1, "AF3"),
}


def compute_aurocs(df, feature_spec):
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
    # Load AF3 static metrics (same for every source)
    af3_csv = ROOT / "results" / "variant_avg_all.csv"
    af3 = None
    if af3_csv.exists():
        af3 = pd.read_csv(af3_csv).rename(columns={"Protein": "variant"})
        af3.loc[af3["variant"] == "C208A", "variant"] = "C212A"

    all_results, all_aurocs = [], []

    for src_name, base_dir, min_mb in SOURCES:
        print(f"\n{'='*60}")
        print(f"SOURCE: {src_name}")
        print(f"{'='*60}")

        rows = []
        for v in VARIANTS:
            dcd_path = None
            for vname in [v, "C208A" if v == "C212A" else v]:
                p = base_dir / vname / "production.dcd"
                if p.exists():
                    dcd_path = p
                    break
            if dcd_path is None:
                print(f"  {v}: no DCD, skip")
                continue
            sz_mb = dcd_path.stat().st_size / 1e6
            if sz_mb < min_mb:
                print(f"  {v}: incomplete ({sz_mb:.0f}MB), skip")
                continue

            print(f"  {v}: {sz_mb:.0f}MB...", end=" ", flush=True)
            traj = load_traj(base_dir, v)
            if traj is None:
                print("FAIL")
                continue
            metrics = compute_all_metrics(traj, v)
            print(f"{metrics['ns_approx']:.0f}ns OK")
            rows.append(metrics)

        if not rows:
            continue

        df = pd.DataFrame(rows)

        # Merge AF3 static features
        if af3 is not None:
            df = df.merge(af3, on="variant", how="left")

        df.insert(0, "source", src_name)
        all_results.append(df)

        # Compute AUROCs for MD features
        all_spec = {**FEATURE_SPEC}
        if af3 is not None:
            all_spec.update(AF3_FEATURES)

        adf = compute_aurocs(df, all_spec)
        if not adf.empty:
            adf.insert(0, "source", src_name)
            adf["n_variants"] = len(df[df["variant"] != "WT"])
            all_aurocs.append(adf)
            print(f"\n  Top 10 AUROCs ({src_name}, {adf['n_variants'].iloc[0]} mut variants):")
            print(adf.head(10)[["feature", "AUROC", "origin"]].to_string(index=False))

    # --- Combine + summarize ---
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(OUT / "multi_replicate_metrics.csv", index=False)
        print(f"\nWrote {OUT / 'multi_replicate_metrics.csv'}")

    if all_aurocs:
        auroc_all = pd.concat(all_aurocs, ignore_index=True)
        auroc_all.to_csv(OUT / "multi_replicate_aurocs.csv", index=False)

        summary = auroc_all.groupby(["feature", "origin"])["AUROC"].agg(
            ["mean", "std", "count"]).reset_index()
        summary.columns = ["feature", "origin", "AUROC_mean", "AUROC_std", "n_sources"]
        summary = summary.sort_values("AUROC_mean", ascending=False)
        summary.to_csv(OUT / "multi_replicate_auroc_summary.csv", index=False)

        print(f"\n{'='*60}")
        print("CROSS-REPLICATE AUROC SUMMARY (mean ± std)")
        print(f"{'='*60}")
        for _, r in summary.iterrows():
            std_str = f"±{r['AUROC_std']:.3f}" if pd.notna(r['AUROC_std']) else "±NA"
            print(f"  {r['feature']:30s} {r['AUROC_mean']:.3f} {std_str}  "
                  f"({int(r['n_sources'])} srcs)  [{r['origin']}]")


if __name__ == "__main__":
    main()
