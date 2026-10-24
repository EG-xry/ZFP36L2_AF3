# Feature-by-feature discrimination of deleterious vs neutral ZFP36L2 nsSNPs

**Figure:** `plots/fig_feature_aurocs.png`
**Source data:** `feature_aurocs.csv`, `af3_vs_md_combined.csv`,
`common139_summary.csv`, `hbonds_per_residue_common139.csv`
**Code:** `make_auroc_figure.py`, `plot_drift_and_hbonds.py`,
`run_full_analysis.py`

---

## 1. What is AUROC and how was it computed

For each *single number* per variant (e.g., "Cα-RMSD at end of run", or
"AF3 pair-ipTM", or "average H-bonds at R160 to RNA"), we asked:

> *If I rank all 9 mutants by this number, how often does a deleterious
> variant rank higher (or lower) than a neutral one?*

This is the **Area Under the Receiver Operating Characteristic curve
(AUROC)**.  Concretely, for our 6 deleterious + 3 neutral variants
there are 6 × 3 = 18 ordered (del, neu) pairs.  AUROC = fraction of
those pairs ordered the way the feature predicts:

- **AUROC = 0.5** — the feature is random (worthless).
- **AUROC = 1.0** — the feature gets every (del, neu) pair right.
- **AUROC = 0.0** — the feature is perfectly *anti*-correlated (a clean
  signal in the wrong direction; can be flipped to ≈1.0 by negating).

We chose the sign of every feature so that "higher score = more
deleterious".  For pair-ipTM (where higher = better quality) we negate
it; for whole-protein RMSD (where higher = more drift) we leave the
sign alone, etc.  Exact signs are in
`run_full_analysis.py:step_4_aurocs()`.

### Pairing of variants used here

| variant | truth     | reason |
|---|---|---|
| Y154H | deleterious | dbSNP rs375096815 + experimental ZF1 binding loss |
| R160W | deleterious | dbSNP rs1183688047 + experimental ZF1 binding loss |
| R184C | deleterious | dbSNP rs1214015428 + experimental ZF2 binding loss |
| G204D | deleterious | dbSNP rs1215671792 + experimental ZF2 binding loss |
| C206F | deleterious | dbSNP rs920398592  + experimental ZF2 binding loss |
| C212A | deleterious (Zn-knockout control) | engineered Zn-coord Cys → Ala |
| S185N | neutral | predicted neutral; experimentally normal binding |
| T187S | neutral | predicted neutral; experimentally normal binding |
| D210E | neutral | predicted neutral; experimentally normal binding (actual mutation in this run is D219E — see VARIANT_LABELS_AUDIT.md) |

WT is excluded because it has no class label to predict.

---

## 2. The features we tested

The full AUROC table is in `feature_aurocs.csv`; here it is grouped by
origin.

### AF3 (single static structure, no MD)

| feature | what it measures | AUROC |
|---|---|---:|
| `pair_ipTM (protein↔RNA)` | AF3's confidence that the protein–RNA interface is correctly predicted (= 0.8·ipTM + 0.2·pTM of the cross-chain pair) | **0.778** |
| `PAE_min (cross-chain mean)` | mean Predicted Aligned Error across the protein↔RNA interface | 0.722 |
| `pLDDT (whole)` | mean per-residue confidence across the whole protein | 0.389 |
| `pLDDT_ZF1` | mean pLDDT restricted to zinc finger 1 | 0.500 |
| `pLDDT_ZF2` | mean pLDDT restricted to zinc finger 2 | 0.389 |
| `pLDDT_min` | the worst-pLDDT residue in the complex | 0.556 |
| `RMSD_whole` | RMSD between AF3 model and 1RGO experimental structure | 0.556 |

### MD bulk (whole-system summaries from the 139 ns common window)

| feature | what it measures | AUROC |
|---|---|---:|
| `rmsd_rnaP_last20pct_mean` | RMSD of all 8 RNA phosphate atoms vs the equilibrated frame, averaged over the last 20 % of the window — i.e., **how far the RNA has slid out of the binding pocket** | 0.556 |
| `contacts_last20pct_mean` | number of protein-heavy / RNA-heavy atom pairs within 4 Å, average | 0.333 |
| `rmsd_zf2_last20pct_mean` | Cα RMSD of zinc finger 2 only (ZFP36L2 res 190–217) | 0.333 |
| `rmsd_zf1_last20pct_mean` | Cα RMSD of zinc finger 1 only (ZFP36L2 res 152–179) | 0.333 |
| `rmsd_ca_last20pct_mean` | Cα RMSD of the entire 70-aa ZFD | 0.278 |

These bulk metrics are uninformative on this dataset — they mix many
small motions of irrelevant residues with a few large motions of
binding residues, so the signal-to-noise is poor.

### MD per-residue H-bonds (this is where the real signal lives)

For each canonical ARE-binding residue of ZFP36L2, we measured the
mean number of hydrogen bonds it makes to *any RNA atom* per simulation
frame, using mdtraj's Baker–Hubbard criterion (donor–acceptor distance
≤ 3.5 Å and donor–H–acceptor angle within standard cutoffs — the same
geometric definition used by GROMACS `gmx hbond`).

| feature | what it measures | AUROC |
|---|---|---:|
| **`hbond_RNA_F156`** | mean H-bonds per frame from Phe156 (ZF1 aromatic-stack residue) to RNA | **0.944** |
| **`hbond_RNA_R160`** | mean H-bonds per frame from Arg160 (ZF1 base-recognition Arg) to RNA | **0.833** |
| `hbond_RNA_H170` | mean H-bonds from His170 (ZF1 Zn-ligand) to RNA | 0.611 |
| `hbond_RNA_S185` | mean H-bonds from Ser185 to RNA | 0.583 |
| `hbond_RNA_R184` | mean H-bonds from Arg184 (ZF2 base-recognition Arg) to RNA | 0.556 |
| `hbond_RNA_Y154` | mean H-bonds from Tyr154 (ZF1 aromatic clamp) to RNA | 0.500 |
| `hbond_RNA_T187` | mean H-bonds from Thr187 to RNA | 0.472 |
| `hbond_RNA_C206` | mean H-bonds from Cys206 (ZF2 Zn-ligand) to RNA | 0.389 |
| `hbond_RNA_G204` | mean H-bonds from Gly204 (ZF2 hinge) to RNA | 0.333 |
| `hbond_RNA_C212` | mean H-bonds from Cys212 (ZF2 Zn-ligand) to RNA | 0.333 |
| `hbond_RNA_C208` | mean H-bonds from Cys208 to RNA | 0.222 |

**Headline:** the **F156 and R160 residue-level H-bond rates** are by
far the strongest discriminators across the 9 mutants, and both beat
the best AF3 feature (`pair_ipTM` = 0.778).

---

## 3. The "tier" picture

Ranked, with the AF3 best of 0.778 marked as the reference:

```
  F156 H-bond to RNA  ████████████████████████████████████████████████████████  0.944   ← MD wins
  R160 H-bond to RNA  ████████████████████████████████████████████████          0.833   ← MD wins
  AF3 pair-ipTM       ███████████████████████████████████████████               0.778   ← best AF3
  AF3 PAE cross-chain ███████████████████████████████████████                   0.722
  H170 H-bond to RNA  ████████████████████████████████                          0.611
  S185 H-bond to RNA  ███████████████████████████                               0.583
  RNA-P RMSD          ███████████████████████                                   0.556
  ...
  AF3 pLDDT-ZF1       ██████████████████                                        0.500   = random
  ...
  total contacts      █████████████                                             0.333
  Cα RMSD             ██████████                                                0.278
```

The two clear winners are MD-derived per-residue features at F156 and
R160; everything else is below the AF3 ipTM baseline. The exact
ranking is in `plots/fig_feature_aurocs.png`.

---

## 4. What each "winning" feature actually means biologically

### F156 — the aromatic stacking clamp of ZF1
Phe156 (with Tyr154 next door) forms the well-known **aromatic-stacking
clamp** that grips a uracil base of the ARE element in TTP-family RNA
binding.  It contributes ≈0.3 H-bonds per frame in WT (a backbone-carbonyl
H-bond to the phosphate backbone); in deleterious variants this drops
by 45 % or more.  F156 is *not* mutated in any of our variants —
the H-bond loss at F156 is a *secondary* allosteric effect of mutations
elsewhere (R184C, R160W).  This is exactly the kind of signal MD can
detect that AF3 cannot.

### R160 — the conserved base-recognition arginine of ZF1
Arg160 makes a salt-bridge / H-bond network to two uracil bases of the
ARE.  In WT it sustains ≈1.0 H-bond per frame to RNA.  The biggest
losses are:
- **R160W** (mutation removes the guanidinium → drops to 0.60/frame)
- **R184C** and **Y154H** (allosteric / fold perturbation; drop to
  0.40–0.47/frame)
- **C206F** drops to 0.88/frame (mild)

This matches the classic TTP/ZFP36L2 RNA-binding literature.

### Cα RMSD, RNA-P RMSD, contact count
These are bulk metrics. They tell you:
- **Cα RMSD (whole-protein)**: how much the entire 70-aa zinc-finger
  domain backbone has drifted from its starting position.
- **ZF1 / ZF2 Cα RMSD**: same idea but restricted to one finger.
- **RNA-P RMSD**: how much the 8 phosphates of the 9-nt RNA have moved
  relative to the start — a direct measure of "did the RNA stay
  bound?".
- **Heavy-atom contacts (4 Å)**: how many close-contact atom pairs
  there are between protein and RNA in the last 20 % of the run.

On their own, these bulk metrics didn't separate del vs neu well
(AUROC 0.28–0.56) because they are dominated by motions of irrelevant
residues at the periphery of the domain.  They become useful when
combined with a residue-specific signal.

---

## 5. How to reproduce

```bash
# from /Users/eric/Desktop/ZFP36L2 AF3
python3 expanse_md_analysis/run_full_analysis.py
python3 expanse_md_analysis/make_auroc_figure.py
```

Outputs:
- `expanse_md_analysis/feature_aurocs.csv`
- `expanse_md_analysis/plots/fig_feature_aurocs.png`
- `expanse_md_analysis/plots/fig_drift_deltas_vs_WT_common139ns.png`
- `expanse_md_analysis/plots/fig_hbonds_per_residue_common139ns.png`

To change the analysis window, edit `COMMON_WINDOW_NS` in
`plot_drift_and_hbonds.py` and rerun the driver.
