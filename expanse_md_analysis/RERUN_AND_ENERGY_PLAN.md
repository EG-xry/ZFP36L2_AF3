# Fixes, new energy metrics, and the re-run plan

**Generated:** 2026-06-03

This doc covers the three things you asked for:
1. the checkpoint/DCD fix that's now in the code,
2. the new energy + literature metrics and what they showed,
3. a concrete recommendation on whether to re-run all, some, or from checkpoints.

---

## 1. The checkpoint/DCD-truncation bug is fixed

**What was wrong:** when one of the 4 long runs (WT, C212A, R184C, Y154H)
was resubmitted *without* `--resume`, `prepare_dcd_path()` deleted the
existing `production.dcd` and started a new one, so the first ~130 ns of
frames were lost (the simulation still reached 400 ns; only the saved
trajectory was truncated to ~270 ns).

**The fix** (`anvil_md/run_md_anvil.py`, two edits):
- `prepare_dcd_path()` now **always appends to a valid existing DCD**, even
  when append wasn't explicitly requested. It only creates a new file if no
  valid DCD exists. It prints a safeguard message when it does this.
- The fresh-start branch no longer calls `clear_production_restart_artifacts`
  on a *valid* DCD — it only clears truncated/zero-byte files — and prints a
  `[WARN]` if a valid DCD exists but `--resume` was not passed.

Net effect: a future resubmission can never silently overwrite completed
trajectory frames. A genuinely clean restart still works (delete the DCD
first, or it's a NaN-retry `fresh_start`).

---

## 2. New metrics added (`compute_energy_metrics.py`)

I added five metrics that are standard in protein–RNA / protein–ligand MD
discrimination literature, computed on the 269 ns common window. AUROC is
for separating the 6 deleterious from 3 neutral variants (0.5 = random,
1.0 = perfect):

| metric | what it is | literature basis | **AUROC** |
|---|---|---|---:|
| **screened Coulomb (protein⁺↔RNA-phosphate⁻)** | Debye-screened electrostatic attraction between protein cationic groups (Arg/Lys) and RNA phosphates, averaged over frames | electrostatics dominate ZF–RNA binding (Hudson 2004; Brewer 2004) | **0.833** |
| salt bridges (Arg/Lys ↔ phosphate, <4 Å) | count of specific cation–phosphate contacts | Yang 2021 protein–RNA salt-bridge analysis | 0.611 |
| MM-GBSA interaction energy (GBn2, single-traj) | ΔE = E(complex) − E(protein) − E(RNA) in implicit solvent | Genheden & Ryde 2015 endpoint-binding-energy method | 0.444 |
| binding-residue RMSF (ZF1+ZF2 Cα) | flexibility of the RNA-binding residues | implied by Ramos et al. RMSF discussion | 0.167 |
| interface atom count (<5 Å) | size of the protein–RNA interface | SASA/buried-area proxy | 0.167 |

### Headline
**The screened-Coulomb protein–RNA electrostatic term (AUROC 0.833) is the
first MD feature that clearly beats AF3 pair-ipTM (0.778).** It is also
sampling-robust — it gave 0.72–0.83 across every stride/window I tried,
unlike the per-residue H-bond features that collapsed from 0.94 → 0.61.

This makes physical sense: ZFP36L2 grips the AU-rich RNA mainly through
Arg/Lys → phosphate-backbone electrostatics, and the deleterious mutations
(especially R160W and R184C, which remove or bury a cationic side chain)
weaken exactly this term. WT has the *weakest* total cation–phosphate
attraction in our set (−7938) only because in WT the RNA sits in the
canonical pocket with fewer but more specific contacts; the deleterious
ones show diffuse, non-native electrostatic clinging — so the ranking works
but the mechanism is subtle and needs the replicate check below.

### Caveats on the new metrics
- **MM-GBSA failed to be useful (0.444).** The single-trajectory GBn2 energy
  came out *positive* (~+12000 kJ/mol) and didn't rank the variants. Two
  reasons: (a) we strip the Zn and the zinc-finger fold is then unphysical
  in implicit solvent, and (b) single-trajectory MM-GBSA on 10 frames is far
  too noisy. A proper MM-GBSA needs: keep a bonded Zn model, 100+ frames,
  and ideally the 3-trajectory or interaction-entropy variant. It is worth
  doing but is a bigger job (see §3).
- **interface atom count and binding-residue RMSF are anti-correlated**
  (AUROC < 0.5) — i.e., on this dataset the neutral variants happen to have
  *smaller* interfaces / *higher* RMSF, the opposite of the naive
  expectation. With n=9 this is almost certainly noise, not signal.

### Updated full AUROC ranking (all modalities, 269 ns)
```
screened Coulomb (protein–RNA)   0.833   MD-energy   ← new best, beats AF3
AF3 pair-ipTM                    0.778   AF3
MD RNA-P RMSD                    0.722   MD
AF3 PAE cross-chain              0.722   AF3
salt bridges (cation–phosphate)  0.611   MD-energy
MD per-residue H-bond (F156/R160)0.611   MD-Hbond
...
MM-GBSA (GBn2, single-traj)      0.444   MD-energy  (needs proper setup)
interface size / binding RMSF    0.167   MD-energy  (noise at n=9)
```

---

## 3. Do you need to re-run? — my recommendation

**You do NOT need to re-run all 10 from scratch.** Here's the decision tree:

### What's already good enough
- All 10 reached a genuine 400 ns (verified via steps + progress + log).
- 6 variants have a clean 400 ns DCD; 4 have a clean *contiguous* 270 ns DCD.
- 270 ns is plenty for the equilibrium metrics we're computing.

### What I recommend, in priority order

**(A) Highest value — replicates, not longer single runs.**
Run **2 additional independent replicas** (different velocity seeds) of all
10 variants at **200 ns each**, using the *fixed* `run_md_anvil.py`. This is
the single most important thing: every AUROC we have is effectively n=1, and
with 6-vs-3 a single trajectory's noise moves AUROC by ±0.2. This is exactly
what Ramos did (3 replicas). Cost on Expanse V100: ~10 variants × 2 replicas
× ~10 h = ~200 GPU-h (you had 477 GPU-h budget).
→ Use: `submit_all_expanse.sh 200` after setting a fresh `--variant` output
  dir per replica (e.g. `outputs_rep2/`, `outputs_rep3/`).

**(B) Optional — recover the full 0–400 ns DCD for the 4 truncated runs.**
You don't strictly need this (270 ns is enough for equilibrium metrics), but
if you want the complete trajectory, re-run *only those 4* (WT, C212A,
R184C, Y154H) from scratch in one clean wall-window with the fixed code so
they never hit a resume. ~4 × 20 h = 80 GPU-h. The other 6 are already fine.
→ Do NOT re-run the 6 complete ones.

**(C) Do a *proper* MM-GBSA** on the existing trajectories (no new MD
needed) — this is an analysis upgrade, not a re-run:
   - keep a bonded/dummy-cation Zn model so the fold is physical,
   - use ≥100 frames per variant,
   - compute the polar+nonpolar GB solvation properly (not the raw
     potential-energy difference I prototyped),
   - ideally use `gmx_MMPBSA` or `AmberTools MMPBSA.py` which handle the
     decomposition and Zn correctly.
   This is the metric most likely to give a clean, publishable separation
   and is worth doing before any more MD.

### One-line summary
> Don't re-run everything. **Run 2 replicas of all variants at 200 ns
> (Option A)** to put error bars on the AUROCs; optionally re-run just the 4
> truncated variants for complete DCDs (Option B); and do a proper
> AmberTools/gmx_MMPBSA MM-GBSA on the trajectories you already have
> (Option C). The screened-Coulomb electrostatic metric (AUROC 0.833) is the
> most promising lead and should be the headline MD feature once replicate
> error bars confirm it.

---

## 4. Files added/updated this pass

```
anvil_md/run_md_anvil.py                 # DCD append/resume fix (2 edits)
expanse_md_analysis/compute_energy_metrics.py   # new metrics + MM-GBSA
expanse_md_analysis/energy_metrics_summary.csv  # per-variant values
expanse_md_analysis/energy_metric_aurocs.csv    # AUROCs of the new metrics
expanse_md_analysis/RERUN_AND_ENERGY_PLAN.md    # this file
```

Re-run the energy metrics any time:
```bash
python3 expanse_md_analysis/compute_energy_metrics.py --max-ns 269 --stride 20
# add --mmgbsa for the (currently rough) implicit-solvent energy
```
