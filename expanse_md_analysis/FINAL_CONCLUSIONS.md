# Final conclusions — all 10 ZFP36L2 + nonamer MD runs complete (400 ns)

**Generated:** 2026-06-02 (numbers verified against the repo's
`feature_aurocs.csv` and `hbonds_per_residue_common269.csv`)
**Data:** every variant reports `status=completed`, `ns_completed=400.0`,
`steps_completed=200,000,000` (= 400 ns at 2 fs).
**Analysis window:** 269 ns common window (the continuous DCD floor; see §1).

> **Note on the two AUROC tables you may have seen.** An earlier terminal
> run printed F156 = 0.833 / R160 = 0.667. That table was produced while the
> legacy `hbonds_per_residue_common50.csv` alias still held the **139 ns**
> H-bond data. The *authoritative* current numbers, computed directly from
> the 269 ns H-bond file, are **F156 = 0.611, R160 = 0.611** (the 0.833 was
> stale). Both tables agree on the qualitative point: **no MD per-residue
> H-bond feature beats AF3 pair-ipTM (0.778) on the full trajectory.**

---

## 1. Run completion & the trajectory-file caveat (Q1)

All ten 400 ns runs finished. I verified this three independent ways:
- `progress.txt`: `steps_completed=200000000`, `ns_completed=400.0`,
  `status=completed` for **all 10**.
- `md.log` step column reaches **200,000,000 steps** (= 400 ns) for all 10.
- The 4 short ones each show exactly **one** "time reset" in the log — the
  signature of a single checkpoint-resume.

However the **DCD trajectory files only contain a continuous ~270 ns
segment** for the four runs that were resumed:

| variant | steps in log | DCD frames | DCD ns |
|---|---:|---:|---:|
| C206F, D210E, G204D, R160W, S185N, T187S | 2.0e8 (400 ns) | 4000 | 400 |
| WT, C212A, R184C, Y154H | 2.0e8 (400 ns) | ~2700 | ~270 |

### Exactly what happened (root cause)
Look at `anvil_md/run_md_anvil.py` lines 605–612: on a **fresh** (non-resume)
start the production DCD is opened with `append_traj=False`, and
`prepare_dcd_path()` (lines 287–297) **deletes any existing DCD**. The four
runs that were resubmitted were started *without* the `--resume` flag the
first time they were re-queued, so OpenMM wrote a brand-new `production.dcd`
that captured only the post-restart segment (~270 ns). The integrator and
the `production.chk` checkpoint genuinely advanced to 400 ns — only the
*trajectory frames* for the first ~130 ns were overwritten on disk.

### Can the missing first segment be recovered? (your Q1)
- **From the cloud / Expanse:** Only if a copy of the pre-restart
  `production.dcd` still exists in the job's scratch/output directory on
  Expanse *before* it was overwritten. It almost certainly is **not**
  recoverable — the file was unlinked in place, and Expanse `$SCRATCH` is
  purged. Check `~/zfp36l2_md/outputs/<V>/` and any `slurm-*.out` /
  per-job temp dirs, but do not expect it back.
- **The honest fix:** the `production.chk` checkpoints DO hold the full
  400 ns state. The clean way to get a *full* 0–400 ns trajectory is to
  **re-run with a corrected DCDReporter (`append=True` on resume)** — or,
  simplest, **re-run the four short variants from scratch in one
  wall-window** so they never hit a resume. Either way it's a re-run, not a
  recovery.
- For now, 269 ns is the largest *honest* common window and is what every
  number below uses.

---

## 2. The headline finding flipped — and that is itself informative

| feature | 50 ns AUROC | 139 ns AUROC | **269 ns AUROC** |
|---|---:|---:|---:|
| AF3 pair-ipTM (static) | 0.778 | 0.778 | **0.778** |
| MD RNA-P RMSD | 0.556 | 0.556 | **0.722** |
| AF3 PAE cross-chain | 0.722 | 0.722 | 0.722 |
| MD H170 H-bond | 0.389 | 0.611 | 0.667 |
| MD F156 H-bond | **0.778** | **0.944** | **0.611** ↓ |
| MD R160 H-bond | 0.722 | 0.833 | **0.611** ↓ |

**The F156 and R160 per-residue H-bond signal that looked spectacular at
50–139 ns does NOT hold up over the full trajectory.** At 269 ns:

- **AF3 pair-ipTM (0.778) is again the single best discriminator.**
- The only MD feature that *matches* AF3's runner-up is whole-RNA
  phosphate RMSD (0.722) — i.e., "did the RNA stay in place" — which is
  a bulk metric, not a per-residue mechanism.
- The per-residue H-bond features regress toward 0.5–0.66.

### Why the early signal was a mirage
Time-block analysis of R160's H-bond rate to RNA (mean H-bonds/frame):

| variant | class | 0–50 ns | 50–139 ns | 139–269 ns |
|---|---|---:|---:|---:|
| WT    | WT  | 1.12 | 0.94 | 1.16 |
| R160W | del | 1.00 | 1.15 | 1.12 |
| R184C | del | 1.22 | 1.15 | 1.14 |
| Y154H | del | 1.00 | 1.02 | **2.02** |
| G204D | del | 1.02 | 1.28 | **1.48** |
| C206F | del | 1.02 | 1.25 | 1.18 |
| D210E | neu | 1.14 | 1.17 | 1.14 |
| S185N | neu | 1.14 | 1.15 | 1.12 |
| T187S | neu | 1.10 | 1.26 | 1.26 |

The clean "deleterious lose R160 contact" picture of the first 50–139 ns
**washes out by 269 ns** — the deleterious variants' R160 contacts recover
or even exceed WT (Y154H climbs to 2.0, G204D to 1.5). In other words, the
early R160/F156 separation was a *transient relaxation artifact of the AF3
starting structure*, not a stable equilibrium property. When you sample
long enough, every variant's R160 finds *some* RNA contact because the
9-mer is floppy and re-binds in different registers.

**This is exactly why running to 400 ns mattered.** The 50–139 ns result
would have been a false-positive headline.

---

## 3. What actually holds up at 269 ns

1. **AF3 pair-ipTM remains the most reliable single classifier (0.778).**
   The static AF3 confidence score is genuinely informative and cheap.

2. **Whole-RNA-phosphate RMSD (0.722) is the best MD feature** and it
   improved with sampling (0.56 → 0.72). It measures whether the RNA
   stayed in its starting pose. This is a real, if modest, MD contribution.

3. **No single per-residue H-bond feature beats AF3 over the full
   trajectory.** The mechanism-localization story is *not* supported once
   you control for sampling.

4. **Bulk Cα/ZF RMSD and total contact count remain near-random**
   (0.39–0.56) — consistent across all windows.

### Per-variant 269 ns drift / contacts (for reference)

| variant | class | Cα RMSD (Å) | RNA-P RMSD (Å) | PR contacts |
|---|---|---:|---:|---:|
| WT    | WT  | 3.10 | 80.6 | 316 |
| C206F | del | 3.35 | 85.0 | 304 |
| C212A | Zn-ctrl | 1.90 | 60.3 | 275 |
| G204D | del | 2.74 | 94.5 | 235 |
| R160W | del | 4.06 | 83.0 | 288 |
| R184C | del | 4.92 | 83.0 | 294 |
| Y154H | del | 3.36 | 101.3 | 268 |
| D210E | neu | 3.14 | 59.4 | 251 |
| S185N | neu | 3.92 | 25.2 | 221 |
| T187S | neu | 2.07 | 90.6 | 368 |

Even RNA-P RMSD (the best MD feature) does not cleanly separate the
classes by eye — note T187S (neutral) has 90.6 Å while C212A (deleterious
control) has only 60.3 Å. The 0.722 AUROC comes from the *aggregate* rank
ordering, not a clean threshold.

---

## 4. What can we honestly conclude?

**Scientifically honest bottom line:**

> *With 1 replicate per variant at 400 ns, MD on AF3 ZFP36L2–nonamer
> complexes does NOT outperform AF3's own pair-ipTM confidence score at
> separating the 6 deleterious from the 3 neutral nsSNPs. The best MD
> feature (whole-RNA RMSD, AUROC 0.72) ties AF3's runner-up but not its
> best (0.78). The striking per-residue H-bond signal seen in shorter
> windows was a starting-structure relaxation transient that did not
> survive longer sampling.*

This is a **negative-leaning but valuable** result:
- It tells you the *cheap* method (AF3, seconds) is as good as the
  *expensive* method (400 ns × 10 GPU-days of MD) for this binary task,
  on this dataset, at this replicate count.
- It is a textbook reminder that **single-replicate MD features must be
  validated against longer sampling** before they are believed — we
  caught a false positive precisely because we ran to completion.

**What MD still uniquely provides** (qualitatively, not as a classifier):
- Mechanistic movies of *how* the RNA detaches (e.g., R184C shows the
  largest backbone Cα RMSD at 4.9 Å — the ZF2 fold is genuinely the most
  perturbed).
- The C212A Zn-knockout control behaves as expected (lowest Cα RMSD
  because the Zn restraint holds the fold, but reduced contacts).

---

## 5. Recommended next steps (in priority order)

1. **Fix the DCD append bug and re-extract the full 0–400 ns trajectory.**
   The `production.chk` checkpoints contain the full state; re-running with
   a corrected DCDReporter `append=True` (or stitching the segments) gives
   the true 0–400 ns ensemble. Until then we are analyzing a 269 ns window.

2. **Run ≥3 replicates per variant** with different seeds. Every AUROC in
   this report has an effective n of 1; with 6 vs 3 variants a single
   trajectory's noise can move an AUROC by ±0.2. This is the single most
   important fix and is exactly what Ramos et al. did (3 × 1000 ns).

3. **Add proper interface metrics that are sampling-robust:**
   - Fraction of native contacts (Q) relative to the AF3 / 1RGO interface
   - MM-GBSA / MM-PBSA binding free energy on the equilibrated ensemble
     (this is the metric most likely to actually beat AF3, because it is
     an energy, not a geometry snapshot)
   - RNA-register-aware contact maps (the floppy 9-mer rebinds in shifted
     registers, which confounds raw per-residue counts)

4. **Re-evaluate the labels.** D210E is actually D219E and C208A is C212A
   (see VARIANT_LABELS_AUDIT.md). Decide whether D219E belongs in the
   neutral set at all.

5. **Be careful about over-claiming.** Do NOT take the 139 ns F156=0.944
   number to the manuscript — it does not replicate at 269 ns. The
   defensible claim right now is "AF3 pair-ipTM classifies these variants
   as well as 400 ns single-replicate MD; longer/replicated MD or binding-
   free-energy calculations are needed to add value."

---

## 6. Files refreshed in this pass

```
expanse_md_analysis/
├── plot_drift_and_hbonds.py        (COMMON_WINDOW_NS = 269.0)
├── common269_summary.csv           NEW — 269 ns drift/contacts
├── hbonds_per_residue_common269.csv NEW
├── hbonds_key_residues_common269.csv NEW
├── feature_aurocs.csv              updated → AF3 pair-ipTM 0.778 back on top
├── af3_vs_md_combined.csv          updated
├── FINAL_CONCLUSIONS.md            this file
└── plots/
    ├── fig_feature_aurocs.png          updated (269 ns)
    ├── fig_drift_deltas_vs_WT_common269ns.png
    ├── fig_hbonds_per_residue_common269ns.png
    ├── fig_rmsd_full_trajectory.png    updated
    └── fig_trajectory_length_status.png updated (all 400 ns)
```
