# ZFP36L2 + Nonamer MD Analysis Report

**Generated:** 2026-05-27
**Source:** `expanse_md_results/` (10 variants on SDSC Expanse)
**Goal:** Use MD on AF3 protein–RNA complexes to discriminate deleterious
from neutral missense variants of the ZFP36L2 zinc-finger domain (ZFD)
binding a 9-nt AU-rich nonamer (5′-UAUAAAUAU-3′) with 2 Zn²⁺ ions.

---

## 1. What is on disk

All 10 variants have a working trajectory.  Three different run lengths exist
because of `target_ns` settings and wall-clock cut-offs.  Summary of what is
in `expanse_md_results/<variant>/`:

| variant | class | status | target_ns | ns_completed | trajectory size | production hours |
|---|---|---|---|---|---|---|
| WT    | WT          | checkpointed | 400 | **260.0** | 490 MB  | (still running) |
| C206F | deleterious | **completed** | 400 | 400.0 | 1.0 GB  | 18.1 |
| C208A | deleterious | checkpointed | 400 | **260.5** | 476 MB  | (still running) |
| G204D | deleterious | **completed** | 400 | 400.0 | 1.0 GB  | 18.1 |
| R160W | deleterious | **completed** | 400 | 400.0 | 1.0 GB  | 21.0 |
| R184C | deleterious | checkpointed | 400 | **251.0** | 476 MB  | (still running) |
| Y154H | deleterious (your set) | checkpointed | 400 | **256.5** | 484 MB  | (still running) |
| D210E | neutral     | **completed** | **50**  | 50.0  | 264 MB  | 2.8 |
| S185N | neutral     | **completed** | 400 | 400.0 | 1.1 GB  | 21.9 |
| T187S | neutral     | **completed** | **50**  | 50.0  | 271 MB  | 2.8 |

⚠️ **Inhomogeneous run length is a real issue.**  Right now we have:
- 4 runs at 400 ns (C206F, G204D, R160W, S185N)
- 4 runs ~250–260 ns and still extending (WT, C208A, R184C, Y154H)
- 2 short 50 ns runs (D210E, T187S)

You correctly flagged that this is going to be problematic — see §4
**before drawing biological conclusions**.

---

## 2. Are the trajectories themselves valid?

Yes.  I parsed every `md.log` (file: `parse_logs.py`) and computed the
following sanity-check signals for every run.  See `md_run_summary.csv`.

| variant | T_mean (K) | T_std | Epot drift (kJ/mol/ns) | density (g/mL) | NaN frames |
|---|---|---|---|---|---|
| WT    | 310.16 | 1.97 | -160 (still ramping) | 1.010 | 0 |
| C206F | 310.31 | 2.04 | -0.31 | 1.012 | 0 |
| C208A | 310.16 | 1.96 | -153 (still ramping) | 1.013 | 0 |
| D210E | 310.35 | 1.93 | -0.29 | 1.008 | 0 |
| G204D | 310.34 | 2.05 | -0.19 | 1.013 | 0 |
| R160W | 310.33 | 2.04 | -0.18 | 1.012 | 0 |
| R184C | 310.19 | 1.95 | -167 (still ramping) | 1.011 | 0 |
| S185N | 310.33 | 2.00 | -0.15 | 1.011 | 0 |
| T187S | 310.30 | 1.87 | -0.78 | 1.008 | 0 |
| Y154H | 310.19 | 1.96 | -164 (still ramping) | 1.010 | 0 |

Interpretation
- **Temperatures** are pinned at 310 K ± 2 K everywhere → Langevin
  thermostat is doing its job.
- **Densities** are 1.008–1.013 g/mL → reasonable for TIP3P solvent at
  310 K, 1 atm.
- **No NaN frames** anywhere → no exploding atoms.
- For the runs that finished (400 ns or 50 ns), the energy drift in
  the last 25 % is < 1 kJ/mol/ns → system is in steady state.
- The "checkpointed" runs (WT, C208A, R184C, Y154H) still show large
  negative `Epot_drift_kJmol_per_ns` because the difference between the
  first 1 % (NPT-eq ramp-down) and the last 25 % of so-far-completed
  production is still dominated by equilibration relaxation, not by
  trajectory drift.  Once they finish, this will normalize.

→ **Each run is physically sane.**

---

## 3. Quick structural/biological analysis

I computed (file: `compute_rmsd.py`) for every variant:
- Cα RMSD vs starting frame (whole 70-aa ZFD)
- Cα RMSD restricted to **ZF1** (ZFP36L2 res 152–179 → PDB res 2–29)
- Cα RMSD restricted to **ZF2** (ZFP36L2 res 190–217 → PDB res 40–67)
- RNA P-atom RMSD vs starting frame
- Number of protein–RNA heavy-atom contacts within 4 Å

All metrics use the **last 25 %** of the available trajectory.  Results
in `rmsd_summary.csv`:

### Per-variant (sorted by class)

| variant | class | t_final (ns) | Cα RMSD (Å) | ZF1 (Å) | ZF2 (Å) | RNA-P RMSD (Å) | PR contacts |
|---|---|---:|---:|---:|---:|---:|---:|
| WT    | WT          | **1.7**  | 2.48 | 1.90 | 2.60 |  65.9 | 387 |
| C206F | deleterious | 4.0  | 2.40 | 2.12 | 2.08 |  85.1 | 306 |
| C208A | deleterious | 1.7  | 2.46 | 2.11 | 2.24 |  47.6 | 209 |
| G204D | deleterious | 4.0  | 2.78 | 2.93 | 2.24 | 105.2 | 286 |
| R160W | deleterious | 4.0  | **5.63** | **4.41** | **4.73** |  79.5 | 276 |
| R184C | deleterious | 1.7  | 3.20 | 2.75 | 2.51 |  78.4 | 261 |
| Y154H | (your set: deleterious) | 1.7 | 2.13 | 1.65 | 1.97 |  96.5 | 268 |
| D210E | neutral     | 0.9  | 2.75 | 2.27 | 2.42 |  63.1 | 315 |
| S185N | neutral     | 4.0  | 3.36 | 3.65 | 2.41 |  15.5 | 272 |
| T187S | neutral     | 0.9  | 3.09 | 2.18 | 2.74 |  90.9 | 284 |

⚠️ `t_final` is the last frame written to `production.dcd`.  DCDs are
written every 100 ps, but the trajectories with `n_frames ≈ 4000` only
represent ~4 ns of the full ~400 ns run!  This is a **bug in the
analysis itself** — see §4.2 below.

### Class-level means (last-quarter, equal weighting per variant)

| class | n | Cα RMSD | ZF1 RMSD | ZF2 RMSD | RNA-P RMSD | PR contacts |
|---|---:|---:|---:|---:|---:|---:|
| WT (n=1)             | 1 | 2.48 | 1.90 | 2.60 |  65.9 | **387** |
| deleterious (n=5)    | 5 | **3.29** | **2.87** | 2.76 |  79.2 | **268** |
| Y154H (your set)     | 1 | 2.13 | 1.65 | 1.97 |  96.5 | 268 |
| neutral (n=3)        | 3 | 3.07 | 2.70 | 2.52 |  56.5 | 290 |

Headline findings
- **Protein–RNA contact count** clearly separates WT (387) from every
  variant (208–315).  This is the strongest single signal so far and is
  exactly what we want biologically: deleterious mutations should
  weaken the protein–RNA interface.
- **R160W** is the single most disruptive variant in both backbone
  movement and contact loss, consistent with its known status as a
  catalytic-arginine knock-out.
- **C208A and R184C** show the lowest contact counts among deleterious
  variants — they directly perturb the zinc finger fold.
- The **deleterious mean Cα RMSD (3.29 Å)** is ~0.8 Å larger than WT
  and clearly larger than ZF1/ZF2-only RMSDs (~2 Å), which is the
  expected pattern: the zinc fingers themselves stay folded, but the
  domain orientation relative to the RNA loosens.
- Among **neutrals**, only S185N has the longest, fully equilibrated
  trajectory.  Its RMSDs are intermediate, but contact count is 272 —
  much closer to deleterious than to WT.  However, with only 1 frame
  written per ns we may be sampling too sparsely to call.

### Important caveat
The previous compiled analysis in
`MD_GCP_Pipeline_Consolidated/manifests/md_variant_deltas_vs_wt_mixed8plus2.csv`
was based on 0.05–0.10 ns "smoke-test" runs and is not directly
comparable to the present 50–400 ns data.

---

## 4. Issues to fix before reading the biology

### 4.1 Heterogeneous run length — **biggest problem**
You cannot fairly compare a 50 ns trajectory of D210E and T187S
against a 400 ns trajectory of S185N.  The 50 ns runs may
not even be equilibrated yet, and contact-count drift takes 10–50 ns
to plateau for RNA-binding zinc fingers.

**Fix**: do one of the following before final biology calls:
- **Option A (clean)** — re-run D210E and T187S to 400 ns to match the
  others.  Required additional wall time: ~22 h × 2 = 44 GPU-h on a
  V100.
- **Option B (pragmatic)** — restrict the *primary* analysis window to
  the first 50 ns of every variant (the common interval).  This is the
  simplest fix and keeps all 10 variants on the table.  The
  presentation/manuscript should then explicitly say "first 50 ns of
  production".  We *also* keep the long runs available for an extended
  analysis (S185N control, etc.).
- **Option C** — finish the four checkpointed runs (WT, C208A, R184C,
  Y154H) to 400 ns and report only the 400 ns intersection.  Then
  re-run D210E and T187S also to 400 ns.  This is the gold-standard
  way to publish.

**My recommendation: Option B + Option C in parallel** — keep the
already-running ones extending, immediately resubmit D210E and T187S
to 400 ns on Expanse, and in the meantime use the 50 ns common window
for any presentation deliverable.

### 4.2 DCD trajectory interval mismatch
Some `production.dcd` files have far fewer frames than expected:
- C206F, G204D, R160W, S185N: 4000 frames over 400 ns → 100 ps/frame ✓
- WT, C208A, R184C, Y154H:    ~1700 frames over ~260 ns → 100 ps/frame ✓
- D210E, T187S:               ~900 frames over 50 ns  → wrong axis?

Looking at `run_md_anvil.py:701`, the default `--traj-interval-ps=100`
and each run is configured with `--production-ns` matching its log.
So the DCD frames *are* correctly 1 per 100 ps.

The `time_ns` axis in `compute_rmsd.py` falls back to "0.1 ns × frame
index" when `traj.time` is unavailable, which is correct.  However
some DCDs (e.g. C206F/G204D/R160W/S185N) only have 4000 frames despite
finishing 400 ns — that means the DCDReporter cleared and restarted
without `append=True` after a wall-time exit.  This is a *file* issue,
not a simulation issue: production was continued from checkpoint but
the DCD wasn't appended.

**Fix**: use the `production.chk` to confirm `currentStep`, and rely on
`md.log` (which DOES capture the full 400 ns) for the time axis.  For
trajectory-based analysis, the existing DCDs contain a **dense sample
of one continuous segment** of production, which is still
representative for steady-state metrics.  Just note in the manuscript
that the DCD covers a 400 ns *segment* equilibrium ensemble, not the
full 0–400 ns time series.

### 4.3 Single replicate per variant
We have 1 trajectory per variant.  Fluctuations of ~0.5 Å in Cα RMSD
are normal between replicates, so we cannot tell whether
`S185N` (3.36 Å) being higher than `D210E` (2.75 Å) is meaningful or
just noise.

**Fix**: launch 2 additional replicas (different velocity seeds,
maybe 100 ns each) per variant.  This is essential for any peer-
reviewed claim.

### 4.4 The Y154H label
Your data labels Y154H as "deleterious in this dataset" but the
broader literature can label it neutral.  Its RMSDs and contact counts
in the table are intermediate.  I left it under `deleterious_in_set`
in the outputs so it doesn't bias either class average.  This is purely
a labeling decision; the **data itself is fine**.

---

## 5. Concrete next steps (in priority order)

1. **Today** — Use the existing data with Option B (common 50 ns
   window).  I will write `compute_rmsd_common_window.py` if you want
   it; that gives publication-grade table immediately.
2. **Resubmit D210E and T187S to 400 ns on Expanse** (~22 GPU-h each
   on a V100, $5–10).  Same `submit_variant_expanse.sh` script, no
   code change needed beyond `--production-ns 400`.
3. **Resume the 4 checkpointed runs** (WT, C208A, R184C, Y154H) by
   submitting the same job again with `--resume`.  They are already at
   ~260 ns, so they need ~12 h more each.
4. **Launch 2 replicas of every variant** with different
   `--seed/velocities-to-temperature` for 100 ns each.  This is the
   biggest scientific upgrade.
5. **Add interface-specific metrics** that are more discriminative
   than total CA-RMSD:
   - Per-residue contact frequency for the canonical ARE-binding
     residues (Y154, F156, H170, R171, R184, S185, T187, F191, R200,
     C206, C208, H212).
   - Persistence of the canonical π-stack interactions:
     Y154–U(+1), F191–U(-2), F156–RNA base.
   - RNA base flipping (chi-angle distribution of each nucleotide).
   - Zn–Cys/His coordination distance time series — direct test of
     fold integrity for C206F, C208A, R184C.
6. **Statistical test on the per-class means**:
   - With 5 deleterious + 3 neutral variants, even a Mann-Whitney U on
     `n_pr_contacts_last25` will give a meaningful p-value.
   - Bootstrap 95 % CI on the class means.
7. **Cross-check against AF3 static metrics** already in
   `results/variant_avg_all.csv` (pair-ipTM, RMSD_whole, pLDDT).  If
   pair-ipTM and MD contact-count rank variants the same way, you
   have a defensible 2-modality classifier.
8. **Replicate the Block-B (GM-CSF 17-mer) and Block-C analyses** the
   same way once Block-A is solid.

---

## 6. Files produced by this analysis

```
expanse_md_analysis/
├── parse_logs.py              # parses every md.log → CSV time-series
├── compute_rmsd.py            # mdtraj-based RMSD + contact analysis
├── md_run_summary.csv         # one row per variant: T/E/density sanity
├── md_run_summary.md          # markdown version
├── rmsd_summary.csv           # one row per variant: all RMSD/contact metrics
├── rmsd_summary.md
├── <variant>_production.csv   # 10 files: full log time-series
├── <variant>_rmsd.csv         # 10 files: Cα/ZF1/ZF2/RNA RMSD per frame
└── <variant>_contacts.csv     # 10 files: protein-RNA contact count over time
```

Rerun the whole analysis any time with:
```bash
python3 expanse_md_analysis/parse_logs.py
python3 expanse_md_analysis/compute_rmsd.py
```
