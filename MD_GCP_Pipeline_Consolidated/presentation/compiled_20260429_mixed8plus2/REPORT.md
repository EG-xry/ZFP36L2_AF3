# AF3 + MD Compiled Results

Bundle folder: `MD_GCP_Pipeline_Consolidated/presentation/compiled_20260429_mixed8plus2`
Date: `2026-04-29`

## Scope

- Included variants: WT + 10 variants (includes controls C208A and R160E).
- Excluded variant: H212A.
- AF3 static source: existing multi-run outputs in `Data/*` and summary tables in `results/*`.
- MD source: cloud minimal-MD outputs in `MD_GCP_Pipeline_Consolidated/manifests/*`.
- MD provenance split used here: `8` from main batch + `2` from fallback short rerun (`S185N, C208A`).

## Key Findings

- AF3 composite deleterious-vs-neutral AUC (variant-level means): `0.857` (95% bootstrap CI `0.571` to `1.000`), permutation p-value `0.0563`.
- Minimal-MD completion status in current table: `10/10` variants marked `complete_minimal_md`.
- AF3 composite score is WT-referenced and directional (higher = worse than WT on interface confidence).

### Most WT-deviated variants by AF3 composite

- C206F: score `3.51` (deleterious)
- R160W: score `2.07` (deleterious)
- C208A: score `1.32` (deleterious)

### Closest to WT by AF3 composite

- Y154H: score `-0.06` (deleterious)
- WT: score `0.00` (WT)
- D210E: score `0.01` (neutral)

## Files

- `tables/af3_variant_summary_selected.csv`: canonical AF3 static means/SD (all models) for selected variants.
- `tables/af3_per_model_pair_metrics_selected.csv`: per-model pair metrics parsed from AF3 summary JSONs.
- `tables/combined_af3_md_selected.csv`: AF3 + MD merged table for cross-modal comparisons.
- `tables/md_run_provenance_selected.csv`: per-variant MD run-source provenance (main batch vs fallback rerun).
- `figures/figure01_af3_pair_metrics_by_variant.png`: interface confidence by variant.
- `figures/figure02_af3_model0_run_stability.png`: run-to-run spread of top AF3 model.
- `figures/figure03_af3_class_separation.png`: deleterious vs neutral separation.
- `figures/figure04_af3_md_bridge.png`: AF3 score vs MD drift.
- `figures/figure05_md_qc_overview.png`: MD runtime QC summary.

## Interpretation Limits

- AF3 classification check uses small N at the variant level (neutral n=3, deleterious n=7).
- Current MD is short-timescale minimal runs and not a full 1000 ns protocol replication.

