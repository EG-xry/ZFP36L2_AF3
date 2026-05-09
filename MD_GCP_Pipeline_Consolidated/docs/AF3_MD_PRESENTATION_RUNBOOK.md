# AF3 + MD Presentation Runbook

## Purpose

Generate a presentation-ready bundle that combines:
- existing AF3 static analyses (`Data/`, `results/`)
- current cloud minimal-MD outputs (`manifests/`)

The bundle includes new figures, merged tables, and two markdown summaries.

## Command

From workspace root:

```bash
cd "/Users/eric/Desktop/ZFP36L2 AF3"
.venv/bin/python MD_GCP_Pipeline_Consolidated/scripts/utils/build_af3_md_presentation_bundle.py
```

## Output location

Default output folder:

```text
/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated/presentation/compiled_YYYYMMDD
```

For this run:

```text
/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated/presentation/compiled_20260428
```

## Main artifacts

- `REPORT.md`
- `AF3_STATIC_EXPANSION_PLAN.md`
- `tables/combined_af3_md_selected.csv`
- `tables/af3_variant_scores_selected.csv`
- `figures/figure01_af3_pair_metrics_by_variant.png`
- `figures/figure02_af3_model0_run_stability.png`
- `figures/figure03_af3_class_separation.png`
- `figures/figure04_af3_md_bridge.png`
- `figures/figure05_md_qc_overview.png`

## Variant scope used

Included:
- WT
- Y154H
- R160W
- R184C
- S185N
- T187S
- G204D
- C206F
- C208A
- D210E

Excluded:
- R160E
- H212A

## Notes

- C208A is retained as a control and has fewer AF3 models (`n=25`) than most others (`n=40`).
- AF3 classification statistics are exploratory due small variant-level sample size.
