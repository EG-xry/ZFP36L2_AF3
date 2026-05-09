# GCP Post-Run Analysis Runbook

This runbook generates all post-cloud outputs as soon as a Batch variant job finishes:

- sync cloud outputs to local workspace
- recompute MD QC summary tables
- recompute structural drift tables (protein + RNA)
- rebuild presentation bundle (figures + report + AF3 static analysis)

## Script

Use:

```bash
bash scripts/gcp/postrun_finalize_and_build_reports.sh
```

## Defaults

- `PROJECT_ID=project-2f27e3aa-de6c-4f65-ad4`
- `REGION=us-west2`
- `JOB_NAME=variants-md-12hmax-20260428-190901`
- `BUCKET_NAME=zfp36l2-md-20260428-12676`
- waits for job completion by default (`WAIT_FOR_JOB=true`)

## Common usage

From:

```bash
cd "/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated"
```

Then run with explicit values (recommended):

```bash
export PROJECT_ID="project-2f27e3aa-de6c-4f65-ad4"
export REGION="us-west2"
export JOB_NAME="variants-md-12hmax-20260428-190901"
export BUCKET_NAME="zfp36l2-md-20260428-12676"
bash scripts/gcp/postrun_finalize_and_build_reports.sh
```

## Optional controls

- `WAIT_FOR_JOB=false` to skip waiting and process current outputs immediately.
- `ALLOW_PARTIAL=true` to continue even if job state is `FAILED`.
- `POLL_SECONDS=120` to adjust polling interval while waiting.
- `PRESENTATION_DIR=/custom/path` to control bundle output location.

## Outputs updated

- `manifests/cloud_variant_analysis.csv`
- `manifests/cloud_variant_structural_drift.csv`
- `manifests/cloud_variant_structural_drift_rna.csv`
- `presentation/compiled_YYYYMMDD/` (report + professional figures)

## Figure set generated

- `figure01_af3_pair_metrics_by_variant.png`
- `figure02_af3_model0_run_stability.png`
- `figure03_af3_class_separation.png`
- `figure04_af3_md_bridge.png`
- `figure05_md_qc_overview.png`
