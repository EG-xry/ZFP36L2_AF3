# Commands Quickstart

Run from:

```bash
cd "/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated"
```

## 1) Rebuild variant manifests and staged inputs

```bash
python3 scripts/utils/build_variant_manifest.py
python3 scripts/gcp/stage_variant_inputs.py --mode copy
python3 scripts/gcp/build_task_environments.py
```

Default manifest excludes `R160E` and `H212A`.

## 2) Check local runtime readiness

```bash
bash scripts/utils/check_local_md_stack.sh
```

## 3) Core artifacts to inspect

```bash
ls manifests
ls runs/variants
head -n 20 manifests/variant_run_manifest.csv
head -n 20 manifests/variant_order.txt
```

## 4) Optional syntax smoke test

```bash
python3 -m py_compile \
  scripts/utils/build_variant_manifest.py \
  scripts/gcp/stage_variant_inputs.py \
  scripts/gcp/build_task_environments.py
```

## 5) Post-run output completeness summary

```bash
python3 scripts/utils/summarize_variant_outputs.py
head -n 20 manifests/variant_output_qc.csv
```

## 6) Cloud auth + WT smoke test

```bash
bash scripts/gcp/gcloud_cli_docker.sh version
bash scripts/gcp/gcloud_cli_docker.sh auth list
bash scripts/gcp/auth_login_via_docker.sh
export PROJECT_ID=\"YOUR_PROJECT_ID\"
export BUCKET_NAME=\"YOUR_UNIQUE_BUCKET_NAME\"
export REGION=\"us-central1\"
bash scripts/gcp/run_wt_smoke_cloud.sh
```

Use the project ID string (for example `project-2f27e3aa-de6c-4f65-ad4`), not the numeric project number.

## 7) Run all selected variants on Batch (10 tasks)

```bash
export PROJECT_ID="YOUR_PROJECT_ID"
export REGION="us-west2"
export BUCKET_NAME="YOUR_BUCKET"

# Optional runtime knobs:
export PARALLELISM="3"
export STAGE1_NS="0.01"
export STAGE2_NS="0.02"
export PRODUCTION_NS="0.05"
export MAX_RUN_DURATION="14400s"

bash scripts/gcp/run_variants_cloud.sh
```

This script automatically:

1. Refreshes manifests and staged variant inputs.
2. Uploads `runs/variants` to `gs://$BUCKET_NAME/variants`.
3. Builds/pushes the image.
4. Renders and submits one Batch job with one task per variant.

## 8) Finalize run + generate QC + professional figures

```bash
export PROJECT_ID="YOUR_PROJECT_ID"
export REGION="us-west2"
export JOB_NAME="YOUR_BATCH_JOB_NAME"
export BUCKET_NAME="YOUR_BUCKET"

bash scripts/gcp/postrun_finalize_and_build_reports.sh
```

This command waits for the Batch job to finish (default), syncs outputs from GCS,
recomputes MD summaries + structural drift tables, and rebuilds the AF3+MD
presentation bundle.
