# GCP Multi-Variant Runbook

This runbook launches one Batch task per selected variant:

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

Excluded by design: `R160E`, `H212A`.

## Prerequisites

Run from:

```bash
cd "/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated"
```

Authenticate if needed:

```bash
bash scripts/gcp/auth_login_via_docker.sh
```

## Launch command

```bash
export PROJECT_ID="project-2f27e3aa-de6c-4f65-ad4"
export REGION="us-west2"
export BUCKET_NAME="zfp36l2-md-20260428-12676"

# Recommended GPU profile
export GPU_TYPE="nvidia-tesla-t4"
export GPU_COUNT="1"
export INSTALL_GPU_DRIVERS="true"
export OPENMM_PLATFORM="CUDA"
export MACHINE_TYPE="n1-standard-4"
export PARALLELISM="1"
export ALLOWED_LOCATIONS="us-west2-b,us-west2-c"

# Optional stage durations / timeout
export STAGE1_NS="0.01"
export STAGE2_NS="0.02"
export PRODUCTION_NS="0.05"
export MAX_RUN_DURATION="14400s"

bash scripts/gcp/run_variants_cloud.sh
```

## Runtime profiles

### GPU profile (recommended)

```bash
export GPU_TYPE="nvidia-tesla-t4"
export GPU_COUNT="1"
export INSTALL_GPU_DRIVERS="true"
export MACHINE_TYPE="n1-standard-4"
export OPENMM_PLATFORM="CUDA"
export PARALLELISM="1"
```

With current quota, `NVIDIA_T4_GPUS=1`, so variant tasks run one at a time unless GPU quota is increased.

### CPU fallback profile

```bash
export GPU_TYPE=""
export OPENMM_PLATFORM="CPU"
export MACHINE_TYPE="n2-standard-2"
export CPU_MILLI="2000"
export MEMORY_MIB="8192"
export PARALLELISM="6"
```

## What this script does

1. Validates Docker + gcloud wrapper access.
2. Enables required APIs.
3. Ensures IAM roles for the Batch runtime service account.
4. Rebuilds:
   - `manifests/variant_run_manifest.csv`
   - `runs/variants/*`
   - `manifests/batch_task_environments.json`
5. Uploads staged inputs to:
   - `gs://$BUCKET_NAME/variants/`
6. Builds and pushes runtime image to Artifact Registry.
7. Renders:
   - `manifests/batch_variants_job.json`
8. Submits one Batch job with one task per variant.

## Monitor

```bash
scripts/gcp/gcloud_cli_docker.sh batch jobs describe JOB_NAME --location "$REGION"
scripts/gcp/gcloud_cli_docker.sh batch tasks list --job JOB_NAME --location "$REGION"
scripts/gcp/gcloud_cli_docker.sh storage ls -r "gs://$BUCKET_NAME/outputs/variants/**"
```

To inspect quota bottlenecks:

```bash
scripts/gcp/gcloud_cli_docker.sh compute project-info describe \
  --project "$PROJECT_ID" \
  --flatten='quotas[]' \
  --format='table(quotas.metric,quotas.limit,quotas.usage)' | \
  rg 'CPUS_ALL_REGIONS|IN_USE_ADDRESSES|NVIDIA_'
```

## Output layout

Each task writes to:

- `gs://$BUCKET_NAME/outputs/variants/<VARIANT>/md.log`
- `gs://$BUCKET_NAME/outputs/variants/<VARIANT>/production.dcd`
- `gs://$BUCKET_NAME/outputs/variants/<VARIANT>/final_frame.pdb`
- `gs://$BUCKET_NAME/outputs/variants/<VARIANT>/task_metadata.env`

## Reliability note

The task runs MD on local VM scratch first and copies final outputs to GCS after completion. This avoids long-write instability from direct trajectory writes to gcsfuse-mounted paths.

Each task also normalizes raw AF3 variant PDBs with `clean_af3_structures.py` before running OpenMM, which avoids RNA terminal template mismatches observed on uncleaned inputs.
