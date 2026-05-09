# GCP WT Smoke Runbook

## What was done autonomously in this workspace

1. Moved source folders into workspace:
   - `/Users/eric/Desktop/ZFP36L2 AF3/attachments`
   - `/Users/eric/Desktop/ZFP36L2 AF3/AlphaFold`
2. Added Google Cloud helper scripts in `scripts/gcp/`.
3. Confirmed local `gcloud` binary is not installed.
4. Enabled a no-install path using Dockerized Google Cloud CLI (`gcloud_cli_docker.sh`).
5. Verified wrapper works by running `gcloud version` in container.
6. Added WT smoke Batch job renderer and end-to-end submit script.
7. Patched project-ID, path-with-spaces, and service-account permission issues in submit flow.
8. Patched Batch GCS mount path from `/mnt/share` to `/mnt/stateful_partition/share` for COS compatibility.
9. Submitted and completed a successful WT cloud smoke run.

## New cloud scripts

- `scripts/gcp/gcloud_cli_docker.sh`
- `scripts/gcp/auth_login_via_docker.sh`
- `scripts/gcp/render_wt_smoke_batch_job.py`
- `scripts/gcp/run_wt_smoke_cloud.sh`
- `scripts/gcp/docker/Dockerfile.wt_smoke`

## Why this flow

You can run Google Cloud from this machine without installing `gcloud` locally by using the official Google Cloud CLI Docker image.

## One-time login

Run from the consolidated folder:

```bash
cd "/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated"
bash scripts/gcp/auth_login_via_docker.sh
```

This opens interactive browser/device login and stores credentials in `~/.config/gcloud`.

## WT smoke test (cloud)

Set environment variables:

```bash
export PROJECT_ID="YOUR_PROJECT_ID"
export BUCKET_NAME="YOUR_UNIQUE_BUCKET_NAME"
export REGION="us-central1"
```

Use the project ID string for `PROJECT_ID` (for example `project-2f27e3aa-de6c-4f65-ad4`), not the numeric project number.

Run:

```bash
bash scripts/gcp/run_wt_smoke_cloud.sh
```

GPU smoke (optional, faster than CPU smoke when quota allows):

```bash
export GPU_TYPE="nvidia-tesla-t4"
export GPU_COUNT="1"
export INSTALL_GPU_DRIVERS="true"
export OPENMM_PLATFORM="CUDA"
export MACHINE_TYPE="n1-standard-4"
export ALLOWED_LOCATIONS="us-west2-b,us-west2-c"

bash scripts/gcp/run_wt_smoke_cloud.sh
```

What it does:

1. Sets active project.
2. Enables APIs (`batch`, `compute`, `logging`, `storage`, `artifactregistry`).
3. Creates Artifact Registry repo if missing.
4. Creates bucket if missing.
5. Builds linux/amd64 WT smoke image and pushes it.
6. Renders `manifests/batch_wt_smoke_job.json`.
7. Submits a one-task Batch job.

When `GPU_TYPE` is set, the script renders a GPU Batch config and switches OpenMM to `CUDA`.

An example config is already generated:

- `manifests/batch_wt_smoke_job.example.json`

## Monitor

```bash
scripts/gcp/gcloud_cli_docker.sh batch jobs list --location "$REGION"
scripts/gcp/gcloud_cli_docker.sh batch jobs describe JOB_NAME --location "$REGION"
scripts/gcp/gcloud_cli_docker.sh batch tasks list --job JOB_NAME --location "$REGION"
```

## If the job succeeds

Outputs are written to the bucket under:

- `gs://$BUCKET_NAME/outputs/wt_smoke/`

## Verified successful run

- Job name: `wt-smoke-20260428-060632`
- State: `SUCCEEDED`
- Runtime: `2915.97s` (~48.6 min)
- Output prefix:
  - `gs://zfp36l2-md-20260428-12676/outputs/wt_smoke/`

## References

- Google Cloud CLI Docker image:
  - https://cloud.google.com/sdk/docs/downloads-docker
- Batch basic jobs and task env vars:
  - https://cloud.google.com/batch/docs/create-run-basic-job
- Batch GCS volume mount (`remotePath`):
  - https://cloud.google.com/batch/docs/create-run-example-job
