#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
WORKSPACE_DIR="$(cd "$ROOT_DIR/.." && pwd)"

GCLOUD_WRAPPER="${GCLOUD_WRAPPER:-$ROOT_DIR/scripts/gcp/gcloud_cli_docker.sh}"
PROJECT_ID="${PROJECT_ID:-project-2f27e3aa-de6c-4f65-ad4}"
REGION="${REGION:-us-west2}"
JOB_NAME="${JOB_NAME:-variants-md-12hmax-20260428-190901}"
BUCKET_NAME="${BUCKET_NAME:-zfp36l2-md-20260428-12676}"
WAIT_FOR_JOB="${WAIT_FOR_JOB:-true}"
ALLOW_PARTIAL="${ALLOW_PARTIAL:-false}"
POLL_SECONDS="${POLL_SECONDS:-120}"
PRESENTATION_DIR="${PRESENTATION_DIR:-$ROOT_DIR/presentation/compiled_$(date +%Y%m%d)}"

PYTHON_BIN="${PYTHON_BIN:-$WORKSPACE_DIR/.venv/bin/python}"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="python3"
fi

log() {
  printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S %Z')" "$*"
}

job_state() {
  "$GCLOUD_WRAPPER" batch jobs describe "$JOB_NAME" \
    --location "$REGION" \
    --project "$PROJECT_ID" \
    --format='value(status.state)'
}

job_counts() {
  "$GCLOUD_WRAPPER" batch jobs describe "$JOB_NAME" \
    --location "$REGION" \
    --project "$PROJECT_ID" \
    --format='value(status.taskGroups.group0.counts.SUCCEEDED,status.taskGroups.group0.counts.RUNNING,status.taskGroups.group0.counts.FAILED,status.taskGroups.group0.counts.PENDING)'
}

if [[ "$WAIT_FOR_JOB" == "true" ]]; then
  log "Waiting for Batch job: $JOB_NAME"
  while true; do
    state="$(job_state)"
    counts="$(job_counts || true)"
    log "Job state=$state counts=$counts"

    case "$state" in
      SUCCEEDED)
        break
        ;;
      FAILED|DELETION_IN_PROGRESS)
        if [[ "$ALLOW_PARTIAL" == "true" ]]; then
          log "Job ended in state=$state; continuing due ALLOW_PARTIAL=true"
          break
        fi
        log "Job ended in state=$state; stopping. Set ALLOW_PARTIAL=true to proceed anyway."
        exit 2
        ;;
      QUEUED|SCHEDULED|RUNNING)
        sleep "$POLL_SECONDS"
        ;;
      *)
        log "Unexpected state: $state"
        sleep "$POLL_SECONDS"
        ;;
    esac
  done
else
  log "Skipping wait (WAIT_FOR_JOB=false)."
fi

log "Syncing outputs from gs://$BUCKET_NAME/outputs/variants -> $ROOT_DIR/runs/cloud_outputs/variants"
mkdir -p "$ROOT_DIR/runs/cloud_outputs/variants"
"$GCLOUD_WRAPPER" storage rsync -r \
  "gs://$BUCKET_NAME/outputs/variants" \
  "$ROOT_DIR/runs/cloud_outputs/variants"

log "Running cloud output QC summary"
"$PYTHON_BIN" "$ROOT_DIR/scripts/utils/analyze_cloud_minimal_md.py" \
  --manifest "$ROOT_DIR/manifests/variant_run_manifest.csv" \
  --outputs-root "$ROOT_DIR/runs/cloud_outputs/variants" \
  --out-csv "$ROOT_DIR/manifests/cloud_variant_analysis.csv"

log "Computing structural drift tables"
"$PYTHON_BIN" "$ROOT_DIR/scripts/utils/compute_cloud_structural_drift.py" \
  --manifest "$ROOT_DIR/manifests/variant_run_manifest.csv" \
  --inputs-root "$ROOT_DIR/runs/variants" \
  --outputs-root "$ROOT_DIR/runs/cloud_outputs/variants" \
  --out-ca-csv "$ROOT_DIR/manifests/cloud_variant_structural_drift.csv" \
  --out-rna-csv "$ROOT_DIR/manifests/cloud_variant_structural_drift_rna.csv"

log "Building presentation/report bundle"
"$PYTHON_BIN" "$ROOT_DIR/scripts/utils/build_af3_md_presentation_bundle.py" \
  --output-dir "$PRESENTATION_DIR"

log "Post-run pipeline complete"
log "QC CSV: $ROOT_DIR/manifests/cloud_variant_analysis.csv"
log "Drift CSV: $ROOT_DIR/manifests/cloud_variant_structural_drift.csv"
log "RNA Drift CSV: $ROOT_DIR/manifests/cloud_variant_structural_drift_rna.csv"
log "Bundle: $PRESENTATION_DIR"
