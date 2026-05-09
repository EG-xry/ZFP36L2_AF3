#!/usr/bin/env bash
set -euo pipefail

# Runs a WT-only cloud smoke test end-to-end using Google Cloud Batch.
# Requires Docker Desktop running and user auth via gcloud container wrapper.

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
GCLOUD_WRAPPER="$ROOT_DIR/scripts/gcp/gcloud_cli_docker.sh"

gcloudw() {
  GCLOUD_WORKSPACE_DIR="$ROOT_DIR" "$GCLOUD_WRAPPER" "$@"
}

PROJECT_ID="${PROJECT_ID:-}"
REGION="${REGION:-us-central1}"
AR_REPO="${AR_REPO:-zfp36l2-md}"
BUCKET_NAME="${BUCKET_NAME:-}"
IMAGE_NAME="${IMAGE_NAME:-zfp36l2-wt-smoke}"
TAG="${TAG:-$(date +%Y%m%d-%H%M%S)}"
JOB_NAME="${JOB_NAME:-wt-smoke-${TAG}}"
JOB_SERVICE_ACCOUNT="${JOB_SERVICE_ACCOUNT:-}"
MACHINE_TYPE="${MACHINE_TYPE:-}"
CPU_MILLI="${CPU_MILLI:-}"
MEMORY_MIB="${MEMORY_MIB:-}"
GPU_TYPE="${GPU_TYPE:-}"
GPU_COUNT="${GPU_COUNT:-1}"
INSTALL_GPU_DRIVERS="${INSTALL_GPU_DRIVERS:-true}"
ALLOWED_LOCATIONS="${ALLOWED_LOCATIONS:-}"
OPENMM_PLATFORM="${OPENMM_PLATFORM:-}"

if [[ -n "$GPU_TYPE" ]]; then
  MACHINE_TYPE="${MACHINE_TYPE:-n1-standard-4}"
  CPU_MILLI="${CPU_MILLI:-4000}"
  MEMORY_MIB="${MEMORY_MIB:-15360}"
  OPENMM_PLATFORM="${OPENMM_PLATFORM:-CUDA}"
else
  MACHINE_TYPE="${MACHINE_TYPE:-e2-standard-4}"
  CPU_MILLI="${CPU_MILLI:-2000}"
  MEMORY_MIB="${MEMORY_MIB:-8192}"
  OPENMM_PLATFORM="${OPENMM_PLATFORM:-CPU}"
fi

if [[ "$OPENMM_PLATFORM" == "CUDA" && -z "$GPU_TYPE" ]]; then
  echo "ERROR: OPENMM_PLATFORM=CUDA requires GPU_TYPE (for example nvidia-tesla-t4)."
  exit 1
fi

if [[ "$OPENMM_PLATFORM" != "CUDA" && -n "$GPU_TYPE" ]]; then
  echo "WARNING: GPU_TYPE is set but OPENMM_PLATFORM='$OPENMM_PLATFORM'."
  echo "         OpenMM may not use the GPU unless platform is CUDA."
fi

if [[ -z "$PROJECT_ID" || -z "$BUCKET_NAME" ]]; then
  echo "ERROR: set PROJECT_ID and BUCKET_NAME before running."
  echo "Example:"
  echo "  export PROJECT_ID='my-gcp-project'"
  echo "  export BUCKET_NAME='my-zfp36l2-batch-bucket'"
  echo "  export REGION='us-central1'"
  echo "  bash scripts/gcp/run_wt_smoke_cloud.sh"
  exit 1
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "ERROR: docker is not installed."
  exit 1
fi

if ! docker info >/dev/null 2>&1; then
  echo "ERROR: docker daemon is not running. Start Docker Desktop first."
  exit 1
fi

if [[ "$PROJECT_ID" =~ ^[0-9]+$ ]]; then
  echo "PROJECT_ID looks like a project number. Resolving to project ID string..."
  RESOLVED_PROJECT_ID="$(gcloudw projects describe "$PROJECT_ID" --format='value(projectId)')"
  if [[ -z "$RESOLVED_PROJECT_ID" ]]; then
    echo "ERROR: Could not resolve project number '$PROJECT_ID' to a project ID."
    exit 1
  fi
  echo "Resolved project number '$PROJECT_ID' -> project ID '$RESOLVED_PROJECT_ID'"
  PROJECT_ID="$RESOLVED_PROJECT_ID"
fi

IMAGE_URI="${REGION}-docker.pkg.dev/${PROJECT_ID}/${AR_REPO}/${IMAGE_NAME}:${TAG}"
JOB_CONFIG="$ROOT_DIR/manifests/batch_wt_smoke_job.json"
JOB_CONFIG_CONTAINER="/workspace/manifests/batch_wt_smoke_job.json"

PROJECT_NUMBER="$(gcloudw projects describe "$PROJECT_ID" --format='value(projectNumber)')"
DEFAULT_COMPUTE_SA="${PROJECT_NUMBER}-compute@developer.gserviceaccount.com"
JOB_SERVICE_ACCOUNT="${JOB_SERVICE_ACCOUNT:-$DEFAULT_COMPUTE_SA}"

echo "[1/8] Checking gcloud wrapper"
gcloudw version >/dev/null

echo "[2/8] Setting active project"
gcloudw config set project "$PROJECT_ID" >/dev/null

echo "[3/8] Enabling required APIs"
gcloudw services enable \
  batch.googleapis.com \
  compute.googleapis.com \
  logging.googleapis.com \
  storage.googleapis.com \
  artifactregistry.googleapis.com

echo "[3.5/8] Ensuring Batch Agent Reporter role for job service account"
HAS_AGENT_REPORTER_ROLE="$(
  gcloudw projects get-iam-policy "$PROJECT_ID" \
    --flatten='bindings[].members' \
    --filter="bindings.role:roles/batch.agentReporter AND bindings.members:serviceAccount:${JOB_SERVICE_ACCOUNT}" \
    --format='value(bindings.role)' || true
)"
if [[ -z "$HAS_AGENT_REPORTER_ROLE" ]]; then
  gcloudw projects add-iam-policy-binding "$PROJECT_ID" \
    --member="serviceAccount:${JOB_SERVICE_ACCOUNT}" \
    --role="roles/batch.agentReporter" \
    --quiet >/dev/null
  echo "Granted roles/batch.agentReporter to ${JOB_SERVICE_ACCOUNT}"
else
  echo "roles/batch.agentReporter already granted to ${JOB_SERVICE_ACCOUNT}"
fi

echo "[3.6/8] Ensuring Logging Writer role for job service account"
HAS_LOG_WRITER_ROLE="$(
  gcloudw projects get-iam-policy "$PROJECT_ID" \
    --flatten='bindings[].members' \
    --filter="bindings.role:roles/logging.logWriter AND bindings.members:serviceAccount:${JOB_SERVICE_ACCOUNT}" \
    --format='value(bindings.role)' || true
)"
if [[ -z "$HAS_LOG_WRITER_ROLE" ]]; then
  gcloudw projects add-iam-policy-binding "$PROJECT_ID" \
    --member="serviceAccount:${JOB_SERVICE_ACCOUNT}" \
    --role="roles/logging.logWriter" \
    --quiet >/dev/null
  echo "Granted roles/logging.logWriter to ${JOB_SERVICE_ACCOUNT}"
else
  echo "roles/logging.logWriter already granted to ${JOB_SERVICE_ACCOUNT}"
fi

echo "[4/8] Ensuring Artifact Registry repository exists"
if ! gcloudw artifacts repositories describe "$AR_REPO" --location "$REGION" >/dev/null 2>&1; then
  gcloudw artifacts repositories create "$AR_REPO" \
    --repository-format=docker \
    --location="$REGION" \
    --description="ZFP36L2 MD images"
fi

echo "[4.5/8] Ensuring Artifact Registry reader access for job service account"
REPO_POLICY_JSON="$(
  gcloudw artifacts repositories get-iam-policy "$AR_REPO" \
    --location "$REGION" \
    --format=json || true
)"
HAS_AR_READER_ROLE="$(
  POLICY_JSON="$REPO_POLICY_JSON" python3 - "$JOB_SERVICE_ACCOUNT" <<'PY'
import json
import os
import sys

sa = f"serviceAccount:{sys.argv[1]}"
try:
    policy = json.loads(os.environ.get("POLICY_JSON", ""))
except Exception:
    print("0")
    raise SystemExit(0)

for binding in policy.get("bindings", []):
    if binding.get("role") == "roles/artifactregistry.reader" and sa in binding.get("members", []):
        print("1")
        raise SystemExit(0)
print("0")
PY
)"
if [[ "$HAS_AR_READER_ROLE" != "1" ]]; then
  gcloudw artifacts repositories add-iam-policy-binding "$AR_REPO" \
    --location "$REGION" \
    --member="serviceAccount:${JOB_SERVICE_ACCOUNT}" \
    --role="roles/artifactregistry.reader" >/dev/null
  echo "Granted roles/artifactregistry.reader on ${AR_REPO} to ${JOB_SERVICE_ACCOUNT}"
else
  echo "roles/artifactregistry.reader already granted on ${AR_REPO}"
fi

echo "[5/8] Ensuring Cloud Storage bucket exists"
if ! gcloudw storage buckets describe "gs://$BUCKET_NAME" >/dev/null 2>&1; then
  gcloudw storage buckets create "gs://$BUCKET_NAME" --location "$REGION"
fi

echo "[5.5/8] Ensuring bucket write access for job service account"
BUCKET_POLICY_JSON="$(gcloudw storage buckets get-iam-policy "gs://$BUCKET_NAME" --format=json || true)"
HAS_BUCKET_WRITE_ROLE="$(
  POLICY_JSON="$BUCKET_POLICY_JSON" python3 - "$JOB_SERVICE_ACCOUNT" <<'PY'
import json
import os
import sys

sa = f"serviceAccount:{sys.argv[1]}"
try:
    policy = json.loads(os.environ.get("POLICY_JSON", ""))
except Exception:
    print("0")
    raise SystemExit(0)

for binding in policy.get("bindings", []):
    if binding.get("role") == "roles/storage.objectAdmin" and sa in binding.get("members", []):
        print("1")
        raise SystemExit(0)
print("0")
PY
)"
if [[ "$HAS_BUCKET_WRITE_ROLE" != "1" ]]; then
  gcloudw storage buckets add-iam-policy-binding "gs://$BUCKET_NAME" \
    --member="serviceAccount:${JOB_SERVICE_ACCOUNT}" \
    --role="roles/storage.objectAdmin" >/dev/null
  echo "Granted roles/storage.objectAdmin on gs://$BUCKET_NAME to ${JOB_SERVICE_ACCOUNT}"
else
  echo "roles/storage.objectAdmin already granted on gs://$BUCKET_NAME"
fi

echo "[6/8] Authenticating docker to Artifact Registry"
ACCESS_TOKEN="$(gcloudw auth print-access-token)"
echo "$ACCESS_TOKEN" | docker login -u oauth2accesstoken --password-stdin "https://${REGION}-docker.pkg.dev" >/dev/null

echo "[7/8] Building and pushing WT smoke image: $IMAGE_URI"
cd "$ROOT_DIR"
docker buildx build \
  --platform linux/amd64 \
  -f scripts/gcp/docker/Dockerfile.wt_smoke \
  -t "$IMAGE_URI" \
  --push \
  .

echo "[8/8] Rendering Batch config + submitting job"
python3 "$ROOT_DIR/scripts/gcp/render_wt_smoke_batch_job.py" \
  --image-uri "$IMAGE_URI" \
  --bucket-name "$BUCKET_NAME" \
  --machine-type "$MACHINE_TYPE" \
  --service-account "$JOB_SERVICE_ACCOUNT" \
  --cpu-milli "$CPU_MILLI" \
  --memory-mib "$MEMORY_MIB" \
  --gpu-type "$GPU_TYPE" \
  --gpu-count "$GPU_COUNT" \
  --install-gpu-drivers "$INSTALL_GPU_DRIVERS" \
  --allowed-locations "$ALLOWED_LOCATIONS" \
  --platform "$OPENMM_PLATFORM" \
  --out-json "$JOB_CONFIG"

gcloudw batch jobs submit "$JOB_NAME" \
  --location "$REGION" \
  --config "$JOB_CONFIG_CONTAINER"

echo

echo "Submitted Batch job: $JOB_NAME"
echo "Region: $REGION"
echo "Image: $IMAGE_URI"
echo "Config: $JOB_CONFIG"
echo "Machine type: $MACHINE_TYPE"
echo "Task resources: cpuMilli=$CPU_MILLI memoryMiB=$MEMORY_MIB"
if [[ -n "$GPU_TYPE" ]]; then
  echo "GPU per instance: type=$GPU_TYPE count=$GPU_COUNT installDrivers=$INSTALL_GPU_DRIVERS"
  if [[ -n "$ALLOWED_LOCATIONS" ]]; then
    echo "Allowed locations: $ALLOWED_LOCATIONS"
  fi
fi
echo "OpenMM platform: $OPENMM_PLATFORM"
echo

echo "Monitor commands:"
echo "  $GCLOUD_WRAPPER batch jobs describe '$JOB_NAME' --location '$REGION'"
echo "  $GCLOUD_WRAPPER batch tasks list --job '$JOB_NAME' --location '$REGION'"
echo "  $GCLOUD_WRAPPER logging read 'resource.type=batch_task AND labels.job_uid:*' --limit 50"
