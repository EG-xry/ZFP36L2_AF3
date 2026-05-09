#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -eq 0 ]; then
  cat <<'USAGE'
Usage:
  scripts/gcp/gcloud_cli_docker.sh <gcloud args...>

Examples:
  scripts/gcp/gcloud_cli_docker.sh version
  scripts/gcp/gcloud_cli_docker.sh auth login
  scripts/gcp/gcloud_cli_docker.sh config list
USAGE
  exit 1
fi

IMAGE="${GCLOUD_CLI_IMAGE:-gcr.io/google.com/cloudsdktool/google-cloud-cli:stable}"
GCLOUD_CONFIG_DIR="${GCLOUD_CONFIG_DIR:-$HOME/.config/gcloud}"
WORKSPACE_DIR="${GCLOUD_WORKSPACE_DIR:-$(pwd)}"

mkdir -p "$GCLOUD_CONFIG_DIR"

if [ -t 0 ] && [ -t 1 ]; then
  exec docker run --rm -it \
    -v "$GCLOUD_CONFIG_DIR:/root/.config/gcloud" \
    -v "$WORKSPACE_DIR:/workspace" \
    -w /workspace \
    "$IMAGE" \
    gcloud "$@"
else
  exec docker run --rm \
    -v "$GCLOUD_CONFIG_DIR:/root/.config/gcloud" \
    -v "$WORKSPACE_DIR:/workspace" \
    -w /workspace \
    "$IMAGE" \
    gcloud "$@"
fi
