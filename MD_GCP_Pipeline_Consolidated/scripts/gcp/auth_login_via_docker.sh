#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
GCLOUD_WRAPPER="$ROOT_DIR/scripts/gcp/gcloud_cli_docker.sh"

echo "Starting interactive gcloud login in a Dockerized gcloud CLI..."
"$GCLOUD_WRAPPER" auth login

echo

echo "Active accounts:"
"$GCLOUD_WRAPPER" auth list

echo

echo "If you also want application default credentials for some tools:"
echo "  $GCLOUD_WRAPPER auth application-default login"
