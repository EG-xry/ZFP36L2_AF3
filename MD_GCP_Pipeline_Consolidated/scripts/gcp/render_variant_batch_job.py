#!/usr/bin/env python3
"""Render a Google Cloud Batch JSON config for multi-variant minimal MD runs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "y", "on"}:
        return True
    if normalized in {"0", "false", "no", "n", "off"}:
        return False
    raise SystemExit(f"Invalid boolean value: {value!r}")


def load_task_environments(path: Path) -> list[dict]:
    if not path.exists():
        raise SystemExit(f"Task environments JSON not found: {path}")

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SystemExit(f"Invalid JSON in task environments file: {exc}") from exc

    if not isinstance(data, list) or not data:
        raise SystemExit("Task environments JSON must be a non-empty list.")
    if len(data) > 200:
        raise SystemExit("Task environments JSON exceeds Batch limit (200 entries).")

    for idx, entry in enumerate(data):
        if not isinstance(entry, dict):
            raise SystemExit(f"Task environment entry #{idx} is not an object.")
        variables = entry.get("variables")
        if not isinstance(variables, dict):
            raise SystemExit(f"Task environment entry #{idx} is missing object field 'variables'.")
        for key in ("VARIANT", "REL_INPUT_PDB"):
            if not variables.get(key):
                raise SystemExit(f"Task environment entry #{idx} missing required variable: {key}")

    return data


def render_command(stage1_ns: float, stage2_ns: float, production_ns: float, platform: str) -> str:
    # Run on local VM scratch first, then copy final artifacts to the GCS mount.
    # This is more robust than writing trajectories directly to gcsfuse paths.
    lines = [
        "set -euo pipefail",
        'if [[ -z "${VARIANT:-}" || -z "${REL_INPUT_PDB:-}" ]]; then',
        '  echo "Missing VARIANT or REL_INPUT_PDB task env vars."',
        "  exit 2",
        "fi",
        'MOUNT_ROOT="/mnt/stateful_partition/share"',
        'IN_PDB="${MOUNT_ROOT}/${REL_INPUT_PDB}"',
        'if [[ ! -f "$IN_PDB" ]]; then',
        '  echo "Input PDB not found: $IN_PDB"',
        '  ls -la "${MOUNT_ROOT}/variants" || true',
        "  exit 3",
        "fi",
        'SCRATCH_BASE="/tmp/zfp36l2_md/${VARIANT}_${BATCH_TASK_INDEX}"',
        'CLEAN_DIR="${SCRATCH_BASE}/cleaned"',
        'LOCAL_OUT="${SCRATCH_BASE}/output"',
        'mkdir -p "$LOCAL_OUT" "$CLEAN_DIR"',
        "python /workspace/templates/zfp36l2_ubuntu_prep_run/clean_af3_structures.py "
        '--output-dir "$CLEAN_DIR" "$IN_PDB"',
        'CLEAN_PDB="$(find "$CLEAN_DIR" -maxdepth 1 -type f -name "*_cleaned.pdb" | head -n 1)"',
        'if [[ -z "${CLEAN_PDB}" || ! -f "${CLEAN_PDB}" ]]; then',
        '  echo "Failed to create cleaned PDB for ${VARIANT}"',
        "  exit 4",
        "fi",
        "python /workspace/templates/zfp36l2_ubuntu_prep_run/run_minimal_md.py "
        '--pdbs "$CLEAN_PDB" '
        '--output-root "$LOCAL_OUT" '
        f"--platform {platform} "
        f"--stage1-ns {stage1_ns} "
        f"--stage2-ns {stage2_ns} "
        f"--production-ns {production_ns}",
        'REMOTE_OUT="${MOUNT_ROOT}/outputs/variants/${VARIANT}"',
        'mkdir -p "$REMOTE_OUT"',
        'cp -r "$LOCAL_OUT"/. "$REMOTE_OUT"/',
        'printf "VARIANT=%s\\nCLASSIFICATION=%s\\nREL_INPUT_PDB=%s\\nBATCH_TASK_INDEX=%s\\nBATCH_TASK_COUNT=%s\\n" '
        '"${VARIANT}" "${CLASSIFICATION:-unknown}" "${REL_INPUT_PDB}" "${BATCH_TASK_INDEX}" "${BATCH_TASK_COUNT}" '
        '> "${REMOTE_OUT}/task_metadata.env"',
        'rm -rf "$SCRATCH_BASE"',
        'echo "Completed variant ${VARIANT}"',
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image-uri", required=True, help="Artifact Registry image URI.")
    parser.add_argument("--bucket-name", required=True, help="GCS bucket name (without gs://).")
    parser.add_argument(
        "--task-environments-json",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "batch_task_environments.json",
        help="JSON array for taskEnvironments.",
    )
    parser.add_argument("--machine-type", default="e2-standard-4", help="Compute Engine machine type.")
    parser.add_argument("--service-account", default=None, help="Service account email for Batch tasks.")
    parser.add_argument("--cpu-milli", type=int, default=2000, help="CPU milli units per task.")
    parser.add_argument("--memory-mib", type=int, default=8192, help="Memory MiB per task.")
    parser.add_argument("--gpu-type", default="", help="GPU accelerator type (for example nvidia-tesla-t4).")
    parser.add_argument("--gpu-count", type=int, default=1, help="GPU count per task instance.")
    parser.add_argument(
        "--install-gpu-drivers",
        default="true",
        help="Whether Batch should install GPU drivers automatically (true/false).",
    )
    parser.add_argument(
        "--allowed-locations",
        default="",
        help="Comma-separated region/zone allow-list (for example us-west2,us-west2-b).",
    )
    parser.add_argument("--max-run-duration", default="14400s", help="Task timeout, for example 14400s.")
    parser.add_argument("--parallelism", type=int, default=3, help="How many variant tasks to run at once.")
    parser.add_argument("--platform", default="CPU", help="OpenMM platform argument.")
    parser.add_argument("--stage1-ns", type=float, default=0.01, help="NVT duration in ns.")
    parser.add_argument("--stage2-ns", type=float, default=0.02, help="NPT equilibration duration in ns.")
    parser.add_argument("--production-ns", type=float, default=0.05, help="Production duration in ns.")
    parser.add_argument(
        "--out-json",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "batch_variants_job.json",
        help="Output Batch config path.",
    )
    args = parser.parse_args()

    task_envs = load_task_environments(args.task_environments_json)
    task_count = len(task_envs)

    if args.parallelism < 1:
        raise SystemExit("parallelism must be >= 1")
    if args.gpu_type and args.gpu_count < 1:
        raise SystemExit("gpu-count must be >= 1 when gpu-type is set")
    parallelism = min(args.parallelism, task_count)

    command = render_command(args.stage1_ns, args.stage2_ns, args.production_ns, args.platform)

    instance_policy = {
        "machineType": args.machine_type,
    }
    instance_entry = {
        "policy": instance_policy,
    }
    if args.gpu_type:
        instance_policy["accelerators"] = [
            {
                "type": args.gpu_type,
                "count": args.gpu_count,
            }
        ]
        # Preferred placement per Batch API: at instances[*] level.
        instance_entry["installGpuDrivers"] = parse_bool(args.install_gpu_drivers)

    allocation_policy = {
        "instances": [instance_entry]
    }
    if args.allowed_locations.strip():
        allowed_locations = []
        for raw in args.allowed_locations.split(","):
            token = raw.strip()
            if not token:
                continue
            if token.startswith(("regions/", "zones/")):
                allowed_locations.append(token)
            elif token.count("-") >= 2 and token[-1].isalpha():
                allowed_locations.append(f"zones/{token}")
            else:
                allowed_locations.append(f"regions/{token}")
        if allowed_locations:
            allocation_policy["location"] = {"allowedLocations": allowed_locations}
    if args.service_account:
        allocation_policy["serviceAccount"] = {"email": args.service_account}

    config = {
        "taskGroups": [
            {
                "taskSpec": {
                    "runnables": [
                        {
                            "container": {
                                "imageUri": args.image_uri,
                                "entrypoint": "/bin/bash",
                                "commands": ["-lc", command],
                            }
                        }
                    ],
                    "computeResource": {
                        "cpuMilli": args.cpu_milli,
                        "memoryMib": args.memory_mib,
                    },
                    "volumes": [
                        {
                            "gcs": {"remotePath": args.bucket_name},
                            "mountPath": "/mnt/stateful_partition/share",
                        }
                    ],
                    "maxRetryCount": 1,
                    "maxRunDuration": args.max_run_duration,
                },
                "taskEnvironments": task_envs,
                "taskCount": task_count,
                "parallelism": parallelism,
            }
        ],
        "allocationPolicy": allocation_policy,
        "labels": {
            "project": "zfp36l2",
            "pipeline": "variants-minimal-md",
        },
        "logsPolicy": {
            "destination": "CLOUD_LOGGING",
        },
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(config, indent=2), encoding="utf-8")
    print(f"Wrote multi-variant Batch config ({task_count} tasks): {args.out_json}")


if __name__ == "__main__":
    main()
