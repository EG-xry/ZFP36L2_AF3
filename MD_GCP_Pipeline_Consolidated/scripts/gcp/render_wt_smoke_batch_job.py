#!/usr/bin/env python3
"""Render a Google Cloud Batch job JSON for WT smoke test."""
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image-uri", required=True, help="Container image URI in Artifact Registry.")
    parser.add_argument(
        "--bucket-name",
        required=True,
        help="Cloud Storage bucket name only (no gs:// prefix) mounted into Batch job.",
    )
    parser.add_argument("--machine-type", default="e2-standard-4", help="Compute Engine machine type.")
    parser.add_argument(
        "--service-account",
        default=None,
        help="Service account email used by Batch VM tasks.",
    )
    parser.add_argument("--gpu-type", default="", help="GPU accelerator type, e.g. nvidia-tesla-t4.")
    parser.add_argument("--gpu-count", type=int, default=1, help="GPU count per VM when gpu-type is set.")
    parser.add_argument(
        "--install-gpu-drivers",
        default="true",
        help="Whether Batch auto-installs GPU drivers when GPU is requested (true/false).",
    )
    parser.add_argument(
        "--allowed-locations",
        default="",
        help="Comma-separated region/zone allow-list (for example us-west2,us-west2-b).",
    )
    parser.add_argument("--cpu-milli", type=int, default=2000, help="vCPU in milliCPU units per task.")
    parser.add_argument("--memory-mib", type=int, default=8192, help="Memory MiB per task.")
    parser.add_argument("--max-run-duration", default="14400s", help="Task max runtime, e.g. 3600s.")
    parser.add_argument("--platform", default="CPU", help="OpenMM platform argument (CPU or CUDA).")
    parser.add_argument(
        "--out-json",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "manifests" / "batch_wt_smoke_job.json",
        help="Output JSON path.",
    )

    parser.add_argument("--stage1-ns", type=float, default=0.01)
    parser.add_argument("--stage2-ns", type=float, default=0.02)
    parser.add_argument("--production-ns", type=float, default=0.05)
    args = parser.parse_args()

    if args.gpu_type and args.gpu_count < 1:
        raise SystemExit("gpu-count must be >= 1 when gpu-type is set")

    cmd = (
        "set -euo pipefail; "
        "mkdir -p /mnt/stateful_partition/share/outputs/wt_smoke; "
        "python /workspace/templates/zfp36l2_ubuntu_prep_run/run_minimal_md.py "
        "--pdbs /workspace/templates/zfp36l2_ubuntu_prep_run/cleaned_pdbs/"
        "fold_zfp36l2_wt_nanomer_9_wt_2zn_model_0_fixed.pdb "
        "--output-root /mnt/stateful_partition/share/outputs/wt_smoke "
        f"--platform {args.platform} "
        f"--stage1-ns {args.stage1_ns} "
        f"--stage2-ns {args.stage2_ns} "
        f"--production-ns {args.production_ns}"
    )

    instance_policy = {"machineType": args.machine_type}
    instance_entry = {"policy": instance_policy}
    if args.gpu_type:
        instance_policy["accelerators"] = [{"type": args.gpu_type, "count": args.gpu_count}]
        instance_entry["installGpuDrivers"] = parse_bool(args.install_gpu_drivers)

    allocation_policy = {"instances": [instance_entry]}
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
                                "commands": ["-lc", cmd],
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
                "taskCount": 1,
                "parallelism": 1,
            }
        ],
        "allocationPolicy": allocation_policy,
        "labels": {
            "project": "zfp36l2",
            "pipeline": "wt-smoke",
        },
        "logsPolicy": {
            "destination": "CLOUD_LOGGING",
        },
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(config, indent=2), encoding="utf-8")
    print(f"Wrote WT smoke Batch config: {args.out_json}")


if __name__ == "__main__":
    main()
