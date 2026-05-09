# Project Status (2026-04-28)

## Latest update (CPU fallback relaunch)

- Paused GPU route due quota/capacity uncertainty:
  - `GPUS_ALL_REGIONS = 0` at project scope
  - repeated `ZONE_RESOURCE_POOL_EXHAUSTED` across tested GPU zones
- Deleted queued GPU variants job:
  - `variants-md-20260428-092316` (us-central1, T4)
- Relaunched full 10-variant run on CPU in `us-west2`:
  - Job: `variants-md-20260428-092928`
  - State after launch: `RUNNING`
  - Parallelism: `6`
  - Task count: `10`
  - Machine type: `e2-standard-4`
  - Per-task resources: `cpuMilli=2000`, `memoryMib=8192`
  - Stage durations: `0.01 ns + 0.02 ns + 0.05 ns`

## Latest update (failure fixes + targeted rerun)

- Diagnosed first 10-variant CPU run failure:
  - `variants-md-20260428-092928` ended `FAILED`
  - `8/10` succeeded; failures were `S185N` and `C208A`
  - Failure signatures:
    - `S185N`: `openmm.OpenMMException: Particle coordinate is NaN`
    - `C208A`: `No template found for residue ... (ZN)` due malformed Zn element alignment in source PDB line formatting
- Implemented runtime fixes:
  - `clean_af3_structures.py`
    - now rewrites all Zn atoms in canonical PDB formatting
    - assigns unique Zn residue numbers
    - avoids element-column drift from raw AF3 lines
  - `run_minimal_md.py`
    - added production-stage NaN fallback:
      - retry with `1.0 fs` timestep
      - add light heavy-atom positional restraints
- Submitted targeted rerun with patched image for failed variants only:
  - Job: `variants-rerun-failed-20260428-143042`
  - Variants: `S185N`, `C208A`
  - Final outcome:
    - `S185N`: `SUCCEEDED`
    - `C208A`: `SUCCEEDED` (used NaN fallback path during production)
- Added cloud-output analysis script:
  - `scripts/utils/analyze_cloud_minimal_md.py`
  - Current analysis snapshot file:
    - `manifests/cloud_variant_analysis.csv`
  - Completed variants in analysis at last sync: `10/10`

## Objective

Run one MD simulation per selected ZFP36L2 variant (10 variants) using the consolidated OpenMM pipeline on Google Cloud Batch.

Selected variants:
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

Explicit exclusions:
- R160E
- H212A

## Completed

- Consolidated messy source folders into this workspace (`attachments`, `AlphaFold`, Ubuntu-transferred scripts/logs).
- Built 10-variant manifest and staged inputs:
  - `manifests/variant_run_manifest.csv`
  - `runs/variants/<VARIANT>/input_af.pdb`
  - `manifests/batch_task_environments.json`
- Implemented cloud auth + submit flow via Dockerized gcloud wrapper (`scripts/gcp/gcloud_cli_docker.sh`).
- Patched path/permission issues in WT cloud smoke submission.
- WT cloud smoke test succeeded:
  - Job: `wt-smoke-20260428-060632`
  - Final state: `SUCCEEDED`
  - Runtime: `2915.97s` (~48.6 min)
  - Completion timestamp: `2026-04-28T13:57:17Z`
  - Output path: `gs://zfp36l2-md-20260428-12676/outputs/wt_smoke/`
- Added full multi-variant cloud launch scripts:
  - `scripts/gcp/render_variant_batch_job.py`
  - `scripts/gcp/run_variants_cloud.sh`
- Improved multi-variant task robustness:
  - run MD to local VM scratch first
  - copy finalized outputs to GCS mount at end
  - avoid direct long trajectory writes to gcsfuse during simulation
- Added per-task AF3 PDB cleaning before MD (`clean_af3_structures.py`) to avoid RNA template mismatches on raw variant PDBs.
- Added GPU-capable Batch config fields to multi-variant renderer:
  - `--gpu-type`, `--gpu-count`, `--install-gpu-drivers`
  - `--machine-type`, `--cpu-milli`, `--memory-mib`, `--allowed-locations`
- Added GPU runtime hardening for Batch JSON:
  - moved `installGpuDrivers` to `allocationPolicy.instances[*]` level
  - kept accelerator definition under `policy.accelerators[]`
- Updated launch scripts for GPU-first behavior:
  - `run_variants_cloud.sh` and `run_wt_smoke_cloud.sh` now auto-switch defaults when `GPU_TYPE` is set:
    - `OPENMM_PLATFORM=CUDA`
    - `MACHINE_TYPE=n1-standard-4`
    - `PARALLELISM=1` (variants)
- Terminated in-flight CPU variants job to avoid spending time/cost on non-GPU runs:
  - Deleted: `variants-md-20260428-074756`
- Verified current accelerator visibility/quotas in `us-west2`:
  - available accelerator type includes `nvidia-tesla-t4` in `us-west2-b` / `us-west2-c`
  - regional quota includes `NVIDIA_T4_GPUS=1`
- Submitted full 10-variant GPU Batch job:
  - Job: `variants-md-20260428-080657`
  - Config highlights: `OPENMM_PLATFORM=CUDA`, `GPU_TYPE=nvidia-tesla-t4`, `taskCount=10`, `parallelism=1`
  - Initial state after submission: `SCHEDULED` with first task `PENDING`
- Diagnosed root cause for repeated `SCHEDULED` jobs:
  - Batch MIG creation repeatedly fails with `ZONE_RESOURCE_POOL_EXHAUSTED` (capacity stockout), not script/auth issues.
  - Verified via `compute instance-groups managed list-errors` on these zones:
    - `us-west2-b`, `us-west2-c`
    - `us-west1-b`, `us-west1-a`
    - `us-central1-f` (latest cross-region attempt)
- Tried multiple GPU/region fallbacks to unblock:
  - `T4` in `us-west2`, `us-west1`, and `us-central1`
  - `P4` in `us-west2-c`
  - `P100` in `us-west1`
  - `L4` in `us-west1-a` (`g2-standard-4`)
- Additional hard quota blocker identified for `L4` path:
  - `QUOTA_EXCEEDED: GPUS_ALL_REGIONS limit 0.0`
- Latest queued fallback job:
  - Job: `variants-md-20260428-091833`
  - Region: `us-central1`
  - GPU: `nvidia-tesla-t4`
  - Allowed zones: `us-central1-a`, `us-central1-c`, `us-central1-f`

## Current blockers

1. GPU capacity shortages (`ZONE_RESOURCE_POOL_EXHAUSTED`) are preventing VM allocation across tested zones.
2. `L4` fallback currently blocked by global quota `GPUS_ALL_REGIONS = 0.0`.
3. Need one successfully running GPU-backed variant job (10/10) before QC.
4. Current cloud flow is a practical minimal-OpenMM pipeline, not yet a strict replication of the paper's GROMACS protocol:
   - paper: GROMACS 2020.3, 310 K, 1000 ns production, 3 replicates
   - current: OpenMM minimal pipeline, single run per variant (durations configurable)
5. Post-run QC still needs to be executed for all 10 variants after first full cloud run.

## Evidence summary

- Prior failure signatures preserved under `evidence/afmi_logs/`.
- Prior local success evidence preserved under `evidence/local_runs/`.
- Cloud smoke success evidence is now available in:
  - Batch job metadata (`wt-smoke-20260428-060632`)
  - bucket outputs under `outputs/wt_smoke/`

## Definition of done for the first cloud variant set

- One successful Batch task for each selected variant (10/10), with:
  - input traceability (`runs/variants/<VARIANT>/metadata.json`)
  - output files per variant in bucket (`outputs/variants/<VARIANT>/...`)
  - no fatal runtime errors in Batch task logs
