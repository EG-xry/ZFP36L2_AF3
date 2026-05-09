# ZFP36L2 MD Pipeline Consolidation (for GCP)

This folder is a cleaned handoff of MD-related work pulled from the messy local/Ubuntu-transferred directories, focused on getting to:

- one MD run per selected variant (10 variants: excludes R160E and H212A)
- reproducible Linux/GCP execution
- clear status of what is already working vs what is still missing

## Start here

1. Read `docs/PROJECT_STATUS.md` for current progress and blockers.
2. Read `docs/GCP_EXECUTION_PLAN.md` for the cloud execution path.
3. Read `docs/SCRIPT_INVENTORY.md` for script ownership and purpose.
4. Read `docs/LOCAL_SMOKE_TESTS.md` for Mac-compatible smoke-test options.
5. Read `docs/COMMANDS_QUICKSTART.md` for executable command sequence.
6. Read `docs/PAPER_NOTES.md` for extracted MD parameter notes from the local ZFP36L2 paper PDF.
7. Read `docs/GCP_WT_SMOKE_RUNBOOK.md` for cloud smoke-test execution.
8. Read `docs/GCP_VARIANT_RUNBOOK.md` for full 10-variant Batch execution.

## Key completed outputs

- Per-variant manifest: `manifests/variant_run_manifest.csv` (10 variants)
- Staged one-run-per-variant folders: `runs/variants/<VARIANT>/input_af.pdb`
- Cloud Batch task mapping helper outputs:
  - `manifests/batch_task_environments.json`
  - `manifests/variant_order.txt`
- Output QC template table:
  - `manifests/variant_output_qc.csv`
- Example WT smoke Batch config:
  - `manifests/batch_wt_smoke_job.example.json`
- WT cloud smoke run completed successfully:
  - job `wt-smoke-20260428-060632`
  - outputs in `gs://zfp36l2-md-20260428-12676/outputs/wt_smoke/`
- Imported local paper PDF:
  - `inputs/reference_paper/Joshua_Project_Silvia.pdf`
- Imported Ubuntu ZFP36L2 prep scripts/logs:
  - `templates/zfp36l2_ubuntu_prep_run/*`
  - `evidence/alphafold_prep_run/*`
- Evidence logs from previous local/Ubuntu runs:
  - `evidence/afmi_logs/*`
  - `evidence/local_runs/*`
- TDP-43 reference prep assets for smoke tests:
  - `templates/tdp43_reference_assets/*`

## Core scripts to use immediately

- `scripts/utils/build_variant_manifest.py`
- `scripts/gcp/stage_variant_inputs.py`
- `scripts/gcp/build_task_environments.py`
- `scripts/gcp/gcloud_cli_docker.sh`
- `scripts/gcp/auth_login_via_docker.sh`
- `scripts/gcp/run_wt_smoke_cloud.sh`
- `scripts/gcp/render_wt_smoke_batch_job.py`
- `scripts/gcp/run_variants_cloud.sh`
- `scripts/gcp/render_variant_batch_job.py`
- `scripts/utils/check_local_md_stack.sh`
- `scripts/utils/summarize_variant_outputs.py`

## Important context

- `attachments/` and `AlphaFold/` now exist under this workspace path and were ingested into this package.
- The local Mac has partial MD tooling; full production runs should be done on Linux (Docker on amd64, GCE, or Batch).
