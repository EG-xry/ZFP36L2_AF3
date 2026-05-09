# Script Inventory

## Newly organized utility scripts

- `scripts/utils/build_variant_manifest.py`
  - Builds `manifests/variant_run_manifest.csv` from run1/model0 ZFP36L2 PDB files.
- `scripts/gcp/stage_variant_inputs.py`
  - Creates `runs/variants/<VARIANT>/` with `input_af.pdb` and `metadata.json`.
- `scripts/gcp/build_task_environments.py`
  - Converts manifest rows into `taskEnvironments` entries for Google Cloud Batch.
- `scripts/gcp/gcloud_cli_docker.sh`
  - Runs `gcloud` from Docker (no local `gcloud` install required).
- `scripts/gcp/auth_login_via_docker.sh`
  - Interactive Google auth helper using dockerized `gcloud`.
- `scripts/gcp/render_wt_smoke_batch_job.py`
  - Renders a one-task WT smoke test Batch job JSON.
- `scripts/gcp/run_wt_smoke_cloud.sh`
  - End-to-end cloud smoke test runner (build, push, submit).
- `scripts/gcp/render_variant_batch_job.py`
  - Renders a multi-task Batch job JSON (one task per variant from `taskEnvironments`).
- `scripts/gcp/run_variants_cloud.sh`
  - End-to-end multi-variant runner (refresh manifests, stage/upload inputs, build/push image, submit Batch).
- `scripts/gcp/docker/Dockerfile.wt_smoke`
  - Lightweight linux/amd64 OpenMM runtime image for WT cloud smoke tests.
- `scripts/utils/check_local_md_stack.sh`
  - Checks binaries and Python modules in both `python3` and `micromamba:omm`.
- `scripts/utils/summarize_variant_outputs.py`
  - Generates post-run completeness table (`manifests/variant_output_qc.csv`).
- `scripts/utils/analyze_cloud_minimal_md.py`
  - Computes cloud-run QC metrics (`manifests/cloud_variant_analysis.csv`).
- `scripts/utils/compute_cloud_structural_drift.py`
  - Recomputes protein/RNA start-vs-final structural drift tables from cloud outputs.
- `scripts/utils/build_af3_md_presentation_bundle.py`
  - Builds presentation-ready AF3+MD tables, figures, and summary report bundle.
- `scripts/gcp/postrun_finalize_and_build_reports.sh`
  - One-command post-run pipeline: wait for job, sync GCS outputs, recompute analyses, regenerate figures/report bundle.

## AF-MI workspace scripts (copied as references)

Path: `scripts/afmi_workspace/`

- `Dockerfile`, `docker-compose.yml`, `run-docker.sh`
  - Linux/amd64 container path for Mac compatibility.
- `run_simulation.sh`
  - Example local wrapper for `mpirun -np 2 python simulate.py 7.4 298`.
- `build_plumed_complete.sh`, `start_plumed_build.sh`, `check_plumed_build.sh`
  - PLUMED build workflow scripts.
- `WRAPPER_FIX_REPORT.md`, `VERIFIED_COMMANDS.md`
  - Prior debugging notes for OpenMM-Plumed wrapper and MPI build.

## AF-MI prep templates

Path: `templates/afmi_prep_run/`

- `simulate.py`, `simulate_em.py`, `simulate_quick.py`
  - Main simulation runners.
- `gen_xml_and_constraints.py`
  - Generates `forcefield.xml` and `r1_excl.pkl` from sequence.
- `make_plumed_distmat.py`, `make_plumed_analysis.py`
  - Build PLUMED input files from AF-derived restraints and domain tables.
- `residues.csv`, `ordered_domains.csv`, `disordered_domains.csv`
  - Required tables for setup.

## OpenMM-Plumed wrapper patch templates

Path: `templates/openmm_plumed_patch/`

- `plumedplugin.i`, `add_plumed_force.cpp`, `setup_add_plumed_force.py`
  - Patch route for SWIG/openmmplumed force insertion issues.

## TDP-43 reference assets (for smoke-test baseline only)

Path: `templates/tdp43_reference_assets/`

- Contains a complete historical prep snapshot (`checkpoint`, `plumed.dat`, `forcefield.xml`, `struct*.pdb`, AF-derived arrays/files).
- This is useful to validate runtime behavior quickly.
- Do not treat these files as ZFP36L2 production inputs; they are reference/example assets.

## ZFP36L2 Ubuntu prep assets

Path: `templates/zfp36l2_ubuntu_prep_run/`

- Contains ZFP36L2-specific scripts and setup files from the Ubuntu workspace.
- Includes `run_zfp36l2_nanomer_2zn_md.py`, `run_minimal_md.py`, and cleaned WT/variant PDBs.
- Treat this folder as the most relevant starting point for project-specific MD runs.

## What to treat as archival evidence

- `evidence/afmi_logs/*`
- `evidence/local_runs/*`

These are useful for debugging and parameter recovery but should not be treated as production pipeline scripts.
