# GCP Execution Plan

## Goal

Run 10 independent MD jobs (WT + selected variants), one run per variant, with reproducible Linux runtime.

## Recommended rollout

### Phase 1: Single-variant cloud smoke test (required)

1. Build/validate Linux runtime image (OpenMM + PLUMED + MPI + wrapper fix).
2. Run one variant end-to-end (recommend `WT`) on GCE or Batch.
3. Confirm outputs and logs are complete and free of SWIG/MPI errors.

### Phase 2: Full 10-variant test

1. Upload staged inputs (`runs/variants`) and required templates to GCS.
2. Use one Batch task per variant.
3. Collect outputs to structured GCS paths.
4. Run lightweight QC summary per task.
5. Use `scripts/gcp/run_variants_cloud.sh` as the default orchestrator.

## Prepared artifacts for cloud scheduling

- Variant manifest: `manifests/variant_run_manifest.csv` (R160E and H212A excluded by default)
- Task mapping file: `manifests/variant_order.txt`
- Batch `taskEnvironments` array: `manifests/batch_task_environments.json`
- WT smoke config template: `manifests/batch_wt_smoke_job.example.json`
- Full variant config renderer: `scripts/gcp/render_variant_batch_job.py`
- Full variant submit orchestrator: `scripts/gcp/run_variants_cloud.sh`

Generate/update them with:

```bash
python3 scripts/utils/build_variant_manifest.py
python3 scripts/gcp/stage_variant_inputs.py --mode copy
python3 scripts/gcp/build_task_environments.py
```

Default exclusions can be overridden if needed:

```bash
python3 scripts/utils/build_variant_manifest.py --exclude H212A --exclude R160E
```

## Accessing Google Cloud from this machine

- Local `gcloud` binary is not installed.
- Use dockerized gcloud wrapper:
  - `scripts/gcp/gcloud_cli_docker.sh`
  - `scripts/gcp/auth_login_via_docker.sh`

## Suggested Batch task contract

Each task should receive at minimum:

- `VARIANT`
- `REL_INPUT_PDB`
- `TEMPERATURE_K`
- `PH`
- `IONIC_STRENGTH_M`

These are already emitted in `batch_task_environments.json`.

## Known engineering risks to harden before scale-up

1. PLUMED MPI compatibility
   - Ensure PLUMED build supports MPI if `WALKERS_MPI` remains enabled.
2. Wrapper/type compatibility
   - Ensure OpenMM/OpenMMPlumed bindings include the force-add workaround.
3. ZFP36L2-specific prep completeness
   - Confirm per-variant preparation of `plumed.dat`, `forcefield.xml`, exclusions, and any required checkpoints/start states.

## Minimal production folder shape per variant

Each task should write to a unique output directory, e.g.:

- `outputs/<VARIANT>/run.log`
- `outputs/<VARIANT>/stats_*.csv`
- `outputs/<VARIANT>/output_*.dcd`
- `outputs/<VARIANT>/COLVAR*`
- `outputs/<VARIANT>/HILLS*`

## Acceptance checklist

- 10/10 tasks completed successfully.
- No fatal OpenMM/PLUMED/MPI exceptions.
- Expected core output files exist for all variants.
- Metadata ties each output back to source input and run parameters.
