# Local Smoke Tests (Mac)

## Current local compatibility snapshot

From `scripts/utils/check_local_md_stack.sh`:

- Binaries present: `python3`, `docker`, `mpirun`, `plumed`
- Missing binary: `gmx`
- Base `python3` is missing MD packages (`openmm`, `openmmplumed`, `mpi4py`, `mdtraj`, `pandas`, `Bio`)
- `micromamba` env `omm` contains all key Python MD packages

## Quick local checks

```bash
bash scripts/utils/check_local_md_stack.sh
python3 -m py_compile scripts/utils/build_variant_manifest.py scripts/gcp/stage_variant_inputs.py scripts/gcp/build_task_environments.py
```

## Workaround path on Mac

Use Linux/amd64 container runtime for reproducibility:

- Reference files:
  - `scripts/afmi_workspace/Dockerfile`
  - `scripts/afmi_workspace/docker-compose.yml`
  - `scripts/afmi_workspace/run-docker.sh`

Why: prior notes indicate plugin/wrapper behavior is more reliable in Linux-style runtime.

## Ready-made smoke-test dataset

- `templates/tdp43_reference_assets/` contains a complete historical prep snapshot.
- Use it to validate that container/runtime plumbing works before wiring ZFP36L2-specific prep.

## Practical recommendation

- Treat local Mac runs as smoke tests only.
- Treat GCP Linux runs as the source of truth for full 12-variant execution.

## Local limitations observed in prior logs

- SWIG cross-module force conversion error occurred in one run history.
- PLUMED-without-MPI error occurred when `WALKERS_MPI` was requested.

These need to be validated in the chosen cloud runtime before large-scale execution.
