# Source Map

This folder was assembled from the following local sources.

## Main source directories

- `/Users/eric/Desktop/ZFP36L2 AF3/AlphaFold`
  - Ubuntu-transferred MD workspace now moved into this project
- `/Users/eric/Desktop/ZFP36L2 AF3/attachments`
  - project paper PDF and older PDB attachments
- `/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main`
  - AF-MI workflow scripts, Docker artifacts, wrapper-fix notes
- `/Users/eric/Desktop/Joshua_Project/AlphaFold-IDP`
  - prep scripts and legacy reference material
- `/Users/eric/Desktop/openmm`
  - prior local run outputs/logs
- `/Users/eric/Desktop/ZFP36L2 AF3/md_simulation`
  - local MD output traces for this project
- `/Users/eric/Desktop/ZFP36L2 AF3/results/converted_pdb`
  - run1/model0 variant PDB sources

## Consolidation structure

- `scripts/` holds executable helpers and migrated workflow scripts.
- `templates/` holds reusable prep/template files.
  - Includes `templates/tdp43_reference_assets/` as a known-good historical baseline.
- `inputs/` holds curated variant/reference inputs.
  - Includes `inputs/reference_paper/Joshua_Project_Silvia.pdf`.
- `evidence/` holds historical logs/outputs for troubleshooting.
- `manifests/` holds run manifests and Batch task mappings.
- `runs/variants/` holds staged per-variant inputs for execution.
