# Next Steps

## Immediate

1. Launch the first full 10-variant cloud run using the new orchestrator:
   - `bash scripts/gcp/run_variants_cloud.sh`
2. Monitor Batch states until all tasks are terminal:
   - `scripts/gcp/gcloud_cli_docker.sh batch jobs describe JOB_NAME --location REGION`
   - `scripts/gcp/gcloud_cli_docker.sh batch tasks list --job JOB_NAME --location REGION`
3. Confirm each variant wrote outputs under:
   - `gs://BUCKET_NAME/outputs/variants/<VARIANT>/`

## Runtime knobs for first 10-variant pass

1. Fast pipeline-validation run (recommended first):
   - `STAGE1_NS=0.01`
   - `STAGE2_NS=0.02`
   - `PRODUCTION_NS=0.05`
2. Closer-to-paper stress test:
   - increase durations and `MAX_RUN_DURATION`
   - reduce `PARALLELISM` if quota or cost becomes limiting

## After first 10-variant run

1. Pull or inspect outputs for all variants.
2. Generate an updated QC table and mark pass/fail per variant.
3. Decide whether to run longer production trajectories for selected variants.

## Longer-term alignment with paper protocol

1. Map from current OpenMM minimal pipeline to paper settings:
   - temperature 310 K
   - longer production (paper reports 1000 ns, 3 replicates)
2. Decide whether to keep OpenMM for throughput or add a dedicated GROMACS branch for strict replication.
