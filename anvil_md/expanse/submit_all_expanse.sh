#!/bin/bash
# ============================================================================
# Submit ALL ZFP36L2 variant MD jobs to SDSC Expanse GPU queue
# Each variant gets its own V100 GPU (parallel execution)
# ============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/zfp36l2_md"
PRODUCTION_NS="${1:-50}"

echo "============================================================"
echo "Submitting ALL ZFP36L2 variant MD jobs to Expanse"
echo "Production: ${PRODUCTION_NS} ns per variant"
echo "============================================================"
echo ""

VARIANTS=(
    "WT"
    "Y154H"
    "R160W"
    "R184C"
    "G204D"
    "C206F"
    "S185N"
    "T187S"
    "D210E"
    "C208A"
)

mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

echo "Variants to submit: ${#VARIANTS[@]}"
echo ""

SUBMITTED=0
for VARIANT in "${VARIANTS[@]}"; do
    echo -n "Submitting ${VARIANT}... "
    JOB_ID=$(sbatch --parsable "${PROJECT_DIR}/scripts/submit_variant_expanse.sh" "$VARIANT" "$PRODUCTION_NS")
    echo "Job ID: ${JOB_ID}"
    SUBMITTED=$((SUBMITTED + 1))
done

echo ""
echo "============================================================"
echo "Submitted ${SUBMITTED} jobs to Expanse."
echo "Monitor with: squeue -u \$USER"
echo "============================================================"
echo ""
echo "V100 performance: ~100-200 ns/day for this system"
echo "  ${PRODUCTION_NS} ns per variant ≈ 6-12h each"
echo "  Total GPU-hours: ~${SUBMITTED} × ~12h = ~$((SUBMITTED * 12))h"
echo "  Budget remaining: 477 GPU-hours"
