#!/bin/bash
# ============================================================================
# Submit ALL ZFP36L2 variant MD jobs to Anvil GPU queue
# Each variant gets its own GPU (parallel execution)
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$HOME/zfp36l2_md"

# Production time per variant (ns)
# With 750 GPU hours total and 10 variants:
#   750 / 10 = 75 GPU-hours per variant maximum
#   A100 gives ~200-400 ns/day for this system size (~15k atoms solvated)
#   In 47.5 hours: ~400-790 ns achievable per job
#   Recommendation: 100 ns per variant (safe, ~12h each = 120h total)
#   Or 50 ns per variant (very safe, ~6h each = 60h total)
PRODUCTION_NS="${1:-100}"

echo "============================================================"
echo "Submitting ALL ZFP36L2 variant MD jobs"
echo "Production: ${PRODUCTION_NS} ns per variant"
echo "============================================================"
echo ""

# All 10 variants (9 from your list + nonamer/C208A)
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

# Create logs directory
mkdir -p "${PROJECT_DIR}/logs"
cd "${PROJECT_DIR}"

echo "Variants to submit: ${#VARIANTS[@]}"
echo ""

SUBMITTED=0
for VARIANT in "${VARIANTS[@]}"; do
    echo -n "Submitting ${VARIANT}... "
    JOB_ID=$(sbatch --parsable "${PROJECT_DIR}/scripts/submit_variant.sh" "$VARIANT" "$PRODUCTION_NS")
    echo "Job ID: ${JOB_ID}"
    SUBMITTED=$((SUBMITTED + 1))
done

echo ""
echo "============================================================"
echo "Submitted ${SUBMITTED} jobs."
echo "Monitor with: squeue -u \$USER"
echo "============================================================"
echo ""
echo "Estimated GPU-hours needed:"
echo "  At ~300 ns/day: ${PRODUCTION_NS} ns takes ~$( echo "scale=1; ${PRODUCTION_NS} / 300 * 24" | bc )h per variant"
echo "  Total: ~$( echo "scale=0; ${PRODUCTION_NS} / 300 * 24 * ${SUBMITTED}" | bc )h for all ${SUBMITTED} variants"
echo "  Budget: 750 GPU-hours"
echo ""
echo "To check progress:"
echo "  cat outputs/*/progress.txt"
