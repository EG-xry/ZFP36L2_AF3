#!/bin/bash
# ============================================================================
# Extend completed 50 ns runs to 400 ns using production.chk (auto --resume).
# Each job uses the 24 h wall limit; resubmit the same command if checkpointed.
# ============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/zfp36l2_md"
TARGET_NS="${1:-400}"
SUBMIT_SCRIPT="${PROJECT_DIR}/submit_variant_expanse.sh"
if [[ ! -f "$SUBMIT_SCRIPT" ]]; then
    SUBMIT_SCRIPT="${PROJECT_DIR}/expanse/submit_variant_expanse.sh"
fi

# Variants that finished 50 ns successfully
VARIANTS=(C206F G204D D210E T187S)

echo "Extending ${#VARIANTS[@]} variants to ${TARGET_NS} ns (resume from checkpoint)..."
for VARIANT in "${VARIANTS[@]}"; do
    OUT="${PROJECT_DIR}/outputs/${VARIANT}"
    if [[ ! -f "${OUT}/production.chk" || ! -f "${OUT}/equilibrated.chk" ]]; then
        echo "  SKIP ${VARIANT}: missing checkpoint(s)"
        continue
    fi
    echo -n "  ${VARIANT}... "
    JOB_ID=$(sbatch --parsable -A sdp166 "$SUBMIT_SCRIPT" "$VARIANT" "$TARGET_NS")
    echo "job ${JOB_ID}"
done

echo ""
echo "Monitor: squeue -u \$USER"
echo "Each job appends to existing production.dcd from ~50 ns → ${TARGET_NS} ns."
