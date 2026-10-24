#!/bin/bash
# ============================================================================
# Fresh 0 → 400 ns runs for specified variants (wipes prior partial outputs).
# Usage: bash submit_fresh_400ns_expanse.sh [VARIANT ...]
# Default: the six variants that failed the 50 ns batch.
# ============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/zfp36l2_md"
TARGET_NS=400
if [[ "${1:-}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    TARGET_NS="$1"
    shift
fi

if [[ $# -gt 0 ]]; then
    VARIANTS=("$@")
else
    VARIANTS=(WT Y154H R160W R184C S185N C208A)
fi

SUBMIT_SCRIPT="${PROJECT_DIR}/submit_variant_expanse.sh"
if [[ ! -f "$SUBMIT_SCRIPT" ]]; then
    SUBMIT_SCRIPT="${PROJECT_DIR}/expanse/submit_variant_expanse.sh"
fi

echo "Fresh ${TARGET_NS} ns runs for: ${VARIANTS[*]}"
for VARIANT in "${VARIANTS[@]}"; do
    OUT="${PROJECT_DIR}/outputs/${VARIANT}"
    if [[ -d "$OUT" ]]; then
        echo "  Clearing ${OUT}..."
        rm -rf "${OUT:?}/"*
    fi
    echo -n "  Submitting ${VARIANT}... "
    JOB_ID=$(sbatch --parsable -A sdp166 "$SUBMIT_SCRIPT" "$VARIANT" "$TARGET_NS")
    echo "job ${JOB_ID}"
done

echo ""
echo "Monitor: squeue -u \$USER"
