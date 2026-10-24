#!/bin/bash
# ============================================================================
# Re-submit variants that did not reach status=completed.
# Uses --skip-equil automatically when equilibrated.pdb already exists.
# ============================================================================

set -euo pipefail

PROJECT_DIR="$HOME/zfp36l2_md"
PRODUCTION_NS="${1:-50}"
SUBMIT_SCRIPT="${PROJECT_DIR}/submit_variant_expanse.sh"

if [[ ! -f "$SUBMIT_SCRIPT" ]]; then
    SUBMIT_SCRIPT="${PROJECT_DIR}/expanse/submit_variant_expanse.sh"
fi

VARIANTS=(WT Y154H R160W R184C G204D C206F S185N T187S D210E C208A)

echo "Checking for incomplete variants..."
TO_RERUN=()
for VARIANT in "${VARIANTS[@]}"; do
    PROGRESS="${PROJECT_DIR}/outputs/${VARIANT}/progress.txt"
    if [[ ! -f "$PROGRESS" ]] || ! grep -q "^status=completed" "$PROGRESS" 2>/dev/null; then
        TO_RERUN+=("$VARIANT")
        echo "  - ${VARIANT} (incomplete or missing)"
    else
        echo "  - ${VARIANT} (completed, skipping)"
    fi
done

if [[ ${#TO_RERUN[@]} -eq 0 ]]; then
    echo "All variants completed. Nothing to rerun."
    exit 0
fi

echo ""
echo "Submitting ${#TO_RERUN[@]} variant(s) at ${PRODUCTION_NS} ns..."
for VARIANT in "${TO_RERUN[@]}"; do
    echo -n "  ${VARIANT}... "
    JOB_ID=$(sbatch --parsable -A sdp166 "$SUBMIT_SCRIPT" "$VARIANT" "$PRODUCTION_NS")
    echo "job ${JOB_ID}"
done

echo ""
echo "Monitor with: squeue -u \$USER"
