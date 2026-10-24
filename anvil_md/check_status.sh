#!/bin/bash
# ============================================================================
# Check status of all ZFP36L2 MD simulation jobs on Anvil
# ============================================================================

PROJECT_DIR="$HOME/zfp36l2_md"
OUTPUT_DIR="${PROJECT_DIR}/outputs"

echo "============================================================"
echo "ZFP36L2 MD Simulation Status Report"
echo "Date: $(date)"
echo "============================================================"
echo ""

# --- SLURM Queue Status ---
echo "=== SLURM Queue Status ==="
squeue -u $USER -o "%.10i %.12j %.10P %.8T %.10M %.10l %.6D %R" 2>/dev/null || echo "Not on Anvil or squeue unavailable"
echo ""

# --- Per-Variant Progress ---
echo "=== Per-Variant Progress ==="
printf "%-10s %-12s %10s %10s %10s %s\n" "Variant" "Status" "ns_done" "target_ns" "%" "Last Update"
echo "----------------------------------------------------------------------"

VARIANTS=("WT" "Y154H" "R160W" "R184C" "G204D" "C206F" "S185N" "T187S" "D210E" "C208A")

TOTAL_NS=0
TOTAL_TARGET=0
COMPLETED=0
RUNNING=0
PENDING=0

for V in "${VARIANTS[@]}"; do
    PROG_FILE="${OUTPUT_DIR}/${V}/progress.txt"
    if [[ -f "$PROG_FILE" ]]; then
        STATUS=$(grep "^status=" "$PROG_FILE" 2>/dev/null | cut -d= -f2 || echo "unknown")
        NS_DONE=$(grep "^ns_completed=" "$PROG_FILE" 2>/dev/null | cut -d= -f2 || echo "0")
        TARGET=$(grep "^target_ns=" "$PROG_FILE" 2>/dev/null | cut -d= -f2 || echo "0")
        TIMESTAMP=$(grep "^timestamp=" "$PROG_FILE" 2>/dev/null | cut -d= -f2 | cut -dT -f1 || echo "?")
        
        if [[ "$TARGET" != "0" ]]; then
            PCT=$(echo "scale=1; $NS_DONE / $TARGET * 100" | bc 2>/dev/null || echo "?")
        else
            PCT="?"
        fi
        
        printf "%-10s %-12s %10s %10s %9s%% %s\n" "$V" "$STATUS" "$NS_DONE" "$TARGET" "$PCT" "$TIMESTAMP"
        
        TOTAL_NS=$(echo "$TOTAL_NS + $NS_DONE" | bc 2>/dev/null || echo "$TOTAL_NS")
        TOTAL_TARGET=$(echo "$TOTAL_TARGET + $TARGET" | bc 2>/dev/null || echo "$TOTAL_TARGET")
        
        if [[ "$STATUS" == "completed" ]]; then
            COMPLETED=$((COMPLETED + 1))
        elif [[ "$STATUS" == "checkpointed" ]]; then
            RUNNING=$((RUNNING + 1))
        fi
    else
        # Check if job is in queue
        JOB_STATE=$(squeue -u $USER -n "zfp36l2_MD" -o "%T" --noheader 2>/dev/null | head -1 || echo "")
        if [[ -n "$JOB_STATE" ]]; then
            printf "%-10s %-12s %10s %10s %10s %s\n" "$V" "queued" "-" "-" "-" "-"
            PENDING=$((PENDING + 1))
        else
            printf "%-10s %-12s %10s %10s %10s %s\n" "$V" "not_started" "-" "-" "-" "-"
        fi
    fi
done

echo "----------------------------------------------------------------------"
echo ""
echo "=== Summary ==="
echo "  Completed:    ${COMPLETED}/10"
echo "  Checkpointed: ${RUNNING}"
echo "  Pending:      ${PENDING}"
echo "  Total ns:     ${TOTAL_NS} / ${TOTAL_TARGET}"
echo ""

# --- GPU Hours Used ---
echo "=== GPU Hours Used (approximate) ==="
sacct -u $USER --format=JobID,JobName,Elapsed,State,AllocGRES --starttime=$(date -d '3 days ago' +%Y-%m-%d 2>/dev/null || date -v-3d +%Y-%m-%d 2>/dev/null || echo "2026-05-22") 2>/dev/null | grep -i "zfp36l2\|gpu" || echo "sacct not available"
echo ""

# --- Disk Usage ---
echo "=== Disk Usage ==="
if [[ -d "$OUTPUT_DIR" ]]; then
    du -sh "${OUTPUT_DIR}"/* 2>/dev/null || echo "No output directories yet"
    echo ""
    du -sh "${OUTPUT_DIR}" 2>/dev/null
fi
echo ""

# --- Check for errors ---
echo "=== Recent Errors (last 20 lines of .err files) ==="
for errfile in ${PROJECT_DIR}/logs/*.err; do
    if [[ -f "$errfile" && -s "$errfile" ]]; then
        echo "--- $(basename $errfile) ---"
        tail -5 "$errfile"
        echo ""
    fi
done 2>/dev/null || echo "No error files found"
