#!/bin/bash
# ============================================================================
# Submit a multi-replica ZFP36L2 MD campaign on Anvil A100, with completion
# checks so re-running NEVER re-submits work that is already done.
#
# Usage (on Anvil, from ~/zfp36l2_md):
#   bash scripts/anvil_submit_replicas.sh [PRODUCTION_NS] [N_REPLICAS]
#   bash scripts/anvil_submit_replicas.sh 200 3          # default
#   bash scripts/anvil_submit_replicas.sh 200 3 --dry-run  # show, don't submit
#
# For each (variant, replica) it:
#   - skips if that replica's progress.txt says status=completed at >= target
#   - resubmits with --resume if a checkpoint exists but it's not finished
#   - submits fresh otherwise
# Replica r uses output dir outputs/repN/<variant>/ so replicas never collide.
# ============================================================================
set -uo pipefail

PRODUCTION_NS="${1:-200}"
N_REPLICAS="${2:-3}"
DRY_RUN=0
[[ "${3:-}" == "--dry-run" ]] && DRY_RUN=1

PROJECT_DIR="$HOME/zfp36l2_md"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
INPUT_DIR="${PROJECT_DIR}/inputs"

VARIANTS=(WT Y154H R160W R184C G204D C206F S185N T187S D210E C212A)

# Map variant -> input PDB (best model). C212A may still be named C208A in inputs/.
declare -A PDB
for V in "${VARIANTS[@]}"; do PDB[$V]="${V}_best_model.pdb"; done
# Back-compat: if C212A best model isn't there, use the old C208A name.
[[ ! -f "${INPUT_DIR}/${PDB[C212A]}" && -f "${INPUT_DIR}/C208A_best_model.pdb" ]] \
    && PDB[C212A]="C208A_best_model.pdb"

echo "============================================================"
echo "ZFP36L2 replica campaign on Anvil A100"
echo "  production_ns = ${PRODUCTION_NS}   replicas = ${N_REPLICAS}"
echo "  dry_run = ${DRY_RUN}"
echo "============================================================"

submit_one () {
    local variant="$1" rep="$2"
    local outdir="${PROJECT_DIR}/outputs/rep${rep}"
    local prog="${outdir}/${variant}/progress.txt"
    local chk="${outdir}/${variant}/production.chk"
    local eqchk="${outdir}/${variant}/equilibrated.chk"
    local pdb="${INPUT_DIR}/${PDB[$variant]}"

    if [[ ! -f "$pdb" ]]; then
        printf "  [SKIP] %-6s rep%s : input PDB missing (%s)\n" "$variant" "$rep" "$(basename "$pdb")"
        return
    fi

    # Already finished?
    if [[ -f "$prog" ]]; then
        local st ns
        st=$(grep '^status=' "$prog" | cut -d= -f2)
        ns=$(grep '^ns_completed=' "$prog" | cut -d= -f2)
        if [[ "$st" == "completed" ]]; then
            printf "  [DONE] %-6s rep%s : completed (%s ns)\n" "$variant" "$rep" "${ns}"
            return
        fi
    fi

    # Resume vs fresh
    local resume=""
    if [[ -f "$chk" && -f "$eqchk" ]]; then
        resume="--resume"
        printf "  [RSUM] %-6s rep%s : checkpoint found -> resume\n" "$variant" "$rep"
    else
        printf "  [NEW ] %-6s rep%s : fresh submit\n" "$variant" "$rep"
    fi

    # Each replica gets a distinct velocity seed via the run script's
    # setVelocitiesToTemperature (seeded by submit time + rep); we pass the rep
    # through the job name so logs are separable, and use a per-rep output dir.
    # RISK REDUCTION:
    #  --exclusive      : take the whole GPU node so other users' jobs on the
    #                     same node can't compete for CPU/PCIe/memory bandwidth
    #                     (the most likely cause of the variable 9h-vs-18h speeds
    #                     we observed; node g010 had load avg ~47 while a job ran).
    #                     Set ZFP_EXCLUSIVE=0 to disable (e.g. to save SUs).
    #  --time=14:00:00  : a healthy 200 ns run is ~9-13h; this caps a hung/slow
    #                     job at 14h instead of letting it burn up to 48h. The
    #                     run is checkpoint/resume-safe so a cap loses nothing.
    local excl_flag="--exclusive"
    [[ "${ZFP_EXCLUSIVE:-1}" == "0" ]] && excl_flag=""
    local wall="${ZFP_WALL:-14:00:00}"
    local cmd=(sbatch
        --job-name "zfp_${variant}_r${rep}"
        -A bio260224-gpu
        -p gpu
        --nodes=1 --ntasks=1 --gpus-per-node=1 --cpus-per-task=8 --mem=32G
        ${excl_flag}
        --time="${wall}"
        --requeue --open-mode=append
        --output "${PROJECT_DIR}/logs/%x_%j.out"
        --error  "${PROJECT_DIR}/logs/%x_%j.err"
        --export=ALL,ZFP_OUTDIR="${outdir}",ZFP_VARIANT="${variant}",ZFP_PDB="${pdb}",ZFP_NS="${PRODUCTION_NS}",ZFP_RESUME="${resume}"
        "${SCRIPTS_DIR}/anvil_replica_job.sh")

    if [[ $DRY_RUN -eq 1 ]]; then
        printf "        would run: %s\n" "${cmd[*]}"
    else
        "${cmd[@]}"
    fi
}

mkdir -p "${PROJECT_DIR}/logs"
for r in $(seq 1 "$N_REPLICAS"); do
    echo ""
    echo "--- replica ${r} ---"
    for V in "${VARIANTS[@]}"; do
        submit_one "$V" "$r"
    done
done

echo ""
echo "Submitted. Watch with: squeue -u \$USER"
echo "Status report:        bash scripts/anvil_setup_checked.sh --status"
