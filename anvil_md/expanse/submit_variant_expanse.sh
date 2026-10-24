#!/bin/bash
#SBATCH --job-name=zfp36l2_MD
#SBATCH -A bio260224
#SBATCH -p gpu-shared

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --no-requeue

# ============================================================================
# ZFP36L2 MD Simulation - SDSC Expanse V100 GPU
# Same architecture as Anvil version, adapted for Expanse SLURM syntax
# ============================================================================

set -euo pipefail

# --- Configuration ---
VARIANT="${1:?Usage: sbatch [--time=HH:MM:SS] submit_variant_expanse.sh <VARIANT> [PRODUCTION_NS] [WALL_HOURS]}"
PRODUCTION_NS="${2:-50}"
WALL_HOURS="${3:-23.5}"
RESUME_FLAG=""

# Project paths on Expanse
PROJECT_DIR="$HOME/zfp36l2_md"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
INPUT_DIR="${PROJECT_DIR}/inputs"
# Use ZFP_OUTDIR if set (for replica runs), otherwise default to outputs/rep2
# (separate from the old 400ns runs in outputs/)
OUTPUT_DIR="${ZFP_OUTDIR:-${PROJECT_DIR}/outputs/rep2}"
CLEAN_DIR="${PROJECT_DIR}/cleaned_pdbs"

mkdir -p "${OUTPUT_DIR}" "${CLEAN_DIR}" logs

# --- Map variant to input PDB ---
declare -A VARIANT_PDBS_BEST
VARIANT_PDBS_BEST[WT]="WT_best_model.pdb"
VARIANT_PDBS_BEST[Y154H]="Y154H_best_model.pdb"
VARIANT_PDBS_BEST[R160W]="R160W_best_model.pdb"
VARIANT_PDBS_BEST[R184C]="R184C_best_model.pdb"
VARIANT_PDBS_BEST[G204D]="G204D_best_model.pdb"
VARIANT_PDBS_BEST[C206F]="C206F_best_model.pdb"
VARIANT_PDBS_BEST[S185N]="S185N_best_model.pdb"
VARIANT_PDBS_BEST[T187S]="T187S_best_model.pdb"
VARIANT_PDBS_BEST[D210E]="D210E_best_model.pdb"
VARIANT_PDBS_BEST[C208A]="C208A_best_model.pdb"

declare -A VARIANT_PDBS_FALLBACK
VARIANT_PDBS_FALLBACK[WT]="fold_zfp36l2_wt_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[Y154H]="fold_zfp36l2_y154h_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[R160W]="fold_zfp36l2_r160w_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[R184C]="fold_zfp36l2_r184c_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[G204D]="fold_zfp36l2_g204d_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[C206F]="fold_zfp36l2_c206f_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[S185N]="fold_zfp36l2_s185n_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[T187S]="fold_zfp36l2_t187s_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[D210E]="fold_zfp36l2_d210e_nonamer_wt_2zn_model_0_run1.pdb"
VARIANT_PDBS_FALLBACK[C208A]="fold_zfp36l2_deleterious_c208a_model_0_run1.pdb"

# Try best model first, fall back to run1
INPUT_PDB="${INPUT_DIR}/${VARIANT_PDBS_BEST[$VARIANT]}"
if [[ ! -f "$INPUT_PDB" ]]; then
    echo "Best model not found, falling back to run1_model0..."
    INPUT_PDB="${INPUT_DIR}/${VARIANT_PDBS_FALLBACK[$VARIANT]}"
fi

if [[ ! -f "$INPUT_PDB" ]]; then
    echo "ERROR: Input PDB not found: $INPUT_PDB"
    exit 1
fi

# --- Check for existing checkpoint (auto-resume) ---
CHECKPOINT="${OUTPUT_DIR}/${VARIANT}/production.chk"
EQUIL_CHK="${OUTPUT_DIR}/${VARIANT}/equilibrated.chk"
EQUIL_PDB="${OUTPUT_DIR}/${VARIANT}/equilibrated.pdb"
PROGRESS="${OUTPUT_DIR}/${VARIANT}/progress.txt"
TRAJ="${OUTPUT_DIR}/${VARIANT}/production.dcd"
RESUME_FLAG=""
SKIP_EQUIL_FLAG=""

if [[ -f "$CHECKPOINT" && -f "$EQUIL_CHK" ]]; then
    echo "Found production checkpoint for ${VARIANT}. Enabling resume."
    RESUME_FLAG="--resume"
elif [[ -f "$EQUIL_PDB" ]]; then
    if [[ ! -f "$PROGRESS" ]] || ! grep -q "^status=completed" "$PROGRESS" 2>/dev/null; then
        echo "Found equilibrated structure for ${VARIANT}. Skipping equilibration, restarting production."
        SKIP_EQUIL_FLAG="--skip-equil"
        rm -f "$CHECKPOINT"
        if [[ -f "$TRAJ" ]] && [[ $(wc -c < "$TRAJ") -lt 100 ]]; then
            rm -f "$TRAJ"
        fi
    fi
fi

# --- Clean the PDB (if not already cleaned) ---
# Use variant-specific output name to avoid cross-contamination
CLEANED_PDB="${CLEAN_DIR}/${VARIANT}_cleaned.pdb"
if [[ ! -f "$CLEANED_PDB" ]]; then
    echo "Cleaning AF3 structure for ${VARIANT}..."
    # Clean to a temp dir to avoid find picking up wrong files
    TEMP_CLEAN="/tmp/clean_${VARIANT}_$$"
    mkdir -p "$TEMP_CLEAN"
    python "${SCRIPTS_DIR}/clean_af3_structures.py" \
        --output-dir "$TEMP_CLEAN" "$INPUT_PDB"
    # Move the cleaned file with correct variant name
    FOUND_CLEAN=$(find "$TEMP_CLEAN" -maxdepth 1 -name "*cleaned.pdb" | head -1)
    if [[ -z "$FOUND_CLEAN" || ! -f "$FOUND_CLEAN" ]]; then
        echo "ERROR: Cleaning failed for ${VARIANT}"
        rm -rf "$TEMP_CLEAN"
        exit 2
    fi
    cp "$FOUND_CLEAN" "$CLEANED_PDB"
    rm -rf "$TEMP_CLEAN"
    echo "Cleaned PDB: $CLEANED_PDB"
fi

# --- Load conda (using Miniconda installed by setup_env_expanse.sh) ---
eval "$($HOME/conda_openmm/bin/conda shell.bash hook)"
conda activate openmm_env

# --- Print diagnostics ---
echo "============================================================"
echo "ZFP36L2 MD Simulation on SDSC Expanse V100"
echo "============================================================"
echo "Variant:        ${VARIANT}"
echo "Input PDB:      ${INPUT_PDB}"
echo "Cleaned PDB:    ${CLEANED_PDB}"
echo "Production:     ${PRODUCTION_NS} ns"
echo "Wall hours:     ${WALL_HOURS} h (python safety margin applied)"
echo "Output:         ${OUTPUT_DIR}/${VARIANT}"
echo "Resume:         ${RESUME_FLAG:-no}"
echo "Skip equil:     ${SKIP_EQUIL_FLAG:-no}"
echo "SLURM Job ID:   ${SLURM_JOB_ID:-local}"
echo "Node:           $(hostname)"
echo "GPU:            $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Date:           $(date)"
echo "============================================================"

# --- Run simulation ---
python "${SCRIPTS_DIR}/run_md_anvil.py" \
    --pdb "$CLEANED_PDB" \
    --variant "$VARIANT" \
    --output-dir "$OUTPUT_DIR" \
    --production-ns "$PRODUCTION_NS" \
    --padding 1.2 \
    --platform CUDA \
    --precision single \
    --wall-hours "$WALL_HOURS" \
    $SKIP_EQUIL_FLAG \
    $RESUME_FLAG

echo ""
echo "Job completed at $(date)"
echo "Output in: ${OUTPUT_DIR}/${VARIANT}"
