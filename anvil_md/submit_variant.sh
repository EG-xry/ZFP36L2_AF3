#!/bin/bash
#SBATCH --job-name=zfp36l2_MD
#SBATCH -A bio260224-gpu
#SBATCH -p gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
# NOTE: default 4 h fits gpu-debug (which has a short QOS wall limit, ~30min–4h).
# For the real long runs, submit to -p gpu and override with: --time=48:00:00
# e.g.  sbatch -p gpu --time=48:00:00 scripts/submit_variant.sh WT 200

#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=YOUR_EMAIL

# ============================================================================
# ZFP36L2 MD Simulation - Anvil A100 GPU
# Single-variant SLURM job (use --array for multiple variants)
# ============================================================================

set -euo pipefail

# --- Configuration ---
VARIANT="${1:?Usage: sbatch submit_variant.sh <VARIANT> [PRODUCTION_NS]}"
PRODUCTION_NS="${2:-100}"
RESUME_FLAG=""

# Project paths on Anvil
PROJECT_DIR="$HOME/zfp36l2_md"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
INPUT_DIR="${PROJECT_DIR}/inputs"
OUTPUT_DIR="${PROJECT_DIR}/outputs"
CLEAN_DIR="${PROJECT_DIR}/cleaned_pdbs"

mkdir -p "${OUTPUT_DIR}" "${CLEAN_DIR}" logs

# --- Map variant to input PDB ---
# Priority: use best-ranked model (from select_best_models.py), fall back to run1_model0
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
    echo "Available variants: ${!VARIANT_PDBS[@]}"
    exit 1
fi

# --- Load modules + activate conda FIRST (cleaning needs python too) ---
module purge
module load modtree/gpu
module load anaconda
module load cuda
CONDA_BASE="$(conda info --base 2>/dev/null)"
if [[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
fi
conda activate openmm_env
PYBIN="$(command -v python || true)"
echo "Using python: ${PYBIN:-<none>}"
if [[ -z "${PYBIN}" ]] || ! python -c "import openmm" 2>/dev/null; then
    echo "ERROR: python/openmm not on PATH after activating openmm_env."
    echo "       conda base = ${CONDA_BASE}"
    exit 3
fi

# --- Check for existing checkpoint (auto-resume) ---
CHECKPOINT="${OUTPUT_DIR}/${VARIANT}/production.chk"
EQUIL_CHK="${OUTPUT_DIR}/${VARIANT}/equilibrated.chk"
if [[ -f "$CHECKPOINT" && -f "$EQUIL_CHK" ]]; then
    echo "Found existing checkpoint for ${VARIANT}. Enabling resume."
    RESUME_FLAG="--resume"
fi

# --- Clean the PDB (if not already cleaned) ---
CLEANED_PDB="${CLEAN_DIR}/${VARIANT}_cleaned.pdb"
if [[ ! -f "$CLEANED_PDB" ]]; then
    echo "Cleaning AF3 structure for ${VARIANT}..."
    python "${SCRIPTS_DIR}/clean_af3_structures.py" \
        --output-dir "${CLEAN_DIR}" "$INPUT_PDB"
    # Find the cleaned file
    CLEANED_PDB=$(find "${CLEAN_DIR}" -maxdepth 1 -name "*cleaned.pdb" -newer "$INPUT_PDB" | head -1)
    if [[ -z "$CLEANED_PDB" || ! -f "$CLEANED_PDB" ]]; then
        echo "ERROR: Cleaning failed for ${VARIANT}"
        exit 2
    fi
fi

# --- Print diagnostics ---
echo "============================================================"
echo "ZFP36L2 MD Simulation on Anvil A100"
echo "============================================================"
echo "Variant:        ${VARIANT}"
echo "Input PDB:      ${INPUT_PDB}"
echo "Cleaned PDB:    ${CLEANED_PDB}"
echo "Production:     ${PRODUCTION_NS} ns"
echo "Output:         ${OUTPUT_DIR}/${VARIANT}"
echo "Resume:         ${RESUME_FLAG:-no}"
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
    --wall-hours 47.5 \
    $RESUME_FLAG

echo ""
echo "Job completed at $(date)"
echo "Output in: ${OUTPUT_DIR}/${VARIANT}"
