#!/bin/bash
# ============================================================================
# Per-job SLURM body for one (variant, replica) ZFP36L2 MD run on Anvil A100.
# Invoked by anvil_submit_replicas.sh via sbatch with these env vars exported:
#   ZFP_OUTDIR   per-replica output dir   (e.g. ~/zfp36l2_md/outputs/rep2)
#   ZFP_VARIANT  variant name             (e.g. R160W)
#   ZFP_PDB      absolute path to input PDB
#   ZFP_NS       production ns            (e.g. 200)
#   ZFP_RESUME   "--resume" or ""
#
# It cleans the PDB if needed, loads modules, activates the env, and runs
# run_md_anvil.py. Auto-resume is honored if a checkpoint exists.
# ============================================================================
set -euo pipefail

PROJECT_DIR="$HOME/zfp36l2_md"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
CLEAN_DIR="${PROJECT_DIR}/cleaned_pdbs"
mkdir -p "${CLEAN_DIR}" "${ZFP_OUTDIR}"

VARIANT="${ZFP_VARIANT:?}"
INPUT_PDB="${ZFP_PDB:?}"
PRODUCTION_NS="${ZFP_NS:?}"
RESUME_FLAG="${ZFP_RESUME:-}"

# --- Modules + env (MUST come before any python call, e.g. PDB cleaning) ---
module purge
module load modtree/gpu
module load anaconda
module load cuda
# Source the conda profile so `conda activate` works in a non-interactive batch
# shell (plain `source activate` can leave python off PATH -> "python: command
# not found"). Then hard-verify python+openmm before doing any work.
CONDA_BASE="$(conda info --base 2>/dev/null)"
if [[ -n "${CONDA_BASE}" && -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
fi
conda activate openmm_env
echo "Using python: $(command -v python || echo none)"
if ! python -c "import openmm" 2>/dev/null; then
    echo "ERROR: python/openmm not on PATH after activating openmm_env (base=${CONDA_BASE})"
    exit 3
fi

# --- Clean the AF3 PDB once per variant (shared across replicas) ---
# NOTE: this runs AFTER conda activate so `python` is guaranteed on PATH.
CLEANED_PDB="${CLEAN_DIR}/${VARIANT}_cleaned.pdb"
if [[ ! -f "$CLEANED_PDB" ]]; then
    echo "[clean] cleaning ${INPUT_PDB} -> ${CLEANED_PDB}"
    python "${SCRIPTS_DIR}/clean_af3_structures.py" --output-dir "${CLEAN_DIR}" "${INPUT_PDB}" || true
    if [[ ! -f "$CLEANED_PDB" ]]; then
        # clean script may name file differently; grab the newest cleaned pdb
        CLEANED_PDB=$(find "${CLEAN_DIR}" -maxdepth 1 -name "*cleaned.pdb" -newer "${INPUT_PDB}" 2>/dev/null | head -1)
    fi
    [[ -f "$CLEANED_PDB" ]] || { echo "ERROR: cleaning failed for ${VARIANT}"; exit 2; }
else
    echo "[clean] reusing existing cleaned PDB: ${CLEANED_PDB}"
fi

echo "============================================================"
echo "ZFP36L2 MD (Anvil A100)"
echo "  variant=${VARIANT}  replica_outdir=${ZFP_OUTDIR}"
echo "  production_ns=${PRODUCTION_NS}  resume=${RESUME_FLAG:-no}"
echo "  node=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
echo "  job=${SLURM_JOB_ID:-local}  date=$(date)"
echo "============================================================"

# --- Run ---
# run_md_anvil.py writes outputs to <output-dir>/<variant>/; we point it at the
# per-replica ZFP_OUTDIR so replicas never overwrite each other.
python "${SCRIPTS_DIR}/run_md_anvil.py" \
    --pdb "${CLEANED_PDB}" \
    --variant "${VARIANT}" \
    --output-dir "${ZFP_OUTDIR}" \
    --production-ns "${PRODUCTION_NS}" \
    --padding 1.2 \
    --platform CUDA \
    --precision single \
    --wall-hours 47.5 \
    ${RESUME_FLAG}

echo "Done: ${VARIANT} (replica dir ${ZFP_OUTDIR}) at $(date)"
