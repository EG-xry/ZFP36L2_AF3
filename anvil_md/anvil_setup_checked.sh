#!/bin/bash
# ============================================================================
# Idempotent Anvil (Purdue, A100) environment + project setup for ZFP36L2 MD.
#
# Safe to re-run: every step CHECKS whether it is already done and skips it.
# Run this on the Anvil login node after `ssh -l YOUR_USER anvil.rcac.purdue.edu`.
#
#   bash anvil_setup_checked.sh            # setup + checks only
#   bash anvil_setup_checked.sh --status   # status report only, change nothing
#
# What it ensures exists:
#   ~/zfp36l2_md/{scripts,inputs,outputs,cleaned_pdbs,logs}
#   conda env 'openmm_env' with a working CUDA OpenMM
#   scripts uploaded (run_md_anvil.py, clean_af3_structures.py, submit_*.sh)
#   input PDBs present (best_models/*.pdb)
# ============================================================================
set -uo pipefail   # NOTE: not -e, so a single failed check doesn't abort the report

PROJECT_DIR="$HOME/zfp36l2_md"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
INPUT_DIR="${PROJECT_DIR}/inputs"
OUTPUT_DIR="${PROJECT_DIR}/outputs"
CLEAN_DIR="${PROJECT_DIR}/cleaned_pdbs"
LOG_DIR="${PROJECT_DIR}/logs"
ENV_NAME="openmm_env"

STATUS_ONLY=0
[[ "${1:-}" == "--status" ]] && STATUS_ONLY=1

ok()   { printf "  \033[32m[OK]\033[0m   %s\n" "$1"; }
todo() { printf "  \033[33m[TODO]\033[0m %s\n" "$1"; }
info() { printf "  \033[36m[..]\033[0m  %s\n" "$1"; }

echo "============================================================"
echo "ZFP36L2 MD — Anvil setup / status check"
echo "User: $(whoami)   Host: $(hostname)   Date: $(date)"
echo "============================================================"

# ----------------------------------------------------------------------------
# 1. Directory tree
# ----------------------------------------------------------------------------
echo ""
echo "[1] Project directories"
for d in "$PROJECT_DIR" "$SCRIPTS_DIR" "$INPUT_DIR" "$OUTPUT_DIR" "$CLEAN_DIR" "$LOG_DIR"; do
    if [[ -d "$d" ]]; then
        ok "exists: $d"
    else
        todo "missing: $d"
        [[ $STATUS_ONLY -eq 0 ]] && mkdir -p "$d" && ok "created: $d"
    fi
done

# ----------------------------------------------------------------------------
# 2. Modules (Anvil GPU stack)
# ----------------------------------------------------------------------------
echo ""
echo "[2] Modules"
if command -v module >/dev/null 2>&1; then
    if [[ $STATUS_ONLY -eq 0 ]]; then
        module purge 2>/dev/null || true
        # modtree/gpu pulls a GPU-aware toolchain on Anvil; anaconda + cuda for OpenMM
        module load modtree/gpu 2>/dev/null || todo "module load modtree/gpu failed (check 'module avail')"
        module load anaconda    2>/dev/null || todo "module load anaconda failed"
        module load cuda        2>/dev/null || todo "module load cuda failed"
        ok "loaded modtree/gpu, anaconda, cuda (see 'module list')"
    else
        info "status-only: skipping module load"
    fi
else
    todo "'module' command not found — are you on an Anvil node?"
fi

# ----------------------------------------------------------------------------
# 3. Conda environment
# ----------------------------------------------------------------------------
echo ""
echo "[3] Conda env '${ENV_NAME}'"
ENV_EXISTS=0
if command -v conda >/dev/null 2>&1; then
    if conda env list 2>/dev/null | grep -qE "^\s*${ENV_NAME}\s|/${ENV_NAME}\$"; then
        ENV_EXISTS=1
        ok "conda env '${ENV_NAME}' already exists"
    else
        todo "conda env '${ENV_NAME}' not found"
        if [[ $STATUS_ONLY -eq 0 ]]; then
            info "accepting conda TOS (idempotent) ..."
            conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
            conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true
            info "creating '${ENV_NAME}' (python 3.11, conda-forge) ..."
            conda create -n "${ENV_NAME}" -y -c conda-forge --override-channels python=3.11 \
                && ENV_EXISTS=1 && ok "created env"
        fi
    fi
else
    todo "'conda' not found — did 'module load anaconda' succeed?"
fi

# ----------------------------------------------------------------------------
# 4. OpenMM + CUDA inside the env  (check before installing)
# ----------------------------------------------------------------------------
echo ""
echo "[4] OpenMM + CUDA in '${ENV_NAME}'"
OPENMM_OK=0
CUDA_OK=0
if [[ $ENV_EXISTS -eq 1 ]]; then
    # shellcheck disable=SC1091
    source activate "${ENV_NAME}" 2>/dev/null || conda activate "${ENV_NAME}" 2>/dev/null || true
    if python -c "import openmm" 2>/dev/null; then
        OPENMM_OK=1
        VER=$(python -c "import openmm; print(openmm.__version__)" 2>/dev/null)
        ok "OpenMM importable (version ${VER})"
        if python -c "from openmm import Platform; Platform.getPlatformByName('CUDA')" 2>/dev/null; then
            CUDA_OK=1
            ok "OpenMM CUDA platform available"
        else
            todo "OpenMM present but CUDA platform NOT available"
        fi
    else
        todo "OpenMM not installed in '${ENV_NAME}'"
    fi

    if [[ $STATUS_ONLY -eq 0 && ( $OPENMM_OK -eq 0 || $CUDA_OK -eq 0 ) ]]; then
        info "installing OpenMM + CUDA-enabled build from conda-forge ..."
        # IMPORTANT: install openmm together with an explicit CUDA toolkit so
        # conda picks the CUDA-enabled openmm build. Installing openmm alone (or
        # pinning cuda-version after the fact) can leave you with a CPU-only
        # build whose CUDA platform is missing (the exact failure we hit).
        # On Anvil A100 the drivers are CUDA 12.x.
        conda install -y -c conda-forge --override-channels \
            openmm "cuda-version>=12.0,<13" numpy mdtraj pdbfixer 2>/dev/null \
        || conda install -y -c conda-forge --override-channels \
            openmm cudatoolkit numpy mdtraj pdbfixer 2>/dev/null \
        || true
        # Re-check
        if python -c "from openmm import Platform; Platform.getPlatformByName('CUDA')" 2>/dev/null; then
            CUDA_OK=1; OPENMM_OK=1; ok "OpenMM CUDA now available"
        else
            todo "CUDA platform still not visible on the LOGIN node."
            info "This is often EXPECTED: login nodes have no GPU, so the CUDA"
            info "plugin can't initialize here. It usually works on a GPU node."
            info "Verify on a GPU node instead (does not change anything):"
            info "  srun -A bio260224-gpu -p gpu-debug --gpus-per-node=1 -t 00:05:00 \\"
            info "       --pty bash -c 'module load modtree/gpu anaconda cuda; \\"
            info "       conda activate openmm_env; python -m openmm.testInstallation'"
        fi
    fi
fi

# ----------------------------------------------------------------------------
# 5. Scripts present
# ----------------------------------------------------------------------------
echo ""
echo "[5] Scripts in ${SCRIPTS_DIR}"
NEED_SCRIPTS=(run_md_anvil.py clean_af3_structures.py submit_variant.sh submit_all_variants.sh check_status.sh)
for s in "${NEED_SCRIPTS[@]}"; do
    if [[ -f "${SCRIPTS_DIR}/${s}" ]]; then
        ok "present: ${s}"
    else
        todo "missing: ${s}  (upload from your Mac — see UPLOAD section below)"
    fi
done

# ----------------------------------------------------------------------------
# 6. Input PDBs present
# ----------------------------------------------------------------------------
echo ""
echo "[6] Input PDBs in ${INPUT_DIR}"
N_PDB=$(find "${INPUT_DIR}" -maxdepth 1 -name "*.pdb" 2>/dev/null | wc -l | tr -d ' ')
if [[ "${N_PDB}" -ge 1 ]]; then
    ok "${N_PDB} PDB file(s) present"
    find "${INPUT_DIR}" -maxdepth 1 -name "*_best_model.pdb" 2>/dev/null \
        | sed 's#.*/#      #' || true
else
    todo "no PDBs found — upload best_models/*.pdb (see UPLOAD section below)"
fi

# ----------------------------------------------------------------------------
# 7. Per-variant run status (what has already completed?)
# ----------------------------------------------------------------------------
echo ""
echo "[7] Existing run progress in ${OUTPUT_DIR}"
VARIANTS=(WT Y154H R160W R184C G204D C206F S185N T187S D210E C212A C208A)
ANY=0
printf "      %-8s %-13s %10s %10s\n" VARIANT STATUS ns target
printf "      %s\n" "------------------------------------------------"
for V in "${VARIANTS[@]}"; do
    P="${OUTPUT_DIR}/${V}/progress.txt"
    if [[ -f "$P" ]]; then
        ANY=1
        st=$(grep '^status=' "$P" 2>/dev/null | cut -d= -f2)
        ns=$(grep '^ns_completed=' "$P" 2>/dev/null | cut -d= -f2)
        tg=$(grep '^target_ns=' "$P" 2>/dev/null | cut -d= -f2)
        printf "      %-8s %-13s %10s %10s\n" "$V" "${st:-?}" "${ns:-?}" "${tg:-?}"
    fi
done
[[ $ANY -eq 0 ]] && info "no runs started yet on Anvil"

# ----------------------------------------------------------------------------
# 8. Allocation sanity
# ----------------------------------------------------------------------------
echo ""
echo "[8] Allocation / quota"
command -v myquota >/dev/null 2>&1 && { myquota 2>/dev/null | head -5 || true; } || info "myquota not available"
command -v showpartitions >/dev/null 2>&1 && info "run 'showpartitions' to see gpu / gpu-debug limits"

# ----------------------------------------------------------------------------
# Summary + next steps
# ----------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "SUMMARY"
echo "============================================================"
[[ $ENV_EXISTS -eq 1 ]] && ok "conda env ready" || todo "create conda env (re-run without --status)"
[[ $OPENMM_OK  -eq 1 ]] && ok "OpenMM installed" || todo "install OpenMM (re-run without --status)"
[[ $CUDA_OK    -eq 1 ]] && ok "CUDA acceleration ready" || todo "fix CUDA platform"
[[ "${N_PDB:-0}" -ge 1 ]] && ok "inputs present" || todo "upload input PDBs"

cat <<'NEXT'

------------------------------------------------------------
UPLOAD (run these on your MAC if [5]/[6] showed TODO):
------------------------------------------------------------
  HOST=YOUR_USER@anvil.rcac.purdue.edu
  scp "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/run_md_anvil.py" \
      "$HOST:~/zfp36l2_md/scripts/"
  scp "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/clean_af3_structures.py" \
      "$HOST:~/zfp36l2_md/scripts/"
  scp "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/submit_variant.sh" \
      "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/submit_all_variants.sh" \
      "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/check_status.sh" \
      "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/anvil_submit_replicas.sh" \
      "$HOST:~/zfp36l2_md/scripts/"
  scp "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/best_models/"*.pdb \
      "$HOST:~/zfp36l2_md/inputs/"

------------------------------------------------------------
RUN (on Anvil, after uploads):
------------------------------------------------------------
  cd ~/zfp36l2_md && chmod +x scripts/*.sh
  # smoke-test one variant on the gpu-debug partition first:
  sbatch -p gpu-debug scripts/submit_variant.sh WT 1
  squeue -u $USER
  # full replica campaign (3 replicas x 200 ns, all variants):
  bash scripts/anvil_submit_replicas.sh 200 3
NEXT
echo ""
echo "Re-run this script any time:  bash anvil_setup_checked.sh --status"
