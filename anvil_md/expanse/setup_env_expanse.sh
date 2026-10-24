#!/bin/bash
# ============================================================================
# Setup conda environment on SDSC Expanse for OpenMM GPU MD
# Run this ONCE after first login to Expanse
# ============================================================================

set -euo pipefail

echo "============================================================"
echo "Setting up OpenMM environment on SDSC Expanse"
echo "============================================================"

# Install Miniconda to a DIFFERENT location to avoid NFS lock issues
# The old ~/miniconda3 has NFS locks from a failed install - just leave it
CONDA_DIR="$HOME/conda_openmm"

if [[ -f "$CONDA_DIR/bin/conda" ]]; then
    echo "Conda already installed at $CONDA_DIR"
else
    echo "Installing Miniconda3 to $CONDA_DIR ..."
    # Download miniconda installer
    curl -sL https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda_install.sh
    # Install to new directory (non-interactive)
    bash /tmp/miniconda_install.sh -b -p "$CONDA_DIR"
    rm -f /tmp/miniconda_install.sh
    echo "Miniconda installed at $CONDA_DIR"
fi

# Initialize conda for this shell session
eval "$($CONDA_DIR/bin/conda shell.bash hook)"
echo "Conda version: $(conda --version)"

# Accept conda TOS if needed
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

# Remove old env if it exists but is broken
conda env remove -n openmm_env 2>/dev/null || true

# Create conda environment
echo "Creating conda environment 'openmm_env'..."
conda create -n openmm_env -y -c conda-forge --override-channels python=3.11

# Activate
conda activate openmm_env

# Install OpenMM + numpy from conda-forge
echo "Installing OpenMM with CUDA support..."
conda install -y -c conda-forge --override-channels openmm numpy cudatoolkit

# Add conda init to bashrc so it auto-loads in job scripts
$CONDA_DIR/bin/conda init bash 2>/dev/null || true

# Verify
echo ""
echo "Verifying OpenMM installation..."
python -c "
import openmm
print(f'OpenMM version: {openmm.__version__}')
from openmm import Platform
for i in range(Platform.getNumPlatforms()):
    p = Platform.getPlatform(i)
    print(f'  Platform: {p.getName()}')
try:
    cuda = Platform.getPlatformByName('CUDA')
    print('CUDA platform found! (may only work on GPU nodes)')
except:
    print('NOTE: CUDA not detected on login node - this is normal.')
    print('      CUDA will be available on GPU compute nodes.')
"

echo ""
echo "============================================================"
echo "Environment setup complete!"
echo "Activate with: conda activate openmm_env"
echo "============================================================"
