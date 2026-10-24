#!/bin/bash
# ============================================================================
# Setup conda environment on Anvil for OpenMM GPU MD
# Run this ONCE after first login to Anvil
# ============================================================================

set -euo pipefail

echo "============================================================"
echo "Setting up OpenMM environment on Anvil"
echo "============================================================"

# Load required modules
module purge
module load modtree/gpu
module load anaconda
module load cuda

# Accept conda Terms of Service (required on Anvil since mid-2025)
echo "Accepting conda Terms of Service..."
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main 2>/dev/null || true
conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r 2>/dev/null || true

# Create conda environment with OpenMM + CUDA support
# Using conda-forge ONLY to avoid Anaconda TOS channel issues
echo "Creating conda environment 'openmm_env'..."
conda create -n openmm_env -y -c conda-forge --override-channels python=3.11

# Activate
source activate openmm_env 2>/dev/null || conda activate openmm_env

# Install OpenMM with CUDA 12 support from conda-forge
# Anvil A100s use CUDA 12.x drivers - need matching openmm-cuda plugin
echo "Installing OpenMM with CUDA support..."
conda install -y -c conda-forge --override-channels openmm numpy

# Install the CUDA platform plugin for OpenMM
# The openmm package from conda-forge includes CUDA support when 
# cudatoolkit is available, but we need the right version matching Anvil's driver
echo "Installing OpenMM CUDA plugin..."
conda install -y -c conda-forge --override-channels cudatoolkit=11.8 2>/dev/null || \
conda install -y -c conda-forge --override-channels cudatoolkit 2>/dev/null || true

# Try installing openmm-cuda (separate package in newer builds)
pip install --no-deps openmm-cuda 2>/dev/null || true

# Verify installation
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
    print(f'CUDA platform found: {cuda.getName()}')
    print('GPU acceleration is available!')
except Exception as e:
    print(f'WARNING: CUDA not available: {e}')
    print('Jobs will fall back to CPU (much slower)')
"

echo ""
echo "============================================================"
echo "Environment setup complete!"
echo "Activate with: conda activate openmm_env"
echo "============================================================"
