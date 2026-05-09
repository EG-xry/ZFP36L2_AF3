#!/bin/bash
# Script to free disk space and install build tools
# Run with: bash install_build_tools.sh

set -e

echo "=========================================="
echo "AlphaFold-IDP: Free Space & Install Build Tools"
echo "=========================================="
echo ""

# Step 1: Clean conda cache
echo "Step 1: Cleaning conda package cache..."
export PATH="$HOME/miniconda3/bin:$PATH"
source $HOME/miniconda3/etc/profile.d/conda.sh 2>/dev/null || true
conda clean --packages -y || echo "Note: Some packages may be in use"

echo ""
echo "Step 2: System cleanup (requires sudo password)..."
echo "You will be prompted for your sudo password."

# Step 2: System cleanup
sudo apt clean
sudo apt autoclean
sudo apt autoremove -y
sudo journalctl --vacuum-time=7d 2>/dev/null || echo "Note: journalctl cleanup skipped"

echo ""
echo "Step 3: Checking disk space..."
df -h / | tail -1

echo ""
echo "Step 4: Installing build tools (requires sudo password)..."
sudo apt update
sudo apt install -y build-essential g++ make

echo ""
echo "Step 5: Verifying installation..."
g++ --version
make --version
cmake --version

echo ""
echo "=========================================="
echo "✓ Build tools installed successfully!"
echo "=========================================="
echo ""
echo "You can now proceed with building OpenMM-Plumed-MPI"
echo "Check FREE_SPACE_AND_BUILD.md for next steps"

