#!/bin/bash
# Complete PLUMED build script with output logging

set -e

CONDA_ROOT="/Users/eric/miniconda3"
CONDA_ENV="$CONDA_ROOT/envs/afmi-tutorial"
BUILD_DIR="/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build"
SOURCE_DIR="$BUILD_DIR/plumed2-2.8.2"
LOG_FILE="$BUILD_DIR/build.log"

mkdir -p "$BUILD_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== PLUMED 2.8.2 Build with MPI ==="
echo "Date: $(date)"
echo "Conda env: $CONDA_ENV"
echo "Build dir: $BUILD_DIR"

# Install MPICH if not present
echo "--- Checking/Installing MPICH ---"
if [ ! -f "$CONDA_ENV/bin/mpicc" ]; then
    echo "MPICH not found, installing..."
    # Note: mpich-devel not available on macOS ARM64, just install mpich
    "$CONDA_ROOT/bin/conda" install -y -c conda-forge mpich
else
    echo "MPICH already installed"
fi

# Set environment
export PATH="$CONDA_ENV/bin:$CONDA_ROOT/bin:$PATH"
export CC="mpicc"
export CXX="mpicxx"

echo "CC: $CC"
echo "CXX: $CXX"
which mpicc
which mpicxx
mpicc --version | head -2
mpicxx --version | head -2

# Extract if needed
cd "$BUILD_DIR"
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Extracting source..."
    tar -xzf plumed-2.8.2.tar.gz
fi

# Configure
cd "$SOURCE_DIR"
echo "Configuring..."
./configure --prefix="$CONDA_ENV" --enable-mpi

# Build
echo "Building..."
make -j4

# Install
echo "Installing..."
make install

# Verify
echo "Verifying..."
if [ -f "$CONDA_ENV/bin/plumed" ]; then
    $CONDA_ENV/bin/plumed --has-mpi && echo "SUCCESS: MPI enabled" || echo "WARNING: MPI not detected"
else
    echo "ERROR: plumed binary not found"
fi

echo "=== Build complete ==="

