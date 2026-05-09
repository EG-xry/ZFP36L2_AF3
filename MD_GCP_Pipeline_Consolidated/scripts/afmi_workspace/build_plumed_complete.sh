#!/bin/bash
# Complete PLUMED build script with monitoring

set -e

CONDA_ROOT="/Users/eric/miniconda3"
CONDA_ENV="$CONDA_ROOT/envs/afmi-tutorial"
BUILD_DIR="/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build"
SOURCE_DIR="$BUILD_DIR/plumed2-2.8.2"
LOG_FILE="$BUILD_DIR/build_complete.log"

mkdir -p "$BUILD_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo "PLUMED 2.8.2 Build with MPI - Complete"
echo "Date: $(date)"
echo "=========================================="

# Activate conda environment
export PATH="$CONDA_ENV/bin:$CONDA_ROOT/bin:$PATH"

# Step 1: Install compilers
echo ""
echo "=== Step 1: Installing C++ compilers ==="
"$CONDA_ROOT/bin/conda" install -y -c conda-forge compilers cxx-compiler 2>&1 | tail -20
echo "Compilers installed"

# Step 2: Install MPICH
echo ""
echo "=== Step 2: Installing MPICH ==="
if [ ! -f "$CONDA_ENV/bin/mpicc" ]; then
    "$CONDA_ROOT/bin/conda" install -y -c conda-forge mpich 2>&1 | tail -20
    echo "MPICH installed"
else
    echo "MPICH already installed"
fi

# Step 3: Verify MPI compilers
echo ""
echo "=== Step 3: Verifying MPI compilers ==="
export CC=mpicc
export CXX=mpicxx
which mpicc mpicxx
mpicc --version | head -2
mpicxx --version | head -2

# Test compiler
echo ""
echo "=== Step 4: Testing compiler ==="
echo 'int main(){return 0;}' > /tmp/test_mpicxx.cpp
TEST_OUTPUT=$(mpicxx /tmp/test_mpicxx.cpp -o /tmp/test_mpicxx 2>&1)
TEST_EXIT=$?
if [ $TEST_EXIT -eq 0 ]; then
    if /tmp/test_mpicxx 2>&1; then
        echo "✓ Compiler test PASSED"
        rm -f /tmp/test_mpicxx.cpp /tmp/test_mpicxx
    else
        echo "✗ Compiler test FAILED - executable doesn't run"
        echo "Output: $TEST_OUTPUT"
        exit 1
    fi
else
    echo "✗ Compiler test FAILED - compilation error"
    echo "Error output: $TEST_OUTPUT"
    echo "Attempting to continue anyway (may work for PLUMED)..."
    rm -f /tmp/test_mpicxx.cpp /tmp/test_mpicxx 2>/dev/null || true
fi

# Step 5: Configure PLUMED
echo ""
echo "=== Step 5: Configuring PLUMED ==="
cd "$SOURCE_DIR"
make distclean 2>/dev/null || true
echo "Running configure..."
./configure --prefix="$CONDA_ENV" --enable-mpi 2>&1 | tee configure_final.log

if [ ! -f Makefile ]; then
    echo "✗ Configuration FAILED"
    echo "Last 50 lines of config.log:"
    tail -50 config.log
    exit 1
fi

echo "✓ Configuration SUCCESSFUL"
grep -E "^prefix|has.*mpi|MPI.*yes" Makefile | head -5

# Step 6: Build PLUMED
echo ""
echo "=== Step 6: Building PLUMED (this takes 10-30 minutes) ==="
echo "Start time: $(date)"
make -j4 2>&1 | tee build_final.log

if [ $? -ne 0 ]; then
    echo "✗ Build FAILED"
    echo "Last 50 lines of build log:"
    tail -50 build_final.log
    exit 1
fi

echo "✓ Build SUCCESSFUL"
echo "End time: $(date)"

# Step 7: Install PLUMED
echo ""
echo "=== Step 7: Installing PLUMED ==="
make install 2>&1 | tee install_final.log
echo "✓ Installation complete"

# Step 8: Verify
echo ""
echo "=== Step 8: Verifying installation ==="
if [ -f "$CONDA_ENV/bin/plumed" ]; then
    echo "PLUMED binary found: $CONDA_ENV/bin/plumed"
    if "$CONDA_ENV/bin/plumed" --has-mpi 2>&1; then
        echo "✓✓✓ SUCCESS: PLUMED has MPI support! ✓✓✓"
    else
        echo "⚠ WARNING: MPI support not detected"
    fi
else
    echo "✗ ERROR: plumed binary not found"
    exit 1
fi

echo ""
echo "=========================================="
echo "PLUMED build complete!"
echo "=========================================="

