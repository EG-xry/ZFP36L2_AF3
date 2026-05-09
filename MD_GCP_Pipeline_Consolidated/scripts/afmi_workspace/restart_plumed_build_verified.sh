#!/bin/bash
# Verified PLUMED build restart script

set -e

BUILD_DIR="/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build/plumed2-2.8.2"
CONDA_ENV="/Users/eric/miniconda3/envs/afmi-tutorial"

cd "$BUILD_DIR" || { echo "Cannot cd to $BUILD_DIR"; exit 1; }

echo "=========================================="
echo "PLUMED Build Restart Script"
echo "=========================================="
echo "Build directory: $BUILD_DIR"
echo ""

# Activate conda environment properly
source /Users/eric/miniconda3/etc/profile.d/conda.sh
conda activate afmi-tutorial
export PATH="$CONDA_ENV/bin:$PATH"

# Ensure stub vim directory exists for optional docs target
if [ ! -d vim ]; then
    mkdir -p vim
    cat > vim/Makefile <<'EOF'
all:
	@true

install:
	@true

clean:
	@true
EOF
fi

# Find MPI compilers
MPICC=$(which mpicc 2>/dev/null)
MPICXX=$(which mpicxx 2>/dev/null)

# Verify MPI compilers
echo "=== Checking MPI compilers ==="
if [ -z "$MPICC" ] || [ -z "$MPICXX" ]; then
    echo "MPI compilers not found. Installing MPICH..."
    conda install -y -c conda-forge mpich 2>&1 | tail -10
    MPICC=$(which mpicc 2>/dev/null)
    MPICXX=$(which mpicxx 2>/dev/null)
fi

if [ -z "$MPICC" ] || [ -z "$MPICXX" ]; then
    echo "✗ ERROR: Cannot find MPI compilers"
    exit 1
fi

export CC="$MPICC"
export CXX="$MPICXX"
echo "✓ mpicc: $CC"
echo "✓ mpicxx: $CXX"
"$CC" --version 2>&1 | head -1
echo ""

# Step 1: Clean and reconfigure
echo "=== Step 1: Configuring PLUMED ==="
if [ -f Makefile ]; then
    echo "Cleaning previous build..."
    make distclean > /dev/null 2>&1 || true
fi

echo "Running configure..."
./configure --prefix="$CONDA_ENV" --enable-mpi 2>&1 | tee configure_restart.log | tail -30

if [ ! -f config.status ] || [ ! -f Makefile ]; then
    echo ""
    echo "✗ Configuration FAILED"
    echo "Last 30 lines:"
    tail -30 configure_restart.log
    exit 1
fi

echo ""
echo "✓ Configuration SUCCESSFUL"
grep -E "^prefix|has.*mpi|MPI.*yes" Makefile | head -3
echo ""

# Step 2: Start build
echo "=== Step 2: Starting build ==="
echo "Start time: $(date)"
echo ""

# Kill any existing make processes
pkill -f "make.*plumed" 2>/dev/null || true
sleep 1

# Start build in background
nohup make -j4 > build_running.log 2>&1 &
BUILD_PID=$!

echo "Build started with PID: $BUILD_PID"
echo "Log file: $BUILD_DIR/build_running.log"
echo ""

# Wait and verify
sleep 5
if ps -p $BUILD_PID > /dev/null 2>&1; then
    echo "✓ Build process is running"
    echo ""
    echo "To monitor: tail -f $BUILD_DIR/build_running.log"
    echo "To check status: bash /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/MONITOR_BUILD.sh"
    echo ""
    echo "Initial log output:"
    tail -15 build_running.log
else
    echo "✗ Build process failed to start"
    echo "Error log:"
    tail -30 build_running.log
    exit 1
fi

