#!/bin/bash
# Start PLUMED build with proper monitoring

set -e

BUILD_DIR="/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build/plumed2-2.8.2"
CONDA_ENV="/Users/eric/miniconda3/envs/afmi-tutorial"

cd "$BUILD_DIR" || exit 1

echo "=========================================="
echo "Starting PLUMED Build"
echo "=========================================="
echo ""

# Set environment
export PATH="$CONDA_ENV/bin:$PATH"
export CC=mpicc
export CXX=mpicxx

# Step 1: Reconfigure if needed
if [ ! -f config.status ] || [ ! -f Makefile ]; then
    echo "=== Step 1: Configuring PLUMED ==="
    ./configure --prefix="$CONDA_ENV" --enable-mpi 2>&1 | tee configure_final.log
    
    if [ ! -f config.status ]; then
        echo "✗ Configuration FAILED"
        tail -30 configure_final.log
        exit 1
    fi
    echo "✓ Configuration complete"
    echo ""
else
    echo "✓ Configuration already exists"
    echo ""
fi

# Step 2: Start build
echo "=== Step 2: Starting build ==="
echo "Start time: $(date)"
echo "Build will run in background"
echo "Monitor with: tail -f build_active.log"
echo ""

# Start build in background
nohup make -j4 > build_active.log 2>&1 &
BUILD_PID=$!

echo "Build started with PID: $BUILD_PID"
echo ""

# Wait a moment and check it started
sleep 3
if ps -p $BUILD_PID > /dev/null; then
    echo "✓ Build process is running (PID: $BUILD_PID)"
    echo ""
    echo "To monitor progress:"
    echo "  tail -f $BUILD_DIR/build_active.log"
    echo ""
    echo "To check status:"
    echo "  bash /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/MONITOR_BUILD.sh"
    echo ""
    echo "Expected build time: 10-30 minutes"
else
    echo "✗ Build process failed to start"
    echo "Check build_active.log for errors"
    exit 1
fi

