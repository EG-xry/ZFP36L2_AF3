#!/bin/bash
# Step-by-step commands to restart PLUMED build
# Run these commands one by one or execute this script

cd /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build/plumed2-2.8.2

# Step 1: Activate conda and find MPI compilers
source /Users/eric/miniconda3/etc/profile.d/conda.sh
conda activate afmi-tutorial

# Step 2: Check MPI compilers
which mpicc mpicxx
mpicc --version | head -2

# Step 3: Set compiler environment
export CC=$(which mpicc)
export CXX=$(which mpicxx)
export PATH=/Users/eric/miniconda3/envs/afmi-tutorial/bin:$PATH

# Step 4: Clean previous build (if exists)
if [ -f Makefile ]; then
    make distclean > /dev/null 2>&1 || true
fi

# Step 5: Configure PLUMED
./configure --prefix=/Users/eric/miniconda3/envs/afmi-tutorial --enable-mpi 2>&1 | tee configure_now.log

# Step 6: Verify configuration succeeded
if [ ! -f config.status ]; then
    echo "ERROR: Configuration failed!"
    tail -30 configure_now.log
    exit 1
fi

# Step 7: Start build in background
nohup make -j4 > build_running.log 2>&1 &
BUILD_PID=$!

echo "Build started with PID: $BUILD_PID"
echo "Monitor with: tail -f build_running.log"
echo "Check status: bash /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/MONITOR_BUILD.sh"

# Step 8: Wait a moment and check it's running
sleep 5
if ps -p $BUILD_PID > /dev/null 2>&1; then
    echo "✓ Build is running"
    tail -20 build_running.log
else
    echo "✗ Build failed to start"
    tail -30 build_running.log
fi

