#!/bin/bash
# Simple PLUMED build restart - use this one

cd /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build/plumed2-2.8.2

# Initialize conda and activate environment
source /Users/eric/miniconda3/etc/profile.d/conda.sh
conda activate afmi-tutorial

# Ensure stub vim directory exists to satisfy optional doc build
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

# Set compilers
export CC=$(which mpicc)
export CXX=$(which mpicxx)
export PATH=/Users/eric/miniconda3/envs/afmi-tutorial/bin:$PATH

echo "=== Checking MPI compilers ==="
which mpicc mpicxx
mpicc --version | head -2
echo ""

# Clean if needed
if [ -f Makefile ]; then
    echo "Cleaning previous build..."
    make distclean > /dev/null 2>&1 || true
fi

# Configure
echo "=== Configuring PLUMED ==="
./configure --prefix=/Users/eric/miniconda3/envs/afmi-tutorial --enable-mpi 2>&1 | tee configure_now.log

if [ ! -f config.status ]; then
    echo "✗ Configuration FAILED"
    tail -30 configure_now.log
    exit 1
fi

echo "✓ Configuration successful"
echo ""

# Start build
echo "=== Starting build ==="
nohup make -j4 > build_running.log 2>&1 &
BUILD_PID=$!

echo "Build started with PID: $BUILD_PID"
echo "Monitor with: tail -f build_running.log"
echo ""

sleep 5
if ps -p $BUILD_PID > /dev/null 2>&1; then
    echo "✓ Build is running"
    tail -15 build_running.log
else
    echo "✗ Build failed to start"
    tail -30 build_running.log
fi

