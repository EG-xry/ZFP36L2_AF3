#!/bin/bash
# Quick script to check PLUMED build progress

CONDA_ENV="/Users/eric/miniconda3/envs/afmi-tutorial"
BUILD_DIR="/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build"
SOURCE_DIR="$BUILD_DIR/plumed2-2.8.2"

echo "=== PLUMED Build Status Check ==="
echo ""

# Check if build script is running
echo "1. Build process status:"
if ps aux | grep -E "build_plumed|make.*plumed" | grep -v grep > /dev/null; then
    echo "   ✓ Build process is RUNNING"
    ps aux | grep -E "build_plumed|make.*plumed" | grep -v grep | head -2
else
    echo "   ✗ Build process NOT running"
fi
echo ""

# Check configuration
echo "2. Configuration status:"
if [ -f "$SOURCE_DIR/Makefile" ]; then
    echo "   ✓ Configuration SUCCESSFUL"
    grep -E "^prefix|has.*mpi" "$SOURCE_DIR/Makefile" | head -3
else
    echo "   ✗ Configuration NOT complete"
fi
echo ""

# Check build progress
echo "3. Build progress:"
if [ -f "$SOURCE_DIR/build_final.log" ]; then
    echo "   Build log exists"
    tail -5 "$SOURCE_DIR/build_final.log" | grep -E "compiling|linking|built|Error" | tail -3
    if [ -f "$SOURCE_DIR/src/lib/plumed" ] || [ -f "$SOURCE_DIR/src/lib/.libs/plumed" ]; then
        echo "   ✓ PLUMED binary built!"
    else
        echo "   ⏳ Build in progress..."
    fi
else
    echo "   Build not started yet"
fi
echo ""

# Check installation
echo "4. Installation status:"
if [ -f "$CONDA_ENV/bin/plumed" ]; then
    echo "   ✓ PLUMED installed at: $CONDA_ENV/bin/plumed"
    if "$CONDA_ENV/bin/plumed" --has-mpi 2>&1 | grep -q "mpi\|MPI"; then
        echo "   ✓✓✓ MPI support ENABLED!"
    else
        echo "   ⚠ MPI support status unknown"
    fi
else
    echo "   ✗ PLUMED not installed yet"
fi
echo ""

# Show recent log entries
echo "5. Recent build log entries:"
if [ -f "$BUILD_DIR/../plumed_build_monitor.log" ]; then
    tail -20 "$BUILD_DIR/../plumed_build_monitor.log" | grep -E "Step|===|SUCCESS|FAILED|✓|✗" | tail -10
else
    echo "   Monitor log not found"
fi

echo ""
echo "=== End Status Check ==="

