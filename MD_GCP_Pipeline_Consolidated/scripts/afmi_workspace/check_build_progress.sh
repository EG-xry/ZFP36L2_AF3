#!/bin/bash
# Comprehensive build progress checker

BUILD_DIR="/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build/plumed2-2.8.2"

echo "=========================================="
echo "PLUMED Build Progress Check"
echo "=========================================="
echo ""

cd "$BUILD_DIR" || exit 1

# Check if build is running
echo "1. Build Process Status:"
if ps aux | grep -E "make.*plumed|make -j4" | grep -v grep > /dev/null; then
    echo "   ✓ Build is RUNNING"
    ps aux | grep -E "make.*plumed|make -j4" | grep -v grep | head -2 | awk '{print "   PID:", $2, "CPU:", $3"%", "Time:", $10}'
else
    echo "   ✗ Build is NOT running"
fi
echo ""

# Check build logs
echo "2. Build Logs:"
for log in build_direct.log build_monitored.log build_final.log build.log; do
    if [ -f "$log" ]; then
        lines=$(wc -l < "$log" 2>/dev/null || echo 0)
        echo "   ✓ $log exists ($lines lines)"
    fi
done
echo ""

# Check build progress from log
if [ -f build_monitored.log ]; then
    echo "3. Build Progress (from build_monitored.log):"
    total_lines=$(wc -l < build_monitored.log)
    compiling=$(grep -c "compiling\|Compiling" build_monitored.log 2>/dev/null || echo 0)
    linking=$(grep -c "linking\|Linking" build_monitored.log 2>/dev/null || echo 0)
    errors=$(grep -c "Error\|error\|FAILED\|failed" build_monitored.log 2>/dev/null || echo 0)
    
    echo "   Total log lines: $total_lines"
    echo "   Compiling operations: $compiling"
    echo "   Linking operations: $linking"
    echo "   Errors found: $errors"
    echo ""
    
    echo "4. Recent Activity (last 20 lines):"
    tail -20 build_monitored.log | head -20
    echo ""
    
    echo "5. Recent Errors (if any):"
    tail -100 build_monitored.log | grep -E "Error|error|FAILED|failed" | tail -5 || echo "   No recent errors"
    echo ""
fi

# Check if PLUMED binary exists
echo "6. PLUMED Binary Status:"
if [ -f src/lib/plumed ] || [ -f src/lib/.libs/plumed ]; then
    echo "   ✓ PLUMED binary FOUND!"
    find src -name "plumed" -type f 2>&1 | head -3 | sed 's/^/   /'
else
    echo "   ✗ PLUMED binary not found yet"
fi
echo ""

# Check object files
echo "7. Build Artifacts:"
obj_files=$(find . -name "*.o" -type f 2>/dev/null | wc -l | tr -d ' ')
echo "   Object files (.o): $obj_files"
lib_files=$(find . -name "*.a" -o -name "*.so" -o -name "*.dylib" 2>/dev/null | wc -l | tr -d ' ')
echo "   Library files: $lib_files"
echo ""

# Check installation
echo "8. Installation Status:"
if [ -f /Users/eric/miniconda3/envs/afmi-tutorial/bin/plumed ]; then
    echo "   ✓ PLUMED installed at: /Users/eric/miniconda3/envs/afmi-tutorial/bin/plumed"
    if /Users/eric/miniconda3/envs/afmi-tutorial/bin/plumed --has-mpi 2>&1 | grep -q "mpi\|MPI"; then
        echo "   ✓✓✓ MPI support ENABLED!"
    else
        echo "   ⚠ MPI status unknown"
    fi
else
    echo "   ✗ PLUMED not installed yet"
fi

echo ""
echo "=========================================="

