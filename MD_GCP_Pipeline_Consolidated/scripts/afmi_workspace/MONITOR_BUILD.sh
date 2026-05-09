#!/bin/bash
# Simple build monitor script

BUILD_DIR="/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build/plumed2-2.8.2"
LOG_FILE="$BUILD_DIR/build_running.log"

cd "$BUILD_DIR" || exit 1

echo "=== PLUMED Build Monitor ==="
echo ""

# Check if build is running
if ps aux | grep -E "[m]ake.*plumed|[m]ake -j4" | grep -v grep > /dev/null; then
    echo "✓ Build is RUNNING"
    ps aux | grep -E "[m]ake.*plumed|[m]ake -j4" | grep -v grep | head -1 | awk '{print "  PID:", $2, "CPU:", $3"%"}'
else
    echo "✗ Build is NOT running"
fi
echo ""

# Check log file
if [ -f "$LOG_FILE" ]; then
    lines=$(wc -l < "$LOG_FILE")
    size=$(du -h "$LOG_FILE" | cut -f1)
    echo "Build log: $LOG_FILE"
    echo "  Lines: $lines"
    echo "  Size: $size"
    echo ""
    
    # Count operations
    compiling=$(grep -c "compiling\|Compiling" "$LOG_FILE" 2>/dev/null || echo 0)
    linking=$(grep -c "linking\|Linking" "$LOG_FILE" 2>/dev/null || echo 0)
    errors=$(grep -c "Error\|error\|FAILED" "$LOG_FILE" 2>/dev/null || echo 0)
    
    echo "Progress:"
    echo "  Compiling: $compiling operations"
    echo "  Linking: $linking operations"
    echo "  Errors: $errors"
    echo ""
    
    echo "Last 10 lines:"
    tail -10 "$LOG_FILE" | sed 's/^/  /'
else
    echo "Build log not found: $LOG_FILE"
fi
echo ""

# Check for PLUMED binary
if [ -f src/lib/plumed ] || [ -f src/lib/.libs/plumed ]; then
    echo "✓✓✓ PLUMED BUILD COMPLETE! ✓✓✓"
    ls -lh src/lib/plumed* src/lib/.libs/plumed* 2>&1 | head -3
else
    echo "PLUMED binary not found yet"
fi

