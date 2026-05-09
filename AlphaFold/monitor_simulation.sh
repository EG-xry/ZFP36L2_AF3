#!/bin/bash
# Script to monitor AF-MI simulation progress and save results

SIM_DIR="/home/hemiola/Desktop/AlphaFold/AlphaFold-IDP/prep_run"
RESULTS_DIR="/home/hemiola/Desktop/AlphaFold/simulation_results"

mkdir -p "$RESULTS_DIR"

echo "=== AlphaFold Metainference Simulation Monitor ==="
echo "Timestamp: $(date)"
echo ""

# Check if simulation is running
if pgrep -f "simulate.py" > /dev/null; then
    echo "✓ Simulation is RUNNING"
    echo "Processes:"
    ps aux | grep "simulate.py" | grep -v grep
else
    echo "✗ Simulation is NOT running"
fi

echo ""
echo "=== Simulation Outputs ==="

# Check for output files
if [ -f "$SIM_DIR/COLVAR.0" ]; then
    echo "✓ COLVAR files found:"
    ls -lh "$SIM_DIR"/COLVAR.* 2>/dev/null | head -5
    echo ""
    echo "Latest COLVAR entries:"
    tail -5 "$SIM_DIR/COLVAR.0" 2>/dev/null
else
    echo "✗ No COLVAR files yet"
fi

echo ""
if [ -f "$SIM_DIR/output_0.dcd" ]; then
    echo "✓ Trajectory files found:"
    ls -lh "$SIM_DIR"/output_*.dcd 2>/dev/null
else
    echo "✗ No trajectory files yet"
fi

echo ""
if [ -f "$SIM_DIR/stats_0.csv" ]; then
    echo "✓ Statistics files found:"
    ls -lh "$SIM_DIR"/stats_*.csv 2>/dev/null
else
    echo "✗ No statistics files yet"
fi

echo ""
echo "=== Log Files ==="
if [ -f "$SIM_DIR/af_mi_simulation.log" ]; then
    echo "Latest log entries:"
    tail -10 "$SIM_DIR/af_mi_simulation.log"
else
    echo "No simulation log found"
fi

echo ""
echo "=== Disk Usage ==="
df -h / | tail -1

echo ""
echo "=== Next Steps ==="
echo "1. Monitor progress: ./monitor_simulation.sh"
echo "2. Check logs: tail -f $SIM_DIR/af_mi_simulation.log"
echo "3. When complete, results will be in: $SIM_DIR"
echo "4. Analysis scripts available in: $SIM_DIR/../scripts_prep"

