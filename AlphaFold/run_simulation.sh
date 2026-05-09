#!/bin/bash
# Script to run AF-MI simulation with proper environment setup

cd /home/hemiola/Desktop/AlphaFold/AlphaFold-IDP/prep_run

# Use the dedicated alphafold-idp environment
ENV_PREFIX="$HOME/miniconda3/envs/alphafold-idp"
export PATH="$ENV_PREFIX/bin:$PATH"
export OPENMM_PLUGIN_DIR="$ENV_PREFIX/lib/plugins"
export LD_LIBRARY_PATH="$ENV_PREFIX/lib:$LD_LIBRARY_PATH"

# Set number of replicas (NR from notebook)
NR=2
pH=7.4
temp=298

echo "=== Starting AlphaFold Metainference Simulation ==="
echo "Replicas: $NR"
echo "pH: $pH"
echo "Temperature: $temp K"
echo "Timestamp: $(date)"
echo ""

# Run simulation
mpirun -np $NR "$ENV_PREFIX/bin/python" simulate.py $pH $temp 2>&1 | tee af_mi_simulation.log

echo ""
echo "=== Simulation Complete ==="
echo "Check outputs: COLVAR.*, output_*.dcd, stats_*.csv"

