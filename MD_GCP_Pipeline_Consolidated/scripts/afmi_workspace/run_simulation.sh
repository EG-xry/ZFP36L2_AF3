#!/bin/bash
# Wrapper script to run AF-MI simulation with correct library paths

export DYLD_LIBRARY_PATH=/Users/eric/micromamba/envs/omm/lib:$DYLD_LIBRARY_PATH

cd /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/AlphaFold-IDP/prep_run

# Activate conda environment and run simulation
source /Users/eric/miniconda3/etc/profile.d/conda.sh
conda activate afmi-tutorial

# Run simulation with MPI
mpirun -np 2 python simulate.py 7.4 298

