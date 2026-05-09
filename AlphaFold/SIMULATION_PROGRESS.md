# AlphaFold Metainference Simulation Progress

## Status: ✅ Energy Minimization Complete | 🔄 AF-MI Simulation Running

### Completed Steps

1. ✅ **Environment Setup**
   - All dependencies installed (OpenMM 8.0, PLUMED 2.8.2, MPI, etc.)
   - OpenMM-Plumed-MPI built and installed
   - Conda environment `alphafold-idp` configured

2. ✅ **Input Preparation**
   - Coarse-grained PDB created (`input_af.pdb` with 414 residues)
   - Terminal residues set to X and Z (required for OpenMM)
   - Forcefield file generated (`forcefield.xml`)
   - PLUMED files created (`plumed.dat`, `plumed_analysis.dat`)
   - Sequence and domain files ready

3. ✅ **Energy Minimization**
   - Successfully completed
   - Checkpoint file created: `checkpoint`
   - Log: `energy_minimization.log`

### Current Step: AF-MI Simulation

**Configuration:**
- Number of replicas: 2
- pH: 7.4
- Temperature: 298 K
- Simulation script: `simulate.py`

**Running Command:**
```bash
export LD_LIBRARY_PATH="/home/hemiola/Desktop/AlphaFold/OpenMM-Plumed-MPI/install/lib:$LD_LIBRARY_PATH"
mpirun -np 2 python simulate.py 7.4 298
```

**Expected Outputs:**
- `COLVAR.0`, `COLVAR.1` - Collective variable trajectories
- `output_0.dcd`, `output_1.dcd` - Trajectory files
- `stats_0.csv`, `stats_1.csv` - Statistics files
- `HILLS.*` - Metadynamics hills files
- `BAYES_*` - Metainference output files

### Monitoring

**Check simulation status:**
```bash
./monitor_simulation.sh
```

**View live log:**
```bash
tail -f /home/hemiola/Desktop/AlphaFold/AlphaFold-IDP/prep_run/af_mi_simulation.log
```

**Check if running:**
```bash
ps aux | grep simulate.py
```

### Next Steps After Simulation Completes

1. **Analysis** (from notebook cells):
   - Generate atomistic structural ensemble
   - Calculate free energy surfaces (FES)
   - Compute RMSF per residue
   - Create visualizations

2. **Analysis Scripts Available:**
   - `script.sh` - Generate structural ensemble
   - `fes2.py` - Calculate FES
   - `resample.py` - Resample by weights
   - `backmap.py` - Backmap to atomistic

### Files and Directories

- **Simulation directory:** `/home/hemiola/Desktop/AlphaFold/AlphaFold-IDP/prep_run`
- **Results:** Output files in simulation directory
- **Scripts:** `/home/hemiola/Desktop/AlphaFold/AlphaFold-IDP/scripts_prep`
- **Monitor script:** `/home/hemiola/Desktop/AlphaFold/monitor_simulation.sh`
- **Run script:** `/home/hemiola/Desktop/AlphaFold/run_simulation.sh`

### Notes

- The simulation may take a long time depending on system size and number of steps
- Monitor disk space: `df -h /`
- Results will be saved incrementally as the simulation progresses
- Use `Ctrl+C` to stop if needed (saves checkpoint)

