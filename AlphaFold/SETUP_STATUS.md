# AlphaFold-IDP Setup Status Report

## ✅ Completed Steps

### 1. Environment Setup
- ✅ Miniconda installed and configured
- ✅ Conda environment `alphafold-idp` created (Python 3.11.6)
- ✅ All required packages installed:
  - OpenMM 8.0.0
  - PLUMED 2.8.2 + py-plumed
  - MPICH 4.3.2 + MPI4py 4.1.1
  - GROMACS 2025.3
  - CMake 4.2.0, SWIG 4.4.0
  - pandas, numpy, matplotlib, mdtraj, biopython
  - ipykernel

### 2. Repository Setup
- ✅ OpenMM-Plumed-MPI cloned
- ✅ AlphaFold-IDP cloned
- ✅ All required zip files present

### 3. Protein System Setup
- ✅ Changed to `prep_run` directory
- ✅ Created `sequence.dat` with TDP-43 sequence
- ✅ Created `ordered_domains.csv` and `disordered_domains.csv`
- ✅ Copied AF data files:
  - `pdb_af.pdb` (AlphaFold PDB structure)
  - `pae.json` (Predicted Aligned Error)
  - `mean_af.csv` (Mean distance matrix)
  - `prob.npy` (Probability distributions)
- ✅ Copied script files:
  - `gen_xml_and_constraints.py`
  - `residues.csv`
- ✅ Generated forcefield XML and constraints files
- ✅ Generated `plumed_analysis.dat`
- ✅ Copied TDP-43 specific PLUMED files:
  - `plumed_TDP-43.dat` → `plumed.dat`
  - `plumed_analysis_TDP-43.dat` → `plumed_analysis.dat`

## ⚠️ Current Blockers

### 1. OpenMM-Plumed-MPI Build Required
**Status**: Not built yet  
**Reason**: Requires C++ compiler and more disk space

**What's needed**:
- C++ compiler (g++ or clang++)
- Build tools (make, cmake - already have cmake)
- ~500MB-1GB free disk space for build

**Current disk status**: 344M free (98% used)

### 2. Missing Module: `openmmplumed`
The simulation scripts require `from openmmplumed import PlumedForce`, which is provided by the OpenMM-Plumed-MPI build.

## 📋 Next Steps Required

### Option 1: Build OpenMM-Plumed-MPI (Recommended)

1. **Free up disk space** (need ~1GB):
   ```bash
   # Run these commands manually (some require sudo):
   sudo apt clean
   sudo apt autoclean
   sudo apt autoremove -y
   sudo journalctl --vacuum-time=7d
   ```

2. **Install build tools** (requires sudo):
   ```bash
   sudo apt update
   sudo apt install -y build-essential g++ make
   ```

3. **Build OpenMM-Plumed-MPI**:
   ```bash
   cd /home/hemiola/Desktop/AlphaFold/OpenMM-Plumed-MPI
   cd build
   
   export PATH="$HOME/miniconda3/bin:$PATH"
   source $HOME/miniconda3/etc/profile.d/conda.sh
   conda activate alphafold-idp
   
   cmake .. \
     -DCMAKE_INSTALL_PREFIX=$(pwd)/../install \
     -DOPENMM_DIR=$(pwd)/../openmm \
     -DPLUMED_INCLUDE_DIR=$(pwd)/../plumed/include \
     -DPLUMED_LIBRARY_DIR=$(pwd)/../plumed/lib \
     -DMPI4PY_DIR=$HOME/miniconda3/envs/alphafold-idp/lib/python3.11/site-packages/mpi4py/include
   
   make -j$(nproc)
   make install
   make PythonInstall
   
   # Copy libraries
   cd ../install/lib
   cp -r * $HOME/miniconda3/envs/alphafold-idp/lib/
   ```

### Option 2: Use Pre-built OpenMM-Plumed (if available)
Check if there's a pre-built version available for ARM64 architecture.

## 📁 Current File Structure

```
/home/hemiola/Desktop/AlphaFold/
├── AlphaFold-IDP/
│   ├── prep_run/
│   │   ├── sequence.dat ✅
│   │   ├── ordered_domains.csv ✅
│   │   ├── disordered_domains.csv ✅
│   │   ├── pdb_af.pdb ✅
│   │   ├── pae.json ✅
│   │   ├── mean_af.csv ✅
│   │   ├── prob.npy ✅
│   │   ├── plumed.dat ✅
│   │   ├── plumed_analysis.dat ✅
│   │   ├── forcefield.xml ✅
│   │   ├── residues.csv ✅
│   │   └── simulate_em.py ✅
│   └── AF_DATA/ ✅
└── OpenMM-Plumed-MPI/
    ├── build/ (created, ready for cmake)
    ├── install/ (created)
    ├── openmm/ (unzipped) ✅
    └── plumed/ (unzipped) ✅
```

## 🔍 Verification Commands

After building OpenMM-Plumed-MPI, verify installation:

```bash
cd /home/hemiola/Desktop/AlphaFold/AlphaFold-IDP/prep_run
export PATH="$HOME/miniconda3/bin:$PATH"
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate alphafold-idp

# Test import
python -c "from openmmplumed import PlumedForce; print('✓ openmmplumed OK')"

# Run energy minimization
python simulate_em.py 7.4 298
```

## 📝 Notes

- All setup files are in place and ready
- The only missing piece is the OpenMM-Plumed-MPI build
- Once built, you can proceed with:
  1. Energy minimization (`simulate_em.py`)
  2. AF-MI simulation (`simulate.py` with mpirun)
  3. Analysis steps

## 🎯 Summary

**Progress**: ~90% complete  
**Remaining**: Build OpenMM-Plumed-MPI (requires disk space and build tools)  
**Ready to run**: Once OpenMM-Plumed-MPI is built, all simulation scripts are ready

