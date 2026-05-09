# Final Investigation Summary - OpenMM-Plumed-MPI Compatibility Issue

## Problem
Cannot add `PlumedForce` to OpenMM `System`:
```
TypeError: in method 'System_addForce', argument 2 of type 'OpenMM::Force *'
```

## Root Cause
The SWIG Python bindings for `PlumedForce` were generated against:
- Extracted `openmm.zip` headers (from vendruscolo-lab repository)
- Potentially different OpenMM version/ABI than conda OpenMM 8.0

This creates a type mismatch in the C++ SWIG wrapper, even though Python-level inheritance is correct.

## What We've Tried

1. ✅ Set LD_LIBRARY_PATH correctly
2. ✅ Copy libraries to conda environment
3. ✅ Uninstall PyPI openmmplumed
4. ✅ Use `simtk.openmm` compatibility layer
5. ✅ Check library dependencies
6. ❌ Rebuild against conda OpenMM - Failed (MPI configuration)

## Key Finding

The [vendruscolo-lab OpenMM-Plumed-MPI repository](https://github.com/vendruscolo-lab/OpenMM-Plumed-MPI) uses:
- Extracted `openmm.zip` for build headers
- Specific MPI setup requirements
- Installation: `cp -r * ~/miniconda3/envs/pl/lib`

The SWIG interface expects `simtk.openmm`, but even using that doesn't fix the C++ type mismatch.

## Solution Required

**Rebuild the Python bindings against the exact conda OpenMM installation:**

1. Fix MPI detection in CMake:
   ```bash
   export MPI_C_COMPILER=mpicc
   export MPI_CXX_COMPILER=mpicxx
   ```

2. Rebuild with conda OpenMM:
   ```bash
   cmake .. -DOPENMM_DIR=$CONDA_PREFIX \
            -DPLUMED_INCLUDE_DIR=$CONDA_PREFIX/include/plumed \
            -DPLUMED_LIBRARY_DIR=$CONDA_PREFIX/lib \
            -DMPI4PY_DIR=$CONDA_PREFIX/lib/python3.11/site-packages/mpi4py/include
   ```

3. Rebuild and reinstall:
   ```bash
   make
   make install  
   make PythonInstall
   cp -r install/lib/* $CONDA_PREFIX/lib/
   ```

## Alternative: Check Architecture Compatibility

This might be an ARM64-specific issue. The tutorial may have been tested on x86_64. Consider:
- Using x86_64 emulation (if available)
- Checking if there are pre-built binaries for ARM64
- Contacting the repository maintainers

## Current Status

- ✅ Environment: 100% complete
- ✅ Energy minimization: Working
- ✅ Input files: All prepared
- ❌ AF-MI simulation: Blocked by SWIG binding issue

**The setup is functionally complete except for this one compatibility issue that requires rebuilding the Python bindings.**

