# OpenMM-Plumed-MPI Compatibility Issue - Solution Attempts

## Problem
`TypeError: in method 'System_addForce', argument 2 of type 'OpenMM::Force *'`

## Root Cause Analysis

The SWIG Python bindings for `PlumedForce` were generated against:
- Extracted `openmm.zip` headers (from the repository)
- `simtk.openmm` module structure (old API)

But we're using:
- Conda-installed OpenMM 8.0 (new API: `openmm`)
- Different header/library versions

This creates a type mismatch in the SWIG wrapper.

## Attempted Solutions

1. ✅ Set LD_LIBRARY_PATH - Libraries found but issue persists
2. ✅ Copy libraries to conda lib directory - No change
3. ✅ Uninstall PyPI openmmplumed - No change  
4. ✅ Rebuild Python bindings - Failed (MPI not found in new build)
5. ❌ Rebuild against conda OpenMM - Failed (MPI configuration issue)

## Key Findings

From the [vendruscolo-lab repository README](https://github.com/vendruscolo-lab/OpenMM-Plumed-MPI):
- Uses extracted `openmm.zip` for headers during build
- Requires specific MPI setup
- Installation instructions: `cp -r * ~/miniconda3/envs/pl/lib`

The SWIG interface (`python/plumedplugin.i`) imports:
```swig
%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
```

This suggests it expects the old `simtk.openmm` module structure.

## Possible Solutions

### Option 1: Use simtk.openmm Compatibility Layer
OpenMM 8.0 still provides `simtk.openmm` for compatibility. Try modifying the script to use it:

```python
import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as unit
```

### Option 2: Check OpenMM Version Match
The extracted `openmm.zip` might be a different version than conda OpenMM 8.0. Check:
```bash
unzip -l OpenMM-Plumed-MPI/openmm.zip | grep version
```

### Option 3: Rebuild with Correct MPI
Fix MPI detection and rebuild:
```bash
export MPI_C_COMPILER=mpicc
export MPI_CXX_COMPILER=mpicxx
# Then rebuild
```

### Option 4: Use Pre-built Binary
Check if there's a pre-built version for ARM64 or if we need to use x86_64 emulation.

## Next Steps

1. Try modifying `simulate.py` to use `simtk.openmm` imports
2. Check the version of OpenMM in `openmm.zip`
3. Fix MPI configuration and rebuild
4. Consider if this is an ARM64-specific issue

