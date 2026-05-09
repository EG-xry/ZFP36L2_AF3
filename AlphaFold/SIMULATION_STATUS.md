# AlphaFold-IDP Simulation Status

## Current Issue: OpenMM Topology Matching Error

### Problem
The energy minimization script (`simulate_em.py`) is failing with the error:
```
ValueError: No template found for residue 1 (MET). The set of atoms matches MET, but the bonds are different. Perhaps the chain is missing a terminal group?
```

### Root Cause
OpenMM's ForceField template matching is checking the topology structure before the script completes its modifications. The coarse-grained forcefield expects:
- Each residue to have one atom (the residue name)
- Two ExternalBond entries (bonds to adjacent residues)
- Terminal residues named X and Z (not MET, etc.)

The script adds bonds between consecutive residues, but OpenMM checks the topology before those bonds are fully recognized in the template matching process.

### Current Setup Status
✅ All dependencies installed (OpenMM, PLUMED, MPI, etc.)
✅ OpenMM-Plumed-MPI built successfully
✅ Coarse-grained PDB file created (`input_af.pdb` with 414 CA atoms)
✅ Forcefield file present (`forcefield.xml`)
✅ PLUMED files generated (`plumed.dat`, `plumed_analysis.dat`)
✅ Sequence and domain files ready

### Next Steps to Resolve

1. **Option 1: Modify the script to handle topology creation differently**
   - Create topology with bonds from the start
   - Ensure terminal residues are properly named before OpenMM checks

2. **Option 2: Check if there's a working example**
   - Look for example input files in the repository
   - Verify the exact PDB format expected

3. **Option 3: Use a different approach**
   - Create the system manually without relying on template matching
   - Or modify the forcefield to be more lenient

### Files Ready for Simulation
- `input_af.pdb` - Coarse-grained structure (414 residues)
- `forcefield.xml` - CALVADOS forcefield
- `plumed.dat` - PLUMED input for metainference
- `plumed_analysis.dat` - PLUMED analysis file
- `sequence.dat` - Protein sequence
- `r1_excl.pkl` - Exclusion list for forces
- `residues.csv` - Residue parameters

### Commands to Run (once issue is resolved)
```bash
cd /home/hemiola/Desktop/AlphaFold/AlphaFold-IDP/prep_run
conda activate alphafold-idp
python simulate_em.py 7.4 298  # Energy minimization
mpirun -np 2 python simulate.py 7.4 298  # AF-MI simulation
```

## AF3 Protein/RNA Prep (2025-12-03)

- Added `scripts_prep/clean_af3_structures.py` to automate chain renaming (protein→A, RNA→B), drop waters/ions, and (re)add Zn²⁺ sites by clustering Cys/His ligands so the Amber14 forcefields load without template errors.
- Generated cleaned inputs for the four GM-CSF 17-nt complexes in `prep_run/cleaned_pdbs/` with `_cleaned.pdb` suffix; each file now contains exactly one protein chain, one RNA chain, and no extraneous heterogens.
- WT, D210E, and GC-mutant RNAs regained two Zn sites (Cys9/18/24 + His28, Cys47/56/62 + His66) while C206F only supports the N-terminal Zn finger because the Cys47/56/62 triad is disrupted by the mutation.
- Zn residues are appended as `HETATM ... ZN A <71/72>` directly after the protein atoms, so they will be parameterized via Amber14 (`amber14-all.xml`, `amber14-nucleic-acid.OL3.xml`, `tip3p.xml`) without extra topology edits.

## Minimal MD Driver (2025-12-03)

- Added `prep_run/run_minimal_md.py`, which solvates each cleaned complex, minimizes, performs the staged NVT/NPT equilibration, and runs a 10 ns production (2 fs timestep, 2 ps DCD stride) with optional Zn–Cys restraints.
- Default targets are the four GM-CSF 17-nt systems; outputs land in `prep_run/minimal_md/<pdb_stem>/` as `md.log`, `production.dcd`, and `final_frame.pdb`.
- The script exposes knobs for platform selection, solvation padding, stride lengths, restraint strengths, and duration overrides so “Step 3-lite” can run with a single CLI call.

