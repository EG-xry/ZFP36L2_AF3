# ZFP36L2 MD Simulations on Purdue Anvil A100

## Overview

This package runs molecular dynamics simulations of ZFP36L2 zinc-finger domain variants complexed with a 9-mer RNA and 2 Zn²⁺ ions using OpenMM on Anvil's NVIDIA A100 GPUs.

**Username:** `YOUR_USER`  
**Allocation:** 750 GPU-hours (ACCESS BIO240054)  
**System:** ZFP36L2 ZFD (70 residues) + RNA 9-mer (UUAUUUAUU) + 2 Zn²⁺ (~770 heavy atoms, ~15,000 solvated)

---

## 🔴 CRITICAL: Is 1000 ns Feasible in 3 Days?

### Short Answer: **Not all 10 variants at 1000 ns, but partial is possible.**

### Performance Estimates (A100 GPU, ~15k solvated atoms):

| Metric | Estimate |
|--------|----------|
| Speed on A100 (mixed precision) | ~200–400 ns/day |
| Best case per variant in 48h | ~800 ns |
| Conservative per variant in 48h | ~400 ns |
| 10 variants × 1000 ns | 10,000 GPU-hours needed |
| **Your budget** | **750 GPU-hours** |

### Recommendations:

| Strategy | Production (ns) | Total GPU-hrs | Feasibility |
|----------|----------------|---------------|-------------|
| **100 ns × 10 variants** | 100 | ~80-120 | ✅ Very safe, fits in budget |
| **200 ns × 10 variants** | 200 | ~160-240 | ✅ Safe |
| **500 ns × 10 variants** | 500 | ~400-600 | ⚠️ Tight but possible |
| 1000 ns × 10 variants | 1000 | ~800-1200 | ❌ Exceeds 750h budget |
| **1000 ns × 5 variants** | 1000 | ~400-600 | ⚠️ Possible if selective |

### My Recommendation: **100 ns per variant (default)**
- Well within your 750 GPU-hour budget (~120h total)
- Scientifically meaningful for comparing variant dynamics
- Each job completes in ~8-12 hours (well under 48h wall limit)
- Room for reruns or extensions if needed
- **Change to 200 ns** if you want more sampling (still safe)

### If you MUST have 1000 ns:
- Run only your most important 5 variants (WT + 2 deleterious + 2 neutral)
- Use `bash submit_all_variants.sh 1000` but only for selected variants
- Each 1000 ns job will need checkpoint-resume across wall-time boundaries

---

## 📁 Files in This Package

| File | Purpose |
|------|---------|
| `run_md_anvil.py` | Main OpenMM MD simulation script (checkpoint-safe) |
| `submit_variant.sh` | SLURM job script for single variant |
| `submit_all_variants.sh` | Submit all 10 variants in parallel |
| `setup_env.sh` | One-time conda environment setup |
| `check_status.sh` | Monitor job progress |
| `ANVIL_README.md` | This file |

---

## 🚀 Step-by-Step Instructions

### Step 1: SSH into Anvil
```bash
ssh YOUR_USER@anvil.rcac.purdue.edu
```

### Step 2: Create project directory on Anvil
```bash
mkdir -p ~/zfp36l2_md/{scripts,inputs,outputs,cleaned_pdbs,logs}
```

### Step 3: Upload files from your Mac to Anvil

**Upload scripts:**
```bash
# From your Mac terminal (NOT on Anvil):
scp -r "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/"* YOUR_USER@anvil.rcac.purdue.edu:~/zfp36l2_md/scripts/
```

**Upload input PDB files:**
```bash
scp "/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated/inputs/variant_models_run1_model0/"*.pdb YOUR_USER@anvil.rcac.purdue.edu:~/zfp36l2_md/inputs/
```

**Upload the cleaning script (needed for AF3 PDB prep):**
```bash
scp "/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated/templates/zfp36l2_ubuntu_prep_run/clean_af3_structures.py" YOUR_USER@anvil.rcac.purdue.edu:~/zfp36l2_md/scripts/
```

### Step 4: Setup environment on Anvil (one-time)
```bash
# On Anvil:
cd ~/zfp36l2_md
chmod +x scripts/*.sh
bash scripts/setup_env.sh
```

### Step 5: Verify your allocation
```bash
# Check ACCESS allocation status
mybalance
# or
accounts
```

### Step 6: Update the SLURM account name
Edit `scripts/submit_variant.sh` and verify the `#SBATCH --account=` line matches your allocation:
```bash
nano ~/zfp36l2_md/scripts/submit_variant.sh
# Look for: #SBATCH --account=BIO240054
# Change BIO240054 to your actual allocation account if different
```

### Step 7: Submit all jobs
```bash
cd ~/zfp36l2_md
bash scripts/submit_all_variants.sh 100    # 100 ns per variant (recommended)
# OR
bash scripts/submit_all_variants.sh 200    # 200 ns per variant
```

### Step 8: Monitor progress
```bash
# Check SLURM queue
squeue -u $USER

# Detailed status report
bash ~/zfp36l2_md/scripts/check_status.sh

# Watch a specific job's log
tail -f ~/zfp36l2_md/outputs/WT/md.log

# Check GPU utilization (if job is running on your node)
srun --jobid=<JOBID> nvidia-smi
```

### Step 9: If a job hits the 48h wall time
The script automatically checkpoints. Simply resubmit:
```bash
# The --resume flag is auto-detected from checkpoint files
sbatch scripts/submit_variant.sh WT 100
```

### Step 10: Download results to your Mac
```bash
# From your Mac:
scp -r YOUR_USER@anvil.rcac.purdue.edu:~/zfp36l2_md/outputs/ "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md_results/"
```

---

## ⚙️ Simulation Parameters

### Changed from GCP runs:

| Parameter | Old (GCP) | New (Anvil) | Reason |
|-----------|-----------|-------------|--------|
| Temperature | 300 K | **310 K** | Physiological temperature |
| Production | 0.23 ns | **100+ ns** | Much longer sampling |
| NVT equil | 0.01 ns | 0.1 ns | Proper heating |
| NPT equil | 0.02 ns | 0.5 ns | Proper density equilibration |
| NPT settling | — | 2 ns | Ramped restraint removal |
| Traj stride | 2 ps | 100 ps | Manageable file sizes |
| Log stride | 1 ps | 10 ps | Reasonable log size |
| Checkpoint | none | every 15 min | Crash safety |

### Kept the same:
- Force field: AMBER14 protein + RNA OL3 + TIP3P
- Integrator: LangevinMiddle, 2 fs timestep
- PME electrostatics, 1.2 nm cutoff
- HBonds constraints
- Zn²⁺ harmonic restraints (ramped: 100 → 500 → 1000 kJ/mol/nm²)
- 1.0 nm solvation padding, 0.15 M NaCl

---

## 🖥️ Anvil A100 GPU Specs & Notes

### Hardware:
- **GPU:** NVIDIA A100 40GB or 80GB HBM2e
- **GPU Memory:** Sufficient for systems up to ~500k atoms
- **Node:** 128 CPU cores, 256 GB RAM per GPU node
- **Interconnect:** HDR InfiniBand (100 Gbps)

### Important Anvil-specific notes:
1. **48-hour wall time limit** — Scripts have auto-checkpoint/resume
2. **GPU partition:** Use `--partition=gpu` in SLURM
3. **Account required:** Must specify `--account=<YOUR_ALLOCATION>`
4. **Module system:** Load `modtree/gpu`, `anaconda`, `cuda`
5. **Mixed precision:** Using mixed (single+double) for best A100 performance
6. **File system:** Use `$HOME` for scripts, `/scratch` for large trajectories if needed
7. **No preemption:** GPU jobs on Anvil are not preempted (unlike some HPC centers)
8. **Job scheduling:** Jobs typically start within minutes to hours depending on queue

### Potential issues:
- **CUDA driver mismatch:** The `setup_env.sh` installs compatible CUDA toolkit via conda
- **Memory:** ~15k atom system fits easily in A100's 40GB VRAM
- **Disk space:** 100 ns trajectory at 100 ps stride = ~100 MB DCD per variant
- **Quota:** Check `myquota` for disk limits

---

## 📊 Expected Output Per Variant

```
outputs/<VARIANT>/
├── md.log              # Energy, temperature, density log (all stages)
├── production.dcd      # Production trajectory (100 ps stride)
├── production.chk      # Production checkpoint (for resume)
├── equilibrated.chk    # Equilibration state checkpoint
├── equilibrated.pdb    # Post-equilibration structure
├── final_frame.pdb     # Last production snapshot
└── progress.txt        # Machine-readable progress tracker
```

---

## 🔧 Troubleshooting

### "CUDA platform not found"
```bash
module load modtree/gpu cuda
conda activate openmm_env
python -c "from openmm import Platform; print(Platform.getPlatformByName('CUDA'))"
```

### "Account not valid"
```bash
accounts     # List your valid accounts
# Update #SBATCH --account= in submit_variant.sh
```

### Job killed at wall time
- This is NORMAL. The script saves checkpoint 30 min before wall time.
- Just resubmit: `sbatch scripts/submit_variant.sh <VARIANT> <NS>`
- Resume is automatic.

### Out of disk space
```bash
myquota
# Consider moving DCD files to /scratch or compressing
```

### NaN coordinates / simulation explodes
- Usually means the PDB structure has clashes
- The cleaning script handles most issues
- Try re-running with the `--padding 1.5` flag for more solvation

---

## 📧 Quick Command Reference

```bash
# SSH to Anvil
ssh YOUR_USER@anvil.rcac.purdue.edu

# Submit all
cd ~/zfp36l2_md && bash scripts/submit_all_variants.sh 100

# Check queue
squeue -u $USER

# Check progress
bash scripts/check_status.sh

# Cancel all jobs
scancel -u $USER

# Cancel specific job
scancel <JOBID>

# Check allocation balance
mybalance

# Download results (from Mac)
scp -r YOUR_USER@anvil.rcac.purdue.edu:~/zfp36l2_md/outputs/ ~/Desktop/anvil_results/
```
