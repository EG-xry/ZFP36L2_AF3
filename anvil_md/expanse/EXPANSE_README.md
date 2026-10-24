# ZFP36L2 MD on SDSC Expanse (V100 GPU) — Quick Start

## Overview

Run the same ZFP36L2 MD simulations on SDSC Expanse V100 GPUs while waiting for Anvil provisioning.

**Username:** `YOUR_USER`  
**Budget:** 477 GPU-hours remaining  
**GPU:** NVIDIA V100 32GB (slower than A100 but fully functional)  
**Account:** `bio260224`  
**Partition:** `gpu-shared` (shared V100 nodes)

---

## Performance: V100 vs A100

| Metric | V100 (Expanse) | A100 (Anvil) |
|--------|---------------|--------------|
| Speed for ~15k atoms | ~100–200 ns/day | ~200–400 ns/day |
| 50 ns per variant | ~6–12 hours | ~3–6 hours |
| 10 variants × 50 ns | ~60–120 GPU-hours | ~30–60 GPU-hours |

**Recommendation:** Run 50 ns per variant on Expanse (~100 GPU-hours total, fits your 477h budget). Extend on Anvil later using checkpoint resume.

---

## Step-by-Step

### Step 1: SSH into Expanse
```bash
ssh YOUR_USER@login.expanse.sdsc.edu
```

### Step 2: Create project directory
```bash
mkdir -p ~/zfp36l2_md/{scripts,inputs,outputs,cleaned_pdbs,logs}
```

### Step 3: Upload files from your Mac

```bash
# Upload scripts (run from Mac terminal)
scp "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/expanse/"*.sh YOUR_USER@login.expanse.sdsc.edu:~/zfp36l2_md/scripts/

# Upload the main MD Python script (same one used for Anvil)
scp "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/run_md_anvil.py" YOUR_USER@login.expanse.sdsc.edu:~/zfp36l2_md/scripts/

# Upload the PDB cleaning script
scp "/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated/templates/zfp36l2_ubuntu_prep_run/clean_af3_structures.py" YOUR_USER@login.expanse.sdsc.edu:~/zfp36l2_md/scripts/

# Upload best-model PDB files
scp "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/best_models/"*.pdb YOUR_USER@login.expanse.sdsc.edu:~/zfp36l2_md/inputs/

# Also upload the run1 fallback PDBs (already on your Mac)
scp "/Users/eric/Desktop/ZFP36L2 AF3/MD_GCP_Pipeline_Consolidated/inputs/variant_models_run1_model0/"*.pdb YOUR_USER@login.expanse.sdsc.edu:~/zfp36l2_md/inputs/

# Upload the status checker
scp "/Users/eric/Desktop/ZFP36L2 AF3/anvil_md/check_status.sh" YOUR_USER@login.expanse.sdsc.edu:~/zfp36l2_md/scripts/
```

### Step 4: Setup environment (one-time, on Expanse)
```bash
cd ~/zfp36l2_md
chmod +x scripts/*.sh
bash scripts/setup_env_expanse.sh
```

### Step 5: Verify allocation
```bash
expanse-client user
# or
sacctmgr show associations where user=YOUR_USER format=Account,Partition,QOS --parsable
```

### Step 6: Test with a quick GPU job
```bash
sbatch -A bio260224 -p gpu-shared --gpus=1 --nodes=1 --ntasks-per-node=1 --time=00:05:00 --wrap="nvidia-smi && hostname" --output=test_gpu.out
squeue -u $USER
cat test_gpu.out
```

### Step 7: Submit all variants (50 ns each)
```bash
cd ~/zfp36l2_md
bash scripts/submit_all_expanse.sh 50
```

### Step 8: Monitor
```bash
# Queue status
squeue -u $USER

# Job progress
cat ~/zfp36l2_md/outputs/*/progress.txt

# Watch a log
tail -f ~/zfp36l2_md/outputs/WT/md.log
```

### Step 9: Download results (from Mac)
```bash
scp -r YOUR_USER@login.expanse.sdsc.edu:~/zfp36l2_md/outputs/ "/Users/eric/Desktop/ZFP36L2 AF3/expanse_md_results/"
```

---

## Key Differences from Anvil

| Setting | Anvil | Expanse |
|---------|-------|---------|
| Account | `bio260224-gpu` | `bio260224` |
| Partition | `gpu` | `gpu-shared` |
| GPU request | `--gpus-per-node=1` | `--gpus=1` |
| GPU type | A100 (40/80GB) | V100 (32GB) |
| Module loads | `modtree/gpu`, `anaconda`, `cuda` | `gpu`, `slurm`, `anaconda3` |
| Username | `YOUR_USER` | `YOUR_USER` |
| Max wall time | 48h | 48h |

The Python simulation script (`run_md_anvil.py`) is **identical** — no changes needed.

---

## Transitioning to Anvil Later

Once Anvil is provisioned, you can continue from Expanse checkpoints:

1. Download Expanse outputs to your Mac
2. Upload to Anvil
3. Re-run with higher `--production-ns` — checkpoint resume will continue from where Expanse left off

```bash
# Download from Expanse (Mac)
scp -r YOUR_USER@login.expanse.sdsc.edu:~/zfp36l2_md/outputs/ /tmp/expanse_outputs/

# Upload to Anvil (Mac)  
scp -r /tmp/expanse_outputs/* YOUR_USER@anvil.rcac.purdue.edu:~/zfp36l2_md/outputs/

# On Anvil, extend to 200 ns (continues from 50 ns checkpoint)
bash scripts/submit_all_variants.sh 200
```

---

## Quick Command Reference

```bash
# SSH to Expanse
ssh YOUR_USER@login.expanse.sdsc.edu

# Submit all 50ns
cd ~/zfp36l2_md && bash scripts/submit_all_expanse.sh 50

# Check queue
squeue -u $USER

# Cancel all
scancel -u $USER

# Check balance
expanse-client user
```
