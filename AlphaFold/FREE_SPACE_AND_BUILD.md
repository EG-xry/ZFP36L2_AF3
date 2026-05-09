# Free Disk Space and Install Build Tools

## Current Status
- **Disk Usage**: 16G/17G (100% full, ~51M free)
- **Need**: ~1GB free for build
- **Conda pkgs cache**: 3.1GB (can be cleaned)

## Step 1: Free Up Disk Space

### A. Clean Conda Package Cache (Can free ~2-3GB)

The conda package cache stores downloaded packages. You can safely remove unused ones:

```bash
export PATH="$HOME/miniconda3/bin:$PATH"
source $HOME/miniconda3/etc/profile.d/conda.sh
conda clean --packages -y
```

**Note**: This will require re-downloading packages if you reinstall them, but it's safe.

### B. System Cleanup (Requires sudo - Run these commands)

Open a terminal and run these commands one by one:

```bash
# 1. Clean APT package cache (~100-200MB)
sudo apt clean
sudo apt autoclean

# 2. Remove unused packages (~200-500MB)
sudo apt autoremove -y

# 3. Clean system logs (keep last 7 days)
sudo journalctl --vacuum-time=7d

# 4. Remove old log files (older than 30 days)
sudo find /var/log -type f -name "*.log" -mtime +30 -delete

# 5. Check for old kernels and remove them
dpkg --list | grep linux-image
sudo apt autoremove --purge
```

### C. Check Disk Space After Cleanup

```bash
df -h /
```

**Target**: Need at least 1GB free (ideally 1.5GB)

## Step 2: Install Build Tools

Once you have enough disk space, install build tools:

```bash
# Update package list
sudo apt update

# Install build essentials (includes g++, make, etc.)
sudo apt install -y build-essential g++ make

# Verify installation
g++ --version
make --version
```

## Step 3: Verify Ready for Build

After cleanup and installation, verify:

```bash
# Check disk space (need >1GB)
df -h /

# Check compiler
g++ --version

# Check make
make --version

# Check cmake (already installed)
cmake --version
```

## Quick All-in-One Script

If you want to run everything at once (copy and paste into terminal):

```bash
#!/bin/bash
echo "=== Step 1: Clean Conda Cache ==="
export PATH="$HOME/miniconda3/bin:$PATH"
source $HOME/miniconda3/etc/profile.d/conda.sh
conda clean --packages -y

echo "=== Step 2: System Cleanup (requires sudo) ==="
sudo apt clean
sudo apt autoclean
sudo apt autoremove -y
sudo journalctl --vacuum-time=7d

echo "=== Step 3: Check Disk Space ==="
df -h /

echo "=== Step 4: Install Build Tools (requires sudo) ==="
sudo apt update
sudo apt install -y build-essential g++ make

echo "=== Step 5: Verify ==="
g++ --version
make --version
```

## Expected Results

After cleanup:
- **Disk space**: Should have 1-2GB free
- **Build tools**: g++, make installed and working
- **Ready**: Can proceed with OpenMM-Plumed-MPI build

## Troubleshooting

### If still low on disk space:
1. Check for large files: `find ~ -type f -size +100M 2>/dev/null`
2. Remove unused conda environments: `conda env list` then `conda env remove -n <env_name>`
3. Check Docker (if installed): `docker system prune -a --volumes`

### If build tools installation fails:
- Check internet connection
- Try: `sudo apt update && sudo apt install -y build-essential`

