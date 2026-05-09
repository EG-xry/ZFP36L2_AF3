# Verified Commands to Restart PLUMED Build

Run these commands **one at a time** in your terminal:

## Step 1: Navigate and activate conda
```bash
cd /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build/plumed2-2.8.2
source /Users/eric/miniconda3/etc/profile.d/conda.sh
conda activate afmi-tutorial
```

## Step 2: Verify MPI compilers are available
```bash
which mpicc mpicxx
mpicc --version | head -2
```

If these commands fail, install MPICH:
```bash
conda install -y -c conda-forge mpich
```

## Step 3: Set compiler environment variables
```bash
export CC=$(which mpicc)
export CXX=$(which mpicxx)
export PATH=/Users/eric/miniconda3/envs/afmi-tutorial/bin:$PATH
```

## Step 4: Clean previous build (if exists)
```bash
if [ -f Makefile ]; then make distclean > /dev/null 2>&1; fi
```

## Step 5: Configure PLUMED
```bash
./configure --prefix=/Users/eric/miniconda3/envs/afmi-tutorial --enable-mpi 2>&1 | tee configure_now.log
```

Wait for this to complete, then verify:
```bash
test -f config.status && echo "Config OK" || echo "Config FAILED"
```

## Step 6: Start the build
```bash
nohup make -j4 > build_running.log 2>&1 &
echo "Build PID: $!"
```

## Step 7: Verify build is running
```bash
sleep 5
ps aux | grep "[m]ake -j4"
tail -20 build_running.log
```

## Monitor progress
```bash
# Watch the log in real-time
tail -f /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/plumed_build/plumed2-2.8.2/build_running.log

# Or use the monitor script
bash /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/MONITOR_BUILD.sh
```

---

**OR** run the automated script:
```bash
bash /Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main/RUN_THESE_COMMANDS.sh
```

