#!/usr/bin/env python
"""
Script to run through the AF-IDP_colab.ipynb notebook steps locally
Adapts Colab notebook for local execution
"""
import os
import sys
import shutil
import subprocess
import csv

# Set library path for OpenMM
os.environ['DYLD_LIBRARY_PATH'] = '/Users/eric/micromamba/envs/omm/lib'

# Base directory
home = os.getcwd()
print(f"Working directory: {home}")

# Step 1: Setup prep_run directory
print("\n=== Step 1: Setting up prep_run directory ===")
prep_run_dir = os.path.join(home, 'AlphaFold-IDP/prep_run')
os.chdir(prep_run_dir)
print(f"Changed to: {prep_run_dir}")

# TDP-43 sequence and parameters (from notebook)
fasta_sequence = "MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVSQCMRGVRLVEGILHAPDAGAGNLVYVVNYPKDNKRKMDETDASSAVKVKRAVQKTSDLIVLGLPAKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRACDCKLPNSKQSQDEPLRSRKVFVGRCTEDMTEDELREFFSQYGDVMDVFIPKPFRAFAFVTFADDQIAQSLCGEDLIIKGISVHISNAEPKHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSAGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGAGSASNAGSGSGFNGGFGSSMDSKSSGAGM"
pH = 7.4
temp = 298
ionic = 0.2
PAE_cut = 5
Pr_cut = 0.2
NR = 2

ordered_domains = {'od1': [3, 79], 'od2': [104,178],'od3':[191,260]}
disordered_domains={'dd1': [1,2],'dd2':[80,103],'dd3':[179,190],'dd4':[261,414]}

# AF input files
pdb_af='tdp43_WtoA_bf4cc_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
json_af='tdp43_WtoA_bf4cc_predicted_aligned_error_v1.json'
npy_af='alphafold2_ptm_model_3_seed_000_prob_distributions.npy'
mean_af='alphafold2_ptm_model_3_seed_000_mean.csv'

# Copy AF data (if not already present)
af_data_dir = os.path.join(home, 'AlphaFold-IDP/AF_DATA')
if os.path.exists(af_data_dir):
    if not os.path.exists('pdb_af.pdb'):
        shutil.copy2(os.path.join(af_data_dir, pdb_af), 'pdb_af.pdb')
        print(f"Copied {pdb_af}")
    if not os.path.exists('pae.json'):
        shutil.copy2(os.path.join(af_data_dir, json_af), 'pae.json')
        print(f"Copied {json_af}")
    if not os.path.exists('mean_af.csv'):
        shutil.copy2(os.path.join(af_data_dir, 'tdp43_WtoA_bf4cc_distmat', mean_af), 'mean_af.csv')
        print(f"Copied {mean_af}")
    if not os.path.exists('prob.npy'):
        shutil.copy2(os.path.join(af_data_dir, 'tdp43_WtoA_bf4cc_distmat', npy_af), 'prob.npy')
        print(f"Copied {npy_af}")

# Write sequence.dat
with open('sequence.dat', 'w') as f:
    f.write(fasta_sequence)

# Write CSV files
with open('ordered_domains.csv', 'w') as f:
    w = csv.DictWriter(f, ordered_domains.keys())
    w.writeheader()
    w.writerow(ordered_domains)

with open('disordered_domains.csv', 'w') as f:
    w = csv.DictWriter(f, disordered_domains.keys())
    w.writeheader()
    w.writerow(disordered_domains)

# Copy scripts if needed
scripts_prep_dir = os.path.join(home, 'AlphaFold-IDP/scripts_prep')
if os.path.exists(scripts_prep_dir):
    if not os.path.exists('gen_xml_and_constraints.py'):
        shutil.copy2(os.path.join(scripts_prep_dir, 'gen_xml_and_constraints.py'), '.')
    if not os.path.exists('residues.csv'):
        shutil.copy2(os.path.join(scripts_prep_dir, 'residues.csv'), '.')

# Step 2: Generate forcefield XML
print("\n=== Step 2: Generating forcefield XML ===")
if os.path.exists('gen_xml_and_constraints.py'):
    cmd = f'python gen_xml_and_constraints.py sequence.dat {pH} {temp} {ionic}'
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: {result.stderr}")
        sys.exit(1)
    print("Forcefield XML generated")
else:
    print("WARNING: gen_xml_and_constraints.py not found")

# Step 3: Generate plumed files
print("\n=== Step 3: Generating plumed files ===")
if os.path.exists(scripts_prep_dir):
    if not os.path.exists('make_plumed_distmat.py'):
        shutil.copy2(os.path.join(scripts_prep_dir, 'make_plumed_distmat.py'), '.')
    if not os.path.exists('make_plumed_analysis.py'):
        shutil.copy2(os.path.join(scripts_prep_dir, 'make_plumed_analysis.py'), '.')

if os.path.exists('make_plumed_distmat.py'):
    cmd = f'python make_plumed_distmat.py sequence.dat {PAE_cut} {Pr_cut}'
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"WARNING: {result.stderr}")
    else:
        print("plumed.dat generated")

if os.path.exists('make_plumed_analysis.py'):
    cmd = f'python make_plumed_analysis.py sequence.dat'
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"WARNING: {result.stderr}")
    else:
        print("plumed_analysis.dat generated")

# Use TDP-43 specific plumed files if available
if os.path.exists(scripts_prep_dir):
    tdp_plumed = os.path.join(scripts_prep_dir, 'plumed_TDP-43.dat')
    tdp_analysis = os.path.join(scripts_prep_dir, 'plumed_analysis_TDP-43.dat')
    if os.path.exists(tdp_plumed):
        shutil.copy2(tdp_plumed, 'plumed.dat')
        print("Using TDP-43 specific plumed.dat")
    if os.path.exists(tdp_analysis):
        shutil.copy2(tdp_analysis, 'plumed_analysis.dat')
        print("Using TDP-43 specific plumed_analysis.dat")

# Step 4: Energy minimization
print("\n=== Step 4: Running energy minimization ===")
if os.path.exists('simulate_em.py'):
    cmd = f'python simulate_em.py {pH} {temp}'
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR: {result.stderr}")
        print(f"STDOUT: {result.stdout}")
        sys.exit(1)
    print("Energy minimization completed")
    if os.path.exists('checkpoint'):
        print("Checkpoint file created")
else:
    print("ERROR: simulate_em.py not found")
    sys.exit(1)

print("\n=== Setup complete! ===")
print("Next steps:")
print(f"1. Run AF-MI simulation: mpirun -np {NR} python simulate.py {pH} {temp}")
print("2. Run analysis scripts")

