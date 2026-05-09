import os
import shutil
import csv

# Set home directory
home = "/Users/eric/Desktop/AlphaFold-MetaInference-Tutorial-main"
os.chdir(home + '/AlphaFold-IDP/prep_run')

dir = os.getcwd()
print(f"Working in: {dir}")

# TDP-43 sequence
fasta_sequence = "MSEYIRVTEDENDEPIEIPSEDDGTVLLSTVTAQFPGACGLRYRNPVSQCMRGVRLVEGILHAPDAGAGNLVYVVNYPKDNKRKMDETDASSAVKVKRAVQKTSDLIVLGLPAKTTEQDLKEYFSTFGEVLMVQVKKDLKTGHSKGFGFVRFTEYETQVKVMSQRHMIDGRACDCKLPNSKQSQDEPLRSRKVFVGRCTEDMTEDELREFFSQYGDVMDVFIPKPFRAFAFVTFADDQIAQSLCGEDLIIKGISVHISNAEPKHNSNRQLERSGRFGGNPGGFGNQGGFGNSRGGGAGLGNNQGSNMGGGMNFGAFSINPAMMAAAQAALQSSAGMMGMLASQQNQSGPSGNNQNQGNMQREPNQAFGSGNNSYSGSNSGAAIGAGSASNAGSGSGFNGGFGSSMDSKSSGAGM"

# Conditions
pH = 7.4
temp = 298
ionic = 0.2
PAE_cut = 5
Pr_cut = 0.2
NR = 2

# Decide the plddt based ordered (od) and disordered (dd) regions
ordered_domains = {'od1': [3, 79], 'od2': [104, 178], 'od3': [191, 260]}
disordered_domains = {'dd1': [1, 2], 'dd2': [80, 103], 'dd3': [179, 190], 'dd4': [261, 414]}

# AF input created in AF-distance map prediction and used in AF-MI
pdb_af = 'tdp43_WtoA_bf4cc_unrelaxed_rank_001_alphafold2_ptm_model_3_seed_000.pdb'
json_af = 'tdp43_WtoA_bf4cc_predicted_aligned_error_v1.json'
npy_af = 'alphafold2_ptm_model_3_seed_000_prob_distributions.npy'
mean_af = 'alphafold2_ptm_model_3_seed_000_mean.csv'

# Copy AF data
print("Copying AF data files...")
shutil.copy2(dir + "/../AF_DATA/" + pdb_af, dir + "/pdb_af.pdb")
shutil.copy2(dir + "/../AF_DATA/" + json_af, dir + "/pae.json")
shutil.copy2(dir + "/../AF_DATA/tdp43_WtoA_bf4cc_distmat/" + mean_af, dir + "/mean_af.csv")
shutil.copy2(dir + "/../AF_DATA/tdp43_WtoA_bf4cc_distmat/" + npy_af, dir + "/prob.npy")
print("AF data files copied successfully")

# Write sequence file
f = open("sequence.dat", "w")
f.write(fasta_sequence)
f.close()
print("Sequence file created")

# Write the csv files
with open('ordered_domains.csv', 'w') as f:
    w = csv.DictWriter(f, ordered_domains.keys())
    w.writeheader()
    w.writerow(ordered_domains)

with open('disordered_domains.csv', 'w') as f:
    w = csv.DictWriter(f, disordered_domains.keys())
    w.writeheader()
    w.writerow(disordered_domains)
print("Domain CSV files created")

# Copy OPENMM calvados forcefield files
shutil.copy2(dir + "/../scripts_prep/gen_xml_and_constraints.py", dir)
shutil.copy2(dir + "/../scripts_prep/residues.csv", dir)
print("Forcefield files copied")

# Generate XML and constraints
path_gen_xml = f'{dir}/gen_xml_and_constraints.py sequence.dat {pH} {temp} {ionic}'
print(f"Running: python {path_gen_xml}")
os.system(f'python {path_gen_xml}')

print("\nSetup complete!")

