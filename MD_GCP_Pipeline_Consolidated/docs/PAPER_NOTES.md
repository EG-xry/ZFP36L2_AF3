# Paper Notes

## Local paper found and ingested

- File: `inputs/reference_paper/Joshua_Project_Silvia.pdf`
- Title: *Identification of deleterious non-synonymous single nucleotide polymorphisms in the mRNA decay activator ZFP36L2*
- Journal: RNA Biology (published online Dec 13, 2024)
- DOI shown in PDF: `10.1080/15476286.2024.2437590`
- Extracted text snippets saved to:
  - `docs/upstream_notes/Joshua_Project_Silvia_MD_extract.txt`

## MD protocol extracted from local paper (Methods text)

From the local PDF molecular dynamics methods section:

- Engine: **GROMACS 2020.3**
- Force field: **AMBER99SB**
- Solvent: **TIP3P** water
- Box: cubic, **10 Å minimum** protein-box edge distance
- Ions: neutralized with **0.15 M NaCl**
- Thermodynamic state: **310 K**, **1 bar**
- Ensembles: **NVT** and **NPT**
- Equilibration: **10 ns**, **2 fs** timestep
- Constraints: **LINCS** (hydrogen bond constraints)
- Electrostatics: **PME**, **12 Å cut-off**
- Production: **1000 ns** per system
- Replicates: **3 independent runs per variant**

## Variants reported in paper MD section

The MD methods paragraph lists WT and variants:

- Y154H
- R160W
- R180C (likely typo for R184C, based on the rest of the paper)
- C206F
- G204D

## Mapping to this project

Current cloud-ready manifest in this workspace is intentionally broader than the paper set:

- Includes controls and additional variants from your project set
- Excludes per your instruction: **R160E** and **H212A**
- Keeps **C208A** for analysis

## Note on AF-MI paper

There is also a separate AF-MI/Nature Communications reference used for method inspiration in earlier notes. That is a different workflow lineage from this RNA Biology paper and should be treated as supplemental methodology context, not the same protocol.
