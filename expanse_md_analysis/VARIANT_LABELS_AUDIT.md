# Audit of variant labels vs. actual mutations in the AF3 structures

**Generated:** 2026-05-27 from `expanse_md_results/<variant>/equilibrated.pdb`,
cross-checked against `Data/<variant>/.../*_job_request.json` and the canonical
sequence `Protein_Sequence.txt`.

The ZFD construct used is residues **151–220** of human ZFP36L2 (UniProt P47974),
sequence `STRYKTELCRPFEESGTCKYGEKCQFAHGFHELRSLTRHPKYKTELCRTFHTIGFCPYGPRCHFIHNADE`.
AF3 outputs use construct-local numbering 1..70 → ZFP36L2 numbering = PDB + 150.

## Verdict

| Folder label | Actual mutation (ZFP36L2 numbering) | Match? | Notes |
|---|---|---|---|
| C206F | **C206F** (C→F at PDB-pos 56 = Z 206) | ✓ | OK |
| **C208A** | **C212A** (C→A at PDB-pos 62 = Z 212) | ✗ | **MISLABELED** — see below |
| **D210E** | **D219E** (D→E at PDB-pos 69 = Z 219) | ✗ | **MISLABELED** — see below |
| G204D | G204D (G→D at PDB-pos 54 = Z 204) | ✓ | OK |
| R160W | R160W (R→W at PDB-pos 10 = Z 160) | ✓ | OK |
| R184C | R184C (R→C at PDB-pos 34 = Z 184) | ✓ | OK |
| S185N | S185N (S→N at PDB-pos 35 = Z 185) | ✓ | OK |
| T187S | T187S (T→S at PDB-pos 37 = Z 187) | ✓ | OK |
| Y154H | Y154H (Y→H at PDB-pos 4 = Z 154) | ✓ | OK |

So your intuition is correct: **"C208A" is actually C212A**.

## Where the mistake came from

The original AF3 job-request JSONs (see
`Data/C208A/fold_zfp36l2_deleterious_c208a/fold_zfp36l2_deleterious_c208a_job_request.json`)
contain the *correct* sequence (mutation at construct position 62, which is C212
in UniProt numbering). The mislabeling is purely in the folder name and in the
downstream pipeline scripts (`anvil_md/submit_all_variants.sh`,
`anvil_md/select_best_models.py`, etc.).

I checked all 9 mutant JSONs:
- 7 of them have the correct mutation at the position implied by their folder name
- **C208A**: actual mutation is **C212A** (C212 is one of the two CCCH Cys-ligands
  of zinc finger 2 — mutating it destroys Zn coordination of ZF2, which is
  *exactly* the right "knockout control" everyone wants for this kind of study)
- **D210E**: actual mutation is **D219E** (D219 is the second-to-last residue
  of the construct, in the disordered C-terminal tail — not a known ARE-binding
  residue. This is probably just a copy/paste error in whoever set up the job.)

## Biological implications

| Label | Real mutation | Why it matters |
|---|---|---|
| C208A → **C212A** | C212 = the third C of the CCCH zinc-coordination motif of ZF2 (C197, C206, C212, H214) | Loses one Zn²⁺ coordinator → ZF2 unfolds → RNA binding abolished. This is a *positive control* for "deleterious," not an SNP from dbSNP. Should be labeled and discussed as a control. |
| D210E → **D219E** | D219 is in the disordered C-terminal tail, away from ZF2 | Probably truly neutral biologically (conservative D→E swap in a flexible region). The "neutral" label is plausibly correct, but **the residue you intended (D210) is not actually mutated** in this run. D210 in ZFP36L2 sits between H214 (ZF2 Zn-coord histidine) and the C-tail, so a real D210E might have *different* effects than what we measured. |

## Practical consequences for the current analysis

- In `RAMOS_COMPARISON.md` and `REPORT.md`, every mention of "C208A" should
  be read as **C212A** and every "D210E" as **D219E**.
- The MD trajectories themselves are valid — the protein in the simulation
  is unambiguously C212A or D219E.
- **C212A behaves like a strong-deleterious variant in our data** (R160
  H-bonds drop from 0.97 → 0.89, ZF2 RMSD stays low, but the contact
  count drops by ~36) which is exactly what you'd expect if you mutate a
  Zn-coordinating Cys. So this run is still informative — just label it
  correctly as a "ZF2 Zn-knockout control" not as the dbSNP C208A.
- **D219E being mislabeled is more annoying**: we have a useless data point
  for "is D210E neutral?" The dbSNP variant the team meant to study is
  ZFP36L2:D210E (which corresponds to PDB construct position 60, residue P
  in our WT sequence — wait, the construct doesn't even *have* a D at
  position 60. There may be no dbSNP variant called D210E in this protein,
  in which case the correct call is to drop this row from the "neutral"
  set entirely, OR rebuild the AF3 input with the correct mutation if
  D210E was intended in some alternative reference sequence.)

## Recommended actions

1. **Rename folders** to reflect the actual mutation:
   - `expanse_md_results/C208A/` → `expanse_md_results/C212A/`
   - `expanse_md_results/D210E/` → `expanse_md_results/D219E/`
   Update `anvil_md/select_best_models.py`, `anvil_md/submit_*.sh`,
   `expanse_md_analysis/*.py`, and any presentation/manifests files.
2. **Reclassify C212A** from "deleterious dbSNP variant" to
   "Zn²⁺-coordination knockout control" in all writeups and Block_Design.txt.
3. **Investigate D210E**: is there a real dbSNP variant of that name? Check
   the original spreadsheet of 5 deleterious + 3 neutrals from the
   bioinformatic pipeline. If `D219E` was intended, fine — relabel. If
   `D210E` was intended (real dbSNP), the AF3 input was generated against
   the wrong residue and that run should be discarded and re-generated.
4. **Decide whether to keep D219E in the "neutral" group** for analysis or
   to drop it. Recommendation: drop it for now; it isn't a published
   dbSNP variant anyway and contaminates the neutral statistics.

After (1)–(4), regenerate Figures 1 and 2 — they will both stay essentially
the same but with corrected x-axis labels.
