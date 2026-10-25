# ZFP36L2 variants: AlphaFold 3 + molecular dynamics

October 2025 - July 2026  

Can molecular dynamics run on AlphaFold 3 models tell harmful ZFP36L2 variants from harmless ones? This repository holds the write up, data and code for a computational project that put that question to eight experimentally labeled missense variants in the ZFP36L2 tandem zinc-finger domain.

**Status: paused (July 2026).** Neither AlphaFold 3 confidence nor MD separates the deleterious variants from the neutral ones. The MD runs fall into two sets: runs that carry each variant, and 19 runs of the wild type protein that test how consistent independent simulations of one complex are. That consistency test shaped the conclusion most: separate WT runs disagree with one another about as much as the variants differ.

![a, Pipeline: ten constructs, AlphaFold 3 prediction, all-atom MD, 45 readouts, AUROC against EMSA labels. b, Deviation of each construct from the wild type in AF3 ipTM, PAE and pLDDT, and their mean. c, AF3 pair-ipTM for a non-binding RNA, the wild type and the eight variants. d, T156 hydrogen bonds in 19 wild-type runs. e, Variant means against the range of the wild-type runs. f, Every exact p-value that eight variants allow, against the Bonferroni line.](figures/fig0_pipeline.png)

**Pipeline and summary.** **a**, Ten constructs of the ZFP36L2 tandem zinc fingers (residues 151-220) were predicted with AlphaFold 3 on the ARE 5′-UUAUUUAUU-3′, simulated with all-atom MD, and 45 readouts were scored against EMSA labels as AUROC; the ROC curve shown is AF3 pair-ipTM. **b**, Each construct's 40 AF3 models against the wild type's, in wild-type standard deviations (+ = worse; asterisks, Welch p < 0.05), and the mean of the three metrics (± s.e.m.). AF3 responds to the mutations, but not along the EMSA labels: the best single metric and the mean both reach an AUROC of 0.73. **c**, AF3 interface confidence falls for a non-binding RNA but barely moves between variants. **d**, T156-RNA hydrogen bonds in 19 runs of the wild-type complex: the two halves of a run mostly agree, separate runs do not. **e**, Variant means from variant runs (dots) against the 19 wild-type runs (squares; shaded range). **f**, With five deleterious and three neutral variants, even a perfect ranking gives p = 0.018, above the Bonferroni threshold for 45 readouts.

## Summary

ZFP36L2 is an RNA-binding protein that sends messenger RNAs to be degraded. Two CCCH zinc fingers grip AU-rich elements (AREs) in the 3′ untranslated region of its targets, and the protein then recruits the decay machinery. Akçeşme et al. (2025) predicted which single-nucleotide variants in that domain should be damaging and tested eight of them by EMSA: five lost RNA binding, three did not. I asked whether physics could recover those labels from structure alone. I predicted each variant RNA complex with AlphaFold 3, simulated it with all-atom MD, and scored 45 structural readouts against the EMSA labels.

- **AlphaFold 3 is almost equally confident about every variant.** Interface confidence (pair-ipTM) spans 0.71 to 0.78, and its best readout separates the classes with an AUROC of 0.733
- **The MD runs form two sets.** The Expanse runs carry each variant: 13 runs across the eight labeled variants, one or two each. Every Anvil run contains the wild-type protein, because the Anvil job script resolved each input to the wild-type structure. With the runs in the wild-type folder, they form a 19-run wild-type test
- **Wild-type runs do not agree with one another.** T156-RNA hydrogen bonds land anywhere from 0.37 to 0.71 across 19 runs of the same protein; pinning one protein's value to ±0.05 would take about 13 runs. Aromatic stacking switches fully on in 6 runs and stays off in 12
- **Measured against that spread, MD does no better than AF3.** The best MD readouts tie pair-ipTM at 0.733. T156 hydrogen bonds score 0.267 on the variant runs; an earlier analysis that pooled runs by folder had given 0.800
- **The strongest lead found does not hold up.** Aromatic stacking reaches 0.933 (uncorrected p = 0.036), but it rests on two runs, fails correction for the 45 readouts tested, and the wild-type protein itself stacks in only 6 of 19 runs
- **Eight variants cap what any method can show.** Even a perfect ranking would give p = 0.018

I paused rather than build on a lead this thin. The WT test is the most transferable result here: it sets the bar any variant effect has to clear, and with one or two runs per variant that bar is high. In the future if I have time, everything needed to pick the project back up is in this repository: inputs, a record of what each MD run simulated, analysis code, and the exact numbers behind every figure.

## The system

![AlphaFold 3 model of the ZFP36L2 tandem zinc fingers bound to UUAUUUAUU RNA with two zinc ions; variant sites shown as colored spheres.](figures/fig1_structure.png)

**AlphaFold 3 starting structure.** Top-ranked AlphaFold 3 model of the wild-type ZFP36L2 tandem zinc fingers (residues 151-220) bound to the ARE 5′-UUAUUUAUU-3′ with two Zn²⁺, the model that started the wild-type MD runs. Spheres mark the Cα of each variant site, colored by EMSA label. † See the naming note under the table.


| Variant | EMSA label                        | Location                   | AF3 pair-ipTM | T156-RNA H-bonds per frame | Aromatic stacking | MD runs |
| ------- | --------------------------------- | -------------------------- | ------------- | -------------------------- | ----------------- | ------- |
| WT      | reference                         |                            | 0.765         | 0.524*                     | 0.260*            | 19*     |
| Y154H   | deleterious                       | RYKTEL lead-in to finger 1 | 0.775         | 0.616                      | 0.575             | 2       |
| R160W   | deleterious                       | finger 1                   | 0.756         | 0.519                      | 0.738             | 2       |
| R184C   | deleterious                       | linker                     | 0.758         | 0.537                      | 0.318             | 1       |
| G204D   | deleterious                       | finger 2                   | 0.767         | 0.489                      | 0.748             | 2       |
| C206F   | deleterious                       | finger 2, zinc ligand      | 0.710         | 0.533                      | 0.826             | 2       |
| S185N   | neutral                           | linker                     | 0.771         | 0.423                      | 0.781             | 1       |
| T187S   | neutral                           | linker                     | 0.762         | 0.445                      | 0.834             | 2       |
| D210E†  | neutral                           | C-terminal tail            | 0.769         | 0.587                      | 0.826             | 1       |
| C208A†  | zinc-knockout control, not scored | finger 2, zinc ligand      | 0.752         | 0.597                      | 0.199             | 1       |


AF3 values are means of 40 models (8 seeds × 5). MD values are means over each variant's own runs; aromatic stacking is the fraction of frames with a ring stacked on a base.  The WT test: 19 complete runs of the wild-type protein. † Names follow Akçeşme et al. 2025. In the construct numbering used here, the AF3 jobs mutate D219 and C212 ([naming notes](#naming-notes)).

## What I did

1. **AlphaFold 3.** I predicted the domain (residues 151-220) with the ARE 5′-UUAUUUAUU-3′ and two Zn²⁺ for 12 constructs on AlphaFold Server: 90 jobs and 450 models. → [AlphaFold 3 results](docs/01_AF3_results.md)
2. **Molecular dynamics.** I started from the top ranked model per variant and ran OpenMM 8.2 with AMBER ff14SB, RNA.OL3 and TIP3P in 0.15 M NaCl at 310 K. I launched 40 production runs of up to 200 ns on SDSC Expanse (V100) and Purdue Anvil (A100). The Expanse runs carry each variant; every Anvil run contains the wild-type protein, so those runs serve as a wild-type consistency test. Of the 33 runs that finished, 14 carry a variant (2.8 µs) and 19 are wild type (3.7 µs). → [MD results](docs/02_MD_results.md)
3. **Benchmark.** I scored 45 readouts against the EMSA labels: 7 AF3 confidence metrics, 25 MD contact and energy metrics, and 13 aromatic-stacking metrics. Each was averaged over a variant's own runs, its direction was fixed in advance, and it was scored as AUROC with an exact permutation p. I also compared each against the spread of the wild-type runs. → [Conclusions](docs/03_conclusions.md)



## Results at a glance


| Readout                                                                         | AUROC | one-sided p       | Where                                                                                  |
| ------------------------------------------------------------------------------- | ----- | ----------------- | -------------------------------------------------------------------------------------- |
| AF3 pair-ipTM (best AF3 readout)                                                | 0.733 | 0.196             | [AF3 results](docs/01_AF3_results.md)                                                  |
| Best MD contact and energy readouts: Coulomb energy, salt bridges, Y170 H-bonds | 0.733 | 0.196             | [MD results](docs/02_MD_results.md#the-screen-on-variant-runs)                         |
| T156-RNA hydrogen bonds on variant runs (0.800 when runs were pooled by folder) | 0.267 | 0.875             | [MD results](docs/02_MD_results.md#variants-against-the-wild-type-spread)              |
| Aromatic stacking on variant runs                                               | 0.933 | 0.036 uncorrected | [MD results](docs/02_MD_results.md#variants-against-the-wild-type-spread)              |
| A perfect ranking of these eight variants                                       | 1.000 | 0.018             | [Conclusions](docs/03_conclusions.md#eight-variants-cannot-carry-a-significant-result) |




## Why the project is paused

The goal was never just a number above 0.5 necessarily. It was a mechanism: a contact that deleterious variants break and neutral ones keep, visible in MD but invisible to AlphaFold 3. The WT testing showed why that was hard to find. Independent runs of one protein disagree by more than the differences I was looking for, so one or two runs per variant cannot pin anything down, and the T156 signal from the pooled analysis disappears once variant and WT runs are kept apart. To be fair, aromatic stacking on the variant runs points the way the textbook mechanism would predict. Yet it depends on two trajectories, it cannot survive correction for the readouts tested, and the wild-type protein, which binds RNA normally, stacks less often than most deleterious variants. The point I conclude then is not that MD failed. It is that eight labeled variants, with one or two runs each, cannot settle the question. More nanoseconds would not change that; more runs per variant and more EMSA validated variants would.

## Read the write-up

1. [AlphaFold 3 results](docs/01_AF3_results.md): inputs, confidence per variant, the non-binding RNA control, and two traps in AF3 summary metrics.
2. [Molecular dynamics results](docs/02_MD_results.md): protocol, the two run sets, the wild-type consistency test, the screen on variant runs, and variants against the wild-type spread.
3. [Conclusions](docs/03_conclusions.md): the n = 8 ceiling, how analysis choices move the answer, lessons, and what would justify resuming.



## Repository map


| Path                                                            | Contents                                                                                                                                                                             |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `docs/`                                                         | The three write-up pages                                                                                                                                                             |
| `figures/`                                                      | Every figure as editable SVG plus PNG; `src/` rebuilds them, `data/` holds derived CSVs, the run inventory (`md_run_inventory.csv`) and `tables.md` (the numbers behind each figure) |
| `Data/`                                                         | AlphaFold Server outputs: 12 constructs × 5-8 seeds (models, confidences, MSAs, templates)                                                                                           |
| `results/`, `analysis/`                                         | AF3 metric tables, including the non-binding RNA control (`analysis/interface_metrics_models.csv`)                                                                                   |
| `anvil_md/`                                                     | Production MD driver (`run_md_anvil.py`), structure prep, HPC submission scripts, MD starting models (`best_models/`)                                                                |
| `expanse_md_analysis/`                                          | MD analysis scripts, per-trajectory metrics, AUROC tables, and the working notes from each stage of the analysis                                                                     |
| `MD_GCP_Pipeline_Consolidated/`, `md_simulation/`, `AlphaFold/` | Early pipeline development and sub-nanosecond smoke tests; not the basis of any result                                                                                               |


The MD trajectories (about 23 GB of DCD files) are too large for GitHub and are not included. The per trajectory summaries that the figures need, and the record of what each run simulated, are included.

## Reproducing the figures

```bash
python -m venv .venv
.venv/bin/pip install numpy pandas matplotlib scipy tabulate
.venv/bin/python figures/src/af3_metrics.py      # per-model AF3 metrics from Data/
.venv/bin/python figures/src/make_figures.py     # every data figure + figures/data/tables.md
```

Two inputs come from resources outside the repository, and their outputs are committed. `figures/src/audit_md_runs.py` reads each run's topology to record what it simulated; it needs the trajectory folders and mdtraj. `figures/src/render_structure.py` needs PyMOL. `figA1_construct_map.svg` is drawn by hand; edit it directly. The figures are set in Source Sans 3 (SIL Open Font License), bundled in `figures/fonts/`; the pages show the PNGs, and the SVGs keep live text for editing once the font is installed.

## Naming notes

- **"F156" is Thr156.** The analysis code labels residue 6 of the construct as F156; the sequence there is R-Y-K-**T**-E-L, the RYKTEL motif that starts finger 1. Five other per-residue columns carry similarly shifted names. The figures and tables use the residue that is actually present ([full mapping](docs/02_MD_results.md#naming-in-the-analysis-code)).
- **C208A and D210E** follow Akçeşme et al. 2025. In the numbering used here (construct residue + 150), the AF3 jobs mutate C212 (a finger 2 zinc ligand) and D219 (the C-terminal tail). The MD code calls the control C212A.
- **MD folder names describe the plan, not always the run.** Every Anvil folder holds a wild-type simulation; `figures/data/md_run_inventory.csv` records what each run contains.



## Acknowledgments

- The project was possible thanks to **Dr. Silvia B. V. Ramos from UNC** and **Dr. Cheng Biao Wu from UCSD**. They provided helpful feedback consistently and fostered discussions throughout the project.
- EMSA labels: Akçeşme B, Hekimoğlu H, Chirasani VR, İş Ş, Atmaca HN, Waldern JM, Ramos SBV. Identification of deleterious non-synonymous single nucleotide polymorphisms in the mRNA decay activator ZFP36L2. *RNA Biology* 22(1):1-15 (2025). [doi:10.1080/15476286.2024.2437590](https://doi.org/10.1080/15476286.2024.2437590)
- Simulations ran on SDSC Expanse and Purdue Anvil through the ACCESS program
- AlphaFold 3 predictions were made with AlphaFold Server and are subject to the AlphaFold Server Output Terms of Use; a copy sits in each `Data/` job folder.

