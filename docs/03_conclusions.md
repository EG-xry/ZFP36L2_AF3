# 3 · Conclusions

[← MD results](02_MD_results.md) · [Back to overview](../README.md)

Structure prediction followed by simulation is often presented as a way to triage variants before anyone runs a gel. This project tested that idea on about as clean a case as exists: one small domain, one RNA, and eight variants whose effects had been measured by EMSA. My answer is a qualified no. MD did find a readout that ranks the variants better than AlphaFold 3's confidence; however, the ranking is not significant, it does not trace to a mechanism, and it moves with ordinary analysis choices. The limiting factor was never simulation length. It was the number of labeled variants.

## What held up

- **AF3 builds a plausible bound complex for every variant.** Its interface confidence also drops sharply for a non-binding RNA (0.765 → 0.376), so the predictor is sensitive to binding, just not to these point mutations.
- **The best MD readout is stable under the checks available.** T156–RNA hydrogen bonds rank the variants with an AUROC of 0.800. That value holds after a standardized re-analysis and after dropping each variant's weakest replicate.
- **The negative results are informative in their own right.** Bulk RMSD, contact counts, and the native-contact Q-score do not separate the classes. Neither does aromatic stacking. MM-GBSA without a zinc model is not physical.

## What did not

- **Significance.** The best readout reaches p = 0.125. Nothing in the 32-readout screen reaches even an uncorrected 0.05.
- **Mechanism.** T156 reads the same base in both classes, stacking does not discriminate, and RYKTEL-region contacts move the wrong way. The signal is a diffuse 0.02 shift with no specific contact behind it.
- **Stability across analysis choices.** Defensible alternatives move the AUROC anywhere from 0.53 to 0.93.

## Eight variants cannot carry a significant result

<p align="center">
  <img src="../figures/figC1_n8_ceiling.png" alt="Left: exact p-value for every attainable AUROC with 5 vs 3 variants; even a perfect ranking gives p = 0.018, above the Bonferroni line. Right: minimum AUROC needed versus number of validated variants; the observed 0.800 needs about 34 variants to clear Bonferroni." width="100%">
</p>

Five deleterious and three neutral variants can be ranked in only 56 distinct ways. A perfect ranking is one of them, so its p-value is 1/56 = 0.018. The project screened 32 readouts, which puts the Bonferroni bar at 0.0016. That bar lies below anything eight variants can produce. More replicates, longer trajectories, or a better readout cannot change this arithmetic.

Panel b turns the question around and asks how many validated variants an AUROC of 0.800 would need:

- **14 variants (7 + 7)** for an uncorrected p < 0.05 on a single readout chosen in advance.
- **About 34 variants (17 + 17)** to survive correction for a 32-readout screen like this one.

## Reasonable analysis choices move the answer

<p align="center">
  <img src="../figures/figC2_forking_paths.png" alt="AUROC for eight defensible analysis variants: T156 at 0.933, 0.533, 0.867 and 0.800 depending on window and replicates; AF3 pair-ipTM at 0.733 or 0.633 and min pLDDT at 0.667 or 0.933 depending on model pooling." width="100%">
</p>

Every row of this figure uses the same eight variants and a choice I could defend in a methods section.

- **Window length.** The same ten single trajectories give 0.933 over their first 139 ns and 0.533 over their first 269 ns. Had I stopped at 139 ns, I would have reported a nominally significant result.
- **Sampling and inclusion rules.** Pooling 2–4 replicates gives 0.867 at 0.1 ns sampling and 0.800 at 1 ns sampling once partial runs are excluded.
- **Model pooling in AF3.** Averaging all 40 models gives pair-ipTM 0.733; using each seed's top model gives 0.633. The same switch lifts min pLDDT from 0.667 to 0.933, a number that measures the disordered C-terminal residue.

With only 15 deleterious–neutral pairs, a single pair changing order moves the AUROC by 0.067. Analysis choices that look minor can therefore produce, or erase, an apparent result. This is the garden of forking paths in miniature (Gelman and Loken 2014), and the honest response is to fix one analysis before looking.

## Why I paused

There is a version of this project that keeps going: more replicates, longer runs, an in-silico ablation of T156 to test causation. I weighed each. More replicates would tighten the T156 estimate but cannot raise the ceiling that eight labels impose. An ablation needs a mechanism to test, and the decomposition found none. To be fair, a null result at n = 8 does not show that MD cannot rank these variants; it shows that this dataset cannot tell either way. That is a reason to stop spending GPU hours on it, not a verdict on the method.

## Lessons

1. **Count labels before nanoseconds.** The number of experimentally validated variants sets the best p-value any method can reach. Check that ceiling before planning the simulations.
2. **Fix the analysis before looking.** Choose one primary readout, window, sampling interval, and pooling rule in advance, then score it once.
3. **Check residue identities against the sequence on day one.** Six of the fourteen per-residue labels in the analysis code were wrong, and an "F156 aromatic stacking" story turned out to be about a threonine.
4. **Read AF3 summary metrics before trusting them.** Here "min pLDDT" was the confidence of the last residue in the construct, and the way models were pooled changed which AF3 metric looked best.
5. **Treat single-replicate MD windows as exploratory.** The same trajectories went from 0.933 to 0.533 as the window grew.
6. **Model the metal if energies matter.** Stripping zinc to fit a standard MM-GBSA workflow gave interaction energies around +12,000 kJ/mol.

## What would justify resuming

- **More EMSA-validated variants.** About 14 balanced variants to test one readout fixed in advance, or about 34 to support a screen. Labels predicted by other tools would not help; scoring against another predictor's calls is circular.
- **One readout fixed in advance.** For example, T156–RNA hydrogen bonds at 1 ns sampling over runs of at least 180 ns, scored once on the new variants.
- **A bonded zinc model** (ZAFF or MCPB.py), so that binding free energies become meaningful.
- **A longer ARE** to pin the RNA register, which should shrink the spread between replicates that currently swamps the signal.

Eight labeled variants could not answer the question I asked, but they did show what an answer would cost. More nanoseconds sharpen each estimate; only more labels can make the comparison mean something.

## References

- Akçeşme B, Hekimoğlu H, Chirasani VR, İş Ş, Atmaca HN, Waldern JM, Ramos SBV. Identification of deleterious non-synonymous single nucleotide polymorphisms in the mRNA decay activator ZFP36L2. *RNA Biology* 22(1):1–15 (2025). doi:10.1080/15476286.2024.2437590
- Gelman A, Loken E. The statistical crisis in science. *American Scientist* 102(6), 460–465 (2014).
- Mann HB, Whitney DR. On a test of whether one of two random variables is stochastically larger than the other. *Ann Math Stat* 18, 50–60 (1947).
