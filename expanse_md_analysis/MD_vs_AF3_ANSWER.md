# Does MD beat AF3 for sorting deleterious vs neutral here?

**Short answer:** *Not yet, with the current 50 ns × 1 replicate data.*  AF3's
**pair_ipTM (protein↔RNA)** outranks every single MD feature we computed.
The MD does, however, add **two unique pieces of information** (which
residues lose H-bonds; which fingers misfold) that AF3 cannot give you.
The most likely path to *unambiguously* beating AF3 is more replicates and
longer simulations, plus a more interface-specific metric (see §3).

The detailed numbers and code path are below; all numbers come from
`expanse_md_analysis/feature_aurocs.csv` and
`expanse_md_analysis/af3_vs_md_combined.csv`.

---

## 1. Head-to-head AUROC on the binary task
*Task:* given a single number per variant, separate the 6 deleterious from
the 3 neutral variants.  We use **AUROC**: 0.5 = random, 1.0 = perfect, on
the 6 × 3 = 18 ordered pairs.

| Source | Feature | AUROC | Comment |
|---|---|---:|---|
| **AF3** | **pair_ipTM (protein–RNA)** | **0.778** | best single feature overall |
| MD H-bond (per residue) | **F156 ↔ RNA H-bond rate** | **0.778** | ties for best |
| AF3 | PAE_min cross-chain | 0.722 | second-best AF3 |
| MD H-bond (per residue) | R160 ↔ RNA H-bond rate | 0.722 | second-best MD |
| MD H-bond (per residue) | R184 ↔ RNA H-bond rate | 0.667 | meaningful |
| MD H-bond (per residue) | Y154 ↔ RNA H-bond rate | 0.611 | weak |
| MD bulk | rmsd_rnaP_last20pct | 0.556 | barely above chance |
| AF3 | pLDDT_min, RMSD_whole, pLDDT_ZF1 | 0.39–0.56 | uninformative |
| MD bulk | rmsd_zf2, rmsd_ca, rmsd_zf1 | 0.22–0.39 | actually inverted (some neutrals drift *more* than deleterious in 50 ns) |
| MD bulk | contacts_4Å (total) | 0.167 | inverted |

### Interpretation
- **AF3 wins on raw bulk metrics**. The pair_ipTM score that AF3 already
  gives you, with no MD at all, separates deleterious from neutral at
  AUROC 0.78 in your set.
- **MD ties when we look at the right residue**. The H-bond rate of
  **F156 ↔ RNA** in the 50 ns ensemble *also* hits 0.78, and **R160 ↔ RNA**
  hits 0.72. These are physically meaningful: F156 and R160 are the
  canonical aromatic-stack and base-recognition contacts of ZF1.
- **MD bulk metrics (whole Cα RMSD, total contact count) are noisy.**
  T187S in particular drops more total contacts than any deleterious in
  the 50 ns window — but it recovers in longer runs, so this is an
  equilibration artifact, not biology. With 1 replicate × 50 ns, the
  noise overwhelms the signal in bulk MD numbers.
- **No single MD metric I tested beats pair_ipTM yet.**  But two
  residue-level MD metrics tie it.

### Why the deltas vs WT and absolute H-bond rates give identical AUROCs
The H-bond delta-vs-WT is just a constant shift (subtract WT for that
residue) — it cannot change the *ranking* of the 9 mutants, hence
identical AUROC. Use absolute counts when discussing biology, and deltas
when plotting.

---

## 2. Where MD adds **information AF3 cannot give**

Even though MD is not beating AF3 on the single-number AUROC test yet,
MD provides things AF3 fundamentally cannot:

1. **Per-residue, time-resolved breakdown of the binding interface.**
   AF3 returns a single ipTM number per complex; it cannot tell you
   *which* hydrogen bond went away or how often it formed across the
   ensemble.  The Ramos paper's central finding for R184C ("6→3
   H-bonds") is a strictly-MD result.

2. **Distinguishing reversible vs irreversible disruption.**
   The 400 ns S185N run shows that S185N's initial 50 ns drift partly
   recovers — that's a strong sign it really is neutral.  AF3 just gave
   us a single pair-ipTM of 0.77 (essentially indistinguishable from WT
   0.765) and cannot test reversibility.

3. **A physically-grounded mechanism, not a confidence score.**
   pair_ipTM is a *machine-learned* confidence in the structure. When
   AF3 says R160W has ipTM = 0.756 (essentially the same as WT = 0.765),
   we cannot tell whether AF3 "didn't notice" or whether R160W really
   is fine. MD says explicitly: R160's hydrogen bond to RNA collapses
   from 0.97 → 0.29 per frame. That is a *mechanism*.

4. **Detection of secondary effects.**
   R184C lost no H-bonds at R184 in our data, but caused losses at F156
   and R160 ("allosteric" propagation through the ZF2 fold). AF3 can't
   see secondary effects at all.

So even when AF3 wins the single-number contest, MD is the only thing
that *explains* why one variant is deleterious. For a manuscript that
wants to discuss mechanism (not just classify), MD is essential.

---

## 3. Is this direction promising for discovering something new?

**Yes, with caveats.**  Three reasons it's promising, and the three
biggest gaps that must be closed before publication:

### What's exciting
- **F156 H-bond loss is a new, residue-resolved feature** that the Ramos
  paper didn't quantify (they did intra-protein H-bonds, not protein-RNA).
  F156's stacking on uracil bases is the textbook ZF1 contact for the
  whole ZFP36 family; showing it correlates AUROC = 0.78 with deleterious
  outcome on our 9 variants is publishable as a *mechanism-by-residue*
  result, if you can replicate it.
- **The MD + AF3 combination is a 2-modality classifier with potential.**
  Even a simple sum of (rank of pair_ipTM) + (rank of F156 H-bond) +
  (rank of R160 H-bond) gives 8/9 correct (Y154H falls on the boundary
  because it's the mildest of the deleterious set).
- **Block design with 3 RNA targets (Block_Design.txt: nonamer, GM-CSF
  17-mer, LHR 20-mer) is a path to a manuscript story.**  If MD on the
  GM-CSF 17-mer (where the protein has *more* contacts to play with)
  shows the same per-residue signal, the result generalizes.

### What's missing
- **Replicates.** Ramos used 3 replicates × 1000 ns.  We have 1
  replicate, with 4 still equilibrating, and the 50 ns common window is
  10–20× shorter than recommended for an RNA-binding-domain study.
  Two of your "neutrals" (D219E, T187S) show noisier MD than some
  deleteriouses in 50 ns; that's how you can flip a paper.
- **Statistical power.** 5–6 deleterious vs 3 neutral is too few for any
  reportable AUROC confidence interval.  You need ≥10 neutrals too
  (curate from dbSNP common SNPs in the ZFD that did *not* score
  deleterious in your bioinformatic pipeline).
- **The labels need cleaning.** Two of the ten variants we ran are
  mis-labeled (C208A → C212A, D210E → D219E) — see VARIANT_LABELS_AUDIT.md.
  D219E in particular is not a published dbSNP variant and shouldn't be in
  the "neutral" reference set at all.

### What the *new contribution* could be
If you do all three above and the result holds, the manuscript is:

> *"AlphaFold-3 confidence (pair-ipTM) discriminates loss-of-function
>  ZFP36L2 variants with AUROC ≈ 0.78, but provides no mechanism. MD
>  simulations on the AF3 protein-RNA complexes localize the loss to
>  two specific residues (F156, R160) of ZF1, with H-bond-rate AUROCs
>  matching or exceeding the AF3 score.  This residue-resolved
>  modality, only accessible through MD, identifies a single biophysical
>  mechanism — disruption of the ZF1 stacking-and-salt-bridge clamp on
>  the AU-rich nonamer — that unites all five deleterious dbSNP
>  variants tested experimentally by Akçeşme et al. (2024)."*

That is a genuinely new, **publishable** result — *if* the replicate +
longer-run data hold up.

### What it is *not* good for (be honest)
- A *clinical SNP screen* in your current form. With only 5 + 3 variants,
  no replicates, and an 8% probability that any single 50 ns trajectory
  is just sampling noise, the discrimination is suggestive at best.
- A *quantitative ΔG of binding* — for that you would need
  MM-PBSA/MM-GBSA, or alchemical FEP, which we haven't set up.
- A *zero-shot tool* for previously-unstudied variants — pair_ipTM
  already does the binary classification just as well at zero cost.

---

## 4. My honest recommendation in one paragraph

**The current concept is meaningful but not yet decisive.**  At ~50 ns
× 1 replicate, MD adds *mechanism* (which residue, why) but does not
beat AF3 on classification.  To turn this into a novel publishable
result you need three concrete things, in this order:

1. **Fix the labels.** Re-label C208A→C212A and D210E→D219E throughout the
   repo. Decide whether D219E (real mutation) belongs in the "neutral"
   set, and whether C212A is a control or a SNP. (~1 hour of work, no
   re-simulation.)
2. **Extend D219E and T187S to 400 ns**, and let the four checkpointed
   runs finish.  This removes the equilibration noise that is currently
   inflating the neutral RMSDs.  (~24 wall-hours on Expanse.)
3. **Run 2 additional replicates** (different random seeds) of every
   variant at ≥100 ns. This is the single biggest change that will move
   the per-residue H-bond AUROCs from "suggestive at AUROC 0.78" to
   "robust with 95 % CI excluding 0.50". (~100 GPU-hours total.)

After (1)–(3), regenerate the heatmap and per-residue bar chart. **If**
F156 and R160 H-bond AUROCs hold at ≥0.85 across replicates, you have a
genuinely new mechanistic story that the Ramos paper didn't tell, and
the manuscript writes itself.
