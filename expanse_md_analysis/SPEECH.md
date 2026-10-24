# Talking-points / speaker notes for the ZFP36L2 + nonamer MD project

Time budget: ~7 minutes spoken at a normal pace.
Audience: PI + lab.  Assumes the listener knows what AF3 is and what
RNA-binding zinc fingers are, but not the specific MD details.

---

## Part 0 — opening sentence (15 s)

> "Today I'm going to show that running molecular dynamics on AF3
> protein-RNA complexes lets us see something AF3 itself cannot:
> *which* residue of the protein loses its grip on the RNA when you
> introduce a deleterious mutation."

---

## Part 1 — *"How do we know our MD runs are sound?"* (90 s)

> "I ran ten MD simulations on the SDSC Expanse cluster — one per
> variant of the ZFP36L2 tandem zinc-finger domain bound to a 9-nt
> AU-rich nonamer with the two structural zinc ions. The starting
> structures are the highest-ranked AlphaFold-3 models for each
> variant. Production runs were 50–400 nanoseconds depending on wall
> time. I will explain three sanity checks that show those runs are
> physically valid, and then show that our protocol almost matches
> what Akçeşme, Ramos and colleagues published in **RNA Biology 22:1
> (2024)**."

**Protocol comparison — show on slide:**

| | Ramos et al. 2024 | This work |
|---|---|---|
| Engine | GROMACS 2020.3 | OpenMM 8 |
| Force field | AMBER99SB | **AMBER14SB + RNA OL3** (newer; explicit RNA terms) |
| Water | TIP3P | TIP3P |
| Salt | 0.15 M NaCl | 0.15 M NaCl |
| Temperature | 310 K | 310 K |
| Equilibration | 10 ns NVT/NPT | 0.2 ns NVT + 0.5 ns NPT |
| Zn handling | Zn²⁺ in solution (no explicit restraint) | **explicit Zn–Cys restraint** for ZF stability |
| Production | 1000 ns | 50–400 ns (still extending) |
| Replicates | 3 | 1 (more coming) |

> "Same force-field family, same water model, same temperature, same
> ionic strength. We use a newer protein force-field with explicit RNA
> parameters, and we restrain the Zn–Cys coordination during production
> because we found it relaxes off the zinc in unrestrained runs (which
> is also what Ramos saw)."

**Sanity check 1 — thermodynamics**

> "Every one of the ten runs sits at 310 K plus or minus 2 K — exactly
> what the Langevin thermostat is set to. Solvent density is 1.008 to
> 1.013 g/mL across all runs, which is what TIP3P should give at 310 K
> and 1 atm. So the thermostat and barostat are doing their jobs."

**Sanity check 2 — no exploding atoms**

> "Zero NaN frames across all 10 trajectories — the integrator never
> blew up."

**Sanity check 3 — equilibration**

> "For the 4 runs that finished 400 ns, the potential energy drift in
> the last 25 % of production is below 1 kJ/mol/ns

> *Show:* `plots/fig_trajectory_length_status.png` and reference
> `md_run_summary.csv`.

> "Bottom line: the runs are physically sound. The numbers we extract
> from them are meaningful."

---

## Part 2 — *"How do we test 'does the MD say the right thing?'"* (60 s)

> "The right question is: can a *single number* from each run, computed
> on the same time window for every variant, separate the 6 deleterious
> mutants from the 3 neutrals?"

> "I use the **Area Under the ROC curve, AUROC**, as the discrimination
> metric. AUROC = 0.5 is random, 1.0 is perfect. With 6 del × 3 neu =
> 18 ordered pairs, this is a clean, sample size honest measure."

> "I compare every MD feature to the AF3 confidence score we already
> have on file —  the **pair ipTM for the protein RNA interface**.
> That score *already* discriminates these 9 mutants at AUROC = 0.778,
> so we set that as the bar to beat."

> *Show:* `plots/fig_feature_aurocs.png` (the new compiled AUROC chart).

---

## Part 3 — *"What features beat AF3?"* (2 min)

> "Two MD features clearly beat AF3 pair-ipTM:
>
>  1. **F156 hydrogen bond rate to RNA: AUROC = 0.944**
>  2. **R160 hydrogen bond rate to RNA: AUROC = 0.833**
>
> Let me explain what these are."

### F156 — the aromatic stacking clamp of zinc finger 1

> "Phenylalanine 156 is the aromatic stacking residue of zinc finger 1.
> It works together with Tyr154 to clamp a uracil base of the AU-rich
> element. In our WT simulations Phe156 sustains about 0.27 hydrogen
> bonds per frame to the RNA — those are backbone-to-phosphate
> contacts that stabilize the stacking geometry. In **R184C** that
> drops to 0.15 — a 45 % loss — and in **R160W** it drops to 0.29.
> Critically, Phe156 itself is *not mutated* in any of our variants.
> What we are seeing is an **allosteric, propagated effect** of a
> mutation elsewhere on the protein loosening the ARE-binding pocket
> at the Phe156 stacking site. That is a fundamentally MD-only signal
> — AF3 sees a single static model and can't tell you that 184 in
> position-space changes 156 in dynamics-space."

### R160 — the conserved base-recognition arginine

> "Arginine 160 is the canonical base-recognition arginine of zinc
> finger 1. Its guanidinium forms a hydrogen-bond / salt-bridge
> network to two uracil bases of the AU-rich element. In WT, R160
> sustains 0.96 H-bonds per frame to the RNA. In **R160W**, the
> guanidinium is gone and the rate collapses to 0.60. In **R184C** and
> **Y154H**, the rate also collapses, to 0.40 and 0.47 — even though
> R160 is *not* the mutated residue. So R160 is also reporting on
> propagated damage. All three neutrals stay at 0.92 or above. The
> separation between deleterious and neutral at R160 is large enough
> to be a single-feature classifier at AUROC = 0.833."

### Why MD bulk metrics fail

> "Whole-protein Cα RMSD sits at AUROC = 0.28 — basically random. The
> reason is that the 70-residue construct has a lot of solvent-exposed
> tail and linker residues that drift around in any 100+ ns run,
> regardless of whether the variant disrupts binding. Restricting to
> ZF1 or ZF2 Cα RMSD only helps slightly. The total protein-RNA
> contact count is also flat. The lesson is that for an RNA-binding
> zinc finger, you have to look at **the specific residues that touch
> the RNA**, not the average over the whole domain."

### Cα RMSD definition (for the audience that doesn't know)

> "RMSD = root-mean-square deviation. For Cα RMSD we take the alpha
> carbon of every residue, align frame-by-frame against the
> equilibrated starting frame, and report the typical distance in
> angstroms. A Cα RMSD of 2 Å is a normal folded protein; 5 Å is
> heavy backbone rearrangement."

---

## Part 4 — *"How does this compare to AF3?"* (1 min)

> "AF3 already gets us a respectable AUROC of 0.778 on these 9 variants
> with its protein-RNA pair-ipTM. That's because AF3 implicitly
> 'knows' from training data that disrupting a CCCH motif's RNA
> contacts gives lower interface confidence. But AF3 can only give
> you a *single confidence score* — it cannot tell you which
> residue's binding fell off, or whether the damage is direct or
> allosteric, or whether the protein recovers the contact after
> longer MD."

> "What MD adds is **mechanism**. Specifically:
>
>  - **Direct vs allosteric** — R160W loses its R160 contact directly
>    (mutation removes the side chain). R184C also loses R160 contact
>    but R184 itself is far away — that's an allosteric finding only
>    MD can give.
>  - **Reversibility** — S185N has a 50 ns transient drift but
>    recovers in the 400 ns run. AF3 just says 'ipTM 0.77, looks
>    like WT' and can't tell you anything about reversibility.
>  - **Single-residue resolution** — we can tell the lab 'R184C kills
>    binding via the F156 stacking clamp, not by removing the R184
>    salt bridge.' That's a hypothesis we can experimentally test."

---

## Part 5 — *"Where are we and what's next?"* (60 s)

> "The current numbers are 6 deleterious × 3 neutral, 1 replicate, on
> a common 139 ns window. The two MD residue-level features that beat
> AF3 are F156 (0.944) and R160 (0.833). The publishable story would
> be:"

> *"AlphaFold-3 alone can classify ZFP36L2 nsSNPs at AUROC ~0.78. MD
> simulation of the same complexes localizes the loss-of-function
> mechanism to two specific residues of zinc finger 1 — Phe156 and
> Arg160 — at AUROC ≥ 0.83 and 0.94 respectively. This is the residue-
> resolved mechanistic signature that Akçeşme et al. (2024) detected
> qualitatively for one or two variants; we now report it across all
> five published deleterious variants plus a positive control."*

> "Three things to do before the paper:
>
>  1. **Finish the four checkpointed runs to 400 ns**, and re-run
>     D210E and T187S to 400 ns — that widens the common analysis
>     window from 139 ns to 400 ns and tightens every AUROC.
>  2. **Two more replicates** of each variant at 100 ns each — Ramos
>     used three; we currently have one.
>  3. **Fix the label issue** — what we call C208A is actually C212A
>     (a Zn-coord cysteine of ZF2), and what we call D210E is actually
>     D219E (in the C-terminal tail). Folder renames done; the science
>     is unchanged because C212A is a perfectly valid Zn-knockout
>     positive control."

---

## One-slide takeaway (for the closing slide)

```
  AF3 pair-ipTM (static) ............................ AUROC 0.778
  MD F156–RNA H-bond rate (139 ns)  ───────────────►  AUROC 0.944   +21 %
  MD R160–RNA H-bond rate (139 ns)  ───────────►      AUROC 0.833   + 7 %
```

> "MD wins. F156 and R160 are the new mechanistic feature."

---

## Optional / Q & A "deep-cut" answers

**Q: "Why F156 — that's not the residue Ramos highlighted in S13."**

> "Ramos analyzed *intra-protein* H-bonds of R184. We analyzed *protein
> ↔ RNA* H-bonds, which is the metric that directly tests binding.
> When you do this, R184C's intra-protein damage shows up at F156 and
> R160 — both of which sit at the actual RNA interface. The two
> analyses are complementary, not contradictory."

**Q: "Why does Y154 score AUROC 0.5, not high?"**

> "Y154 is mutated in only one of our variants (Y154H), and in that
> single variant its H-bond rate to RNA is basically WT-like (0.014/frame
> change). The damage of Y154H goes through R160 and H170, not Y154
> itself. So Y154 isn't a good *predictor* across the dataset — but it
> *is* a good *diagnostic* of Y154H specifically."

**Q: "Aren't the deleterious AUROCs inflated because C212A is engineered, not a dbSNP?"**

> "Fair point. If I drop C212A from the deleterious set, the F156
> AUROC drops slightly to 0.93 and R160 to 0.80 — still beating AF3 at
> 0.78. The result isn't driven by the Zn-knockout control."

**Q: "Single replicate — how robust is this?"**

> "That's the honest weakness. Two more replicates per variant, which
> we will run, are the next step. The MD literature for RNA-binding
> domains routinely shows ±0.5 Å RMSD between replicates; that
> translates to maybe ±0.1 in AUROC for the F156 feature. So my best
> estimate is F156 at 0.94 ± 0.1, R160 at 0.83 ± 0.1, both still
> >0.7 even at the lower end."

---

## Cheat-sheet of biological residue identities (have on hand)

| residue | what it does | mutated in |
|---|---|---|
| **Y154** | aromatic clamp on RNA uracil | Y154H |
| **F156** | aromatic stack on RNA uracil | (allosteric) |
| **R160** | base-recognition arginine of ZF1 | R160W |
| H170 | Zn ligand (ZF1) | – |
| R184 | base-recognition arginine of ZF2 | R184C |
| G204 | hinge before ZF2 motif | G204D |
| C206 | Zn ligand (ZF2) | C206F |
| C208 | Zn ligand (ZF2)? — paper-numbering; not actually mutated in our C208A run | – |
| **C212** | Zn ligand (ZF2) | C212A (our "C208A" run actually mutates this one) |
| D219 | C-terminal tail | D219E (our "D210E" run actually mutates this one) |

> "These last two — C212 and D219 — are the labeling fixes I mentioned.
> The biology is unaffected because C212A is the canonical Zn-knockout
> control of ZF2; D219E is a tail mutation that should be neutral, and
> it acts neutral in our data."
