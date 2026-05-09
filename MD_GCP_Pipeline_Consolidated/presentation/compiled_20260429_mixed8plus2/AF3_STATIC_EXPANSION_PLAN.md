# AF3 Static Expansion Plan (No New MD Required)

## Priority 1: High-value analyses on existing runs

1. Bootstrap stability profiles per variant
Run bootstrap CIs for pair ipTM and pair PAE across existing model pools to quantify uncertainty, not just means.
2. Model-rank sensitivity
Compare conclusions using model_0 only versus all-model pooled metrics to test rank robustness.
3. Run-level variance decomposition
Split variance into between-run and within-run components to identify variants with unstable AF3 behavior.
4. Paired metric discordance map
Flag variants where ipTM improves but PAE worsens (or vice versa), suggesting fragile interfaces.

## Priority 2: New AF3 static runs to request

1. Add replicate runs where any variant has lower run depth than the rest
Balancing run counts improves fair comparisons and statistical power.
2. Sequence-context perturbation panel
For each variant, run at least one altered RNA context (mutant/non-binder control) while keeping Zn entities constant.
3. Zn-entity perturbation sanity checks
Run a limited set with modified ion setup to test whether deleterious/neutral trends are Zn-model dependent.
4. Single-mutation neighborhood scans
Run nearby substitutions around top signals (e.g., R160, C206-210 region) to map local static sensitivity gradients.

## Priority 3: AF3-derived features for stronger classification

1. Contact-map overlap to WT
Compute residue-by-residue protein-RNA contact retention from static structures.
2. Per-residue confidence drop signatures
Track whether confidence loss localizes to zinc fingers, linker, or RNA-contact residues.
3. Interface electrostatics proxy panel
Add charge-change and aromaticity-change descriptors alongside AF3 confidence metrics.

## Suggested execution order (24-hour window)

1. Finish Priority 1 analyses first (no new runs required).
2. Launch only C208A balancing runs next (fastest improvement in dataset quality).
3. Use updated balanced dataset for final presentation plots and conclusions.

