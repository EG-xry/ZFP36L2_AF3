# H-bond Pilot Summary (5 Paper Variants)

This is a pilot comparison on short available trajectories, not the paper's 1000 ns GROMACS setup.

## Directional comparison vs WT counterpart residues

- Y154H: mutant mean `1.20` vs WT `0.56` (Δ `+0.64`, higher)
- R160W: mutant mean `0.92` vs WT `1.28` (Δ `-0.36`, lower)
- R184C: mutant mean `0.76` vs WT `1.48` (Δ `-0.72`, lower)
- G204D: mutant mean `0.00` vs WT `0.32` (Δ `-0.32`, lower)
- C206F: mutant mean `0.04` vs WT `0.28` (Δ `-0.24`, lower)

## Caveats

- Uses short trajectories from current OpenMM cloud outputs.
- Uses geometric H-bond detection in MDAnalysis, not `gmx hbond`.
- No replicate-level uncertainty estimates yet.

