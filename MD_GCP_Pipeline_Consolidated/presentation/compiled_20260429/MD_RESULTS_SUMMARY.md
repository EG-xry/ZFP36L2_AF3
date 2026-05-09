# MD Results Summary (Mixed 8+2 Dataset)

## Data used

- 8 variants from main cloud batch output.
- 2 variants from fallback short rerun: S185N and C208A.
- WT baseline is the WT run in the same synced cloud output set.

## WT baseline

- WT protein drift (CA RMSD): `1.9605 A`
- WT RNA drift (P RMSD): `2.3229 A`

## Mean deltas vs WT (negative = lower drift than WT)

- Deleterious, protein delta: `-0.4613 A`
- Neutral, protein delta: `-0.3540 A`
- Deleterious, RNA delta: `-0.2905 A`
- Neutral, RNA delta: `-0.4692 A`

## Neutral vs Deleterious class comparison

- Welch t-test (protein delta): `p=0.5344`
- Welch t-test (RNA delta): `p=0.3573`

Interpretation: with this small mixed dataset, there is no statistically strong class separation between neutral and deleterious in MD drift deltas.

