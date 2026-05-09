# ZFP36L2 AF3 — Averaged metrics (all models; counts by variant: WT=40, C206F=40, C208A=40, G204D=40, R160W=40, R160E=25, R184C=40, Y154H=40, H212A=25, D210E=40, S185N=40, T187S=40)

- Variants analysed: 12
- Models per variant: WT=40, C206F=40, C208A=40, G204D=40, R160W=40, R160E=25, R184C=40, Y154H=40, H212A=25, D210E=40, S185N=40, T187S=40
- all models; counts by variant: WT=40, C206F=40, C208A=40, G204D=40, R160W=40, R160E=25, R184C=40, Y154H=40, H212A=25, D210E=40, S185N=40, T187S=40
- `_sd` columns are standard deviations (ddof = 1)

## Metric definitions

| Metric | Description |
|--------|-------------|
| pair_ipTM_protein_rna | Pairwise protein↔RNA ipTM (interface confidence; higher = better) |
| pair_pae_min_cross_mean | Mean of PAE(protein→RNA) and PAE(RNA→protein); lower = better |
| pLDDT | Mean per-residue Cα confidence (0–100) |
| pLDDT_ZF1 / ZF2 | Mean Cα pLDDT for zinc finger 1 / 2 |
| pLDDT_min | Minimum per-residue Cα pLDDT across the protein domain (matches AF3 server coloring) |
| RMSD_whole | Backbone RMSD vs WT reference (Å) |

| Protein   | RNA                    | ProteinExpClass   | RNAExpClass   |   pair_ipTM_protein_rna |   pair_pae_min_cross_mean |   pLDDT |   pLDDT_ZF1 |   pLDDT_ZF2 |   pLDDT_min |   RMSD_whole |   pair_ipTM_protein_rna_sd |   pair_pae_min_cross_mean_sd |   pLDDT_sd |   pLDDT_ZF1_sd |   pLDDT_ZF2_sd |   pLDDT_min_sd |   RMSD_whole_sd |   n_models |
|:----------|:-----------------------|:------------------|:--------------|------------------------:|--------------------------:|--------:|------------:|------------:|------------:|-------------:|---------------------------:|-----------------------------:|-----------:|---------------:|---------------:|---------------:|----------------:|-----------:|
| WT        | Nonamer 9nt (WT) + 2Zn | WT                | WT            |                 0.76525 |                   1.0755  | 96.9287 |     98.3815 |     97.6739 |     76.2472 |   1.8665e-07 |                 0.0113199  |                    0.0462269 |   0.34846  |      0.211938  |       0.2372   |        7.01296 |     5.57693e-08 |         40 |
| C206F     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.70975 |                   1.27225 | 95.3965 |     97.9216 |     94.6162 |     66.5977 |   1.13568    |                 0.00999679 |                    0.0807493 |   0.272419 |      0.238255  |       0.477254 |        2.83029 |     0.299613    |         40 |
| C208A     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.7515  |                   1.19362 | 97.0297 |     98.3515 |     97.8359 |     74.8405 |   0.752401   |                 0.0114466  |                    0.0770613 |   0.24718  |      0.11454   |       0.191374 |        5.01329 |     0.16797     |         40 |
| G204D     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.76675 |                   1.11962 | 96.8363 |     98.2307 |     97.8112 |     71.3162 |   0.674892   |                 0.0107148  |                    0.0575057 |   0.200742 |      0.143627  |       0.153468 |        2.77663 |     0.199795    |         40 |
| H212A     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.768   |                   1.0536  | 96.8053 |     98.2034 |     97.4646 |     78.9624 |   1.22331    |                 0.0135401  |                    0.0306703 |   0.364391 |      0.193942  |       0.415835 |        3.14671 |     0.254006    |         25 |
| R160E     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.7616  |                   1.1064  | 96.8469 |     98.3637 |     97.7194 |     72.6672 |   0.798752   |                 0.017      |                    0.0390384 |   0.235101 |      0.0539791 |       0.271878 |        3.15993 |     0.241455    |         25 |
| R160W     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.75575 |                   1.32012 | 97.0451 |     98.3867 |     97.8933 |     75.7165 |   0.757798   |                 0.0125856  |                    0.11297   |   0.194624 |      0.132467  |       0.13453  |        5.14606 |     0.249819    |         40 |
| R184C     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.75825 |                   1.103   | 96.6938 |     98.2151 |     97.7135 |     70.8413 |   0.712438   |                 0.0127877  |                    0.0380755 |   0.17571  |      0.119996  |       0.174822 |        2.99159 |     0.20167     |         40 |
| Y154H     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.77525 |                   1.08975 | 97.044  |     98.3523 |     97.9768 |     73.559  |   0.75049    |                 0.00876693 |                    0.0383297 |   0.218987 |      0.129097  |       0.150589 |        5.6034  |     0.207656    |         40 |
| D210E     | Nonamer 9nt (WT) + 2Zn | neutral           | WT            |                 0.76925 |                   1.0825  | 96.7838 |     98.1449 |     97.7115 |     74.547  |   0.89334    |                 0.0111832  |                    0.0314398 |   0.315351 |      0.325416  |       0.265798 |        3.83086 |     0.237761    |         40 |
| S185N     | Nonamer 9nt (WT) + 2Zn | neutral           | WT            |                 0.77075 |                   1.12163 | 97.0246 |     98.3649 |     97.9495 |     72.9352 |   0.710874   |                 0.0114102  |                    0.0716508 |   0.208801 |      0.167247  |       0.145739 |        3.71344 |     0.205594    |         40 |
| T187S     | Nonamer 9nt (WT) + 2Zn | neutral           | WT            |                 0.76225 |                   1.11088 | 96.7984 |     98.3062 |     97.6835 |     73.361  |   0.722231   |                 0.0164063  |                    0.054076  |   0.272699 |      0.148622  |       0.261032 |        3.54875 |     0.199376    |         40 |
