# ZFP36L2 AF3 — Averaged metrics (model_0 only; counts by variant: WT=8, C206F=8, C208A=8, G204D=8, R160W=8, R160E=5, R184C=8, Y154H=8, H212A=5, D210E=8, S185N=8, T187S=8)

- Variants analysed: 12
- Models per variant: WT=8, C206F=8, C208A=8, G204D=8, R160W=8, R160E=5, R184C=8, Y154H=8, H212A=5, D210E=8, S185N=8, T187S=8
- model_0 only; counts by variant: WT=8, C206F=8, C208A=8, G204D=8, R160W=8, R160E=5, R184C=8, Y154H=8, H212A=5, D210E=8, S185N=8, T187S=8
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
| WT        | Nonamer 9nt (WT) + 2Zn | WT                | WT            |                 0.7775  |                   1.05063 | 96.9882 |     98.4575 |     97.8272 |     72.47   |  1.92297e-07 |                 0.0046291  |                    0.0351971 |  0.189203  |      0.069126  |      0.124834  |        4.15843 |     5.57881e-08 |          8 |
| C206F     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.72125 |                   1.22188 | 95.5677 |     98.084  |     94.8099 |     66.1012 |  1.27787     |                 0.00834523 |                    0.068917  |  0.283976  |      0.0975728 |      0.501956  |        2.04893 |     0.257448    |          8 |
| C208A     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.75875 |                   1.16125 | 97.1424 |     98.4249 |     97.9561 |     73.2013 |  0.801052    |                 0.00834523 |                    0.0467325 |  0.130862  |      0.0619068 |      0.0528696 |        4.06894 |     0.17277     |          8 |
| G204D     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.77875 |                   1.09563 | 96.9344 |     98.2832 |     97.9484 |     70.1012 |  0.700097    |                 0.0064087  |                    0.062016  |  0.0878357 |      0.0726874 |      0.0247724 |        1.38786 |     0.149932    |          8 |
| H212A     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.776   |                   1.039   | 96.915  |     98.2572 |     97.6183 |     78.812  |  1.10184     |                 0.00547723 |                    0.0357771 |  0.215456  |      0.0882872 |      0.141659  |        3.88982 |     0.176976    |          5 |
| R160E     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.774   |                   1.08    | 97.0357 |     98.4068 |     97.9468 |     72.102  |  0.606333    |                 0.00547723 |                    0.027157  |  0.0905451 |      0.020548  |      0.0584055 |        1.44325 |     0.138947    |          5 |
| R160W     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.76375 |                   1.26625 | 97.0879 |     98.4261 |     97.9643 |     74.2262 |  0.766822    |                 0.0176777  |                    0.13695   |  0.168063  |      0.0696657 |      0.104469  |        4.06317 |     0.27055     |          8 |
| R184C     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.77    |                   1.06812 | 96.7951 |     98.2492 |     97.8457 |     70.5763 |  0.656906    |                 0.00755929 |                    0.0352478 |  0.206626  |      0.0973457 |      0.177137  |        2.85077 |     0.0937631   |          8 |
| Y154H     | Nonamer 9nt (WT) + 2Zn | deleterious       | WT            |                 0.785   |                   1.06563 | 97.0031 |     98.3152 |     98.0325 |     70.6037 |  0.740372    |                 0.00755929 |                    0.0271816 |  0.110944  |      0.174016  |      0.0668506 |        2.61122 |     0.159069    |          8 |
| D210E     | Nonamer 9nt (WT) + 2Zn | neutral           | WT            |                 0.77625 |                   1.08062 | 96.9536 |     98.272  |     97.8251 |     77.1325 |  0.86591     |                 0.00916125 |                    0.0437066 |  0.404741  |      0.115685  |      0.325139  |        4.59857 |     0.321206    |          8 |
| S185N     | Nonamer 9nt (WT) + 2Zn | neutral           | WT            |                 0.78125 |                   1.07875 | 97.1186 |     98.4241 |     98.0479 |     73.605  |  0.711024    |                 0.00991031 |                    0.0510427 |  0.155146  |      0.0754355 |      0.0836495 |        4.54792 |     0.112434    |          8 |
| T187S     | Nonamer 9nt (WT) + 2Zn | neutral           | WT            |                 0.77    |                   1.09313 | 96.9991 |     98.366  |     97.8429 |     74.5187 |  0.703551    |                 0.0106904  |                    0.041052  |  0.182092  |      0.123761  |      0.130786  |        2.1528  |     0.134716    |          8 |
