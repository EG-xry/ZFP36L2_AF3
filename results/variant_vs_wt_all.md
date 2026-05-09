# ZFP36L2 AF3 — Welch's t-test: variant vs WT (all models; counts by variant: WT=40, C206F=40, C208A=40, G204D=40, R160W=40, R160E=25, R184C=40, Y154H=40, H212A=25, D210E=40, S185N=40, T187S=40)

## Methodology

For each metric, an independent two-sample Welch's t-test (unequal
variance) compares each variant's available values against WT's
available values for the same analysis set. Significance: α = 0.05.
For hypothesis-testing only, each group is normalized to assumed n=40
(bootstrap upsampling when fewer values exist; deterministic downsampling when more exist).

## Metric interpretation

| Metric | Higher is … | Meaning |
|--------|-------------|---------|
| pair ipTM (protein–RNA) | **Better** | Pairwise protein↔RNA predicted TM-score; higher = more confident interface |
| PAE cross-mean | **Worse** | Predicted aligned error between protein ↔ RNA (Å); higher = less confident |
| pLDDT | **Better** | Mean per-residue Cα confidence (0–100); higher = more confident structure |
| pLDDT min (Cα) | **Better** | Lowest per-residue Cα pLDDT in the protein domain; higher = fewer disordered regions |

RMSD is excluded from this test because WT RMSD is ~0 by definition
(self-aligned reference), making the comparison trivially significant
for every variant.

## Results

- **"Better than WT = Yes *"** means the variant is *significantly* better than WT for that metric (p < 0.05).
- **"-"** means the difference is not statistically significant.

| Protein   | ProteinExpClass   | Metric                  |   WT_mean |   Variant_mean |   Difference |   t_statistic |     p_value | Significant (p<0.05)   | Direction vs WT   | Better than WT   |   n_assumed_each_group |
|:----------|:------------------|:------------------------|----------:|---------------:|-------------:|--------------:|------------:|:-----------------------|:------------------|:-----------------|-----------------------:|
| C206F     | deleterious       | pair ipTM (protein–RNA) |    0.7652 |         0.7097 |      -0.0555 |      -23.2425 | 1.69771e-36 | Yes                    | lower             | No               |                     40 |
| C206F     | deleterious       | PAE cross-mean (Å)      |    1.0755 |         1.2723 |       0.1967 |       13.3737 | 6.09012e-20 | Yes                    | higher            | No               |                     40 |
| C206F     | deleterious       | pLDDT                   |   96.9287 |        95.3965 |      -1.5322 |      -21.9094 | 5.25919e-34 | Yes                    | lower             | No               |                     40 |
| C206F     | deleterious       | pLDDT min (Cα)          |   76.2472 |        66.5977 |      -9.6495 |       -8.0699 | 1.07369e-10 | Yes                    | lower             | No               |                     40 |
| C208A     | deleterious       | pair ipTM (protein–RNA) |    0.7652 |         0.7515 |      -0.0138 |       -5.4019 | 6.93225e-07 | Yes                    | lower             | No               |                     40 |
| C208A     | deleterious       | PAE cross-mean (Å)      |    1.0755 |         1.1936 |       0.1181 |        8.3136 | 9.18613e-12 | Yes                    | higher            | No               |                     40 |
| C208A     | deleterious       | pLDDT                   |   96.9287 |        97.0297 |       0.1009 |        1.4944 | 0.139535    | No                     | higher            | -                |                     40 |
| C208A     | deleterious       | pLDDT min (Cα)          |   76.2472 |        74.8405 |      -1.4068 |       -1.0321 | 0.305561    | No                     | lower             | -                |                     40 |
| G204D     | deleterious       | pair ipTM (protein–RNA) |    0.7652 |         0.7668 |       0.0015 |        0.6086 | 0.544534    | No                     | higher            | -                |                     40 |
| G204D     | deleterious       | PAE cross-mean (Å)      |    1.0755 |         1.1196 |       0.0441 |        3.7824 | 0.000311305 | Yes                    | higher            | No               |                     40 |
| G204D     | deleterious       | pLDDT                   |   96.9287 |        96.8363 |      -0.0924 |       -1.4528 | 0.151301    | No                     | lower             | -                |                     40 |
| G204D     | deleterious       | pLDDT min (Cα)          |   76.2472 |        71.3163 |      -4.931  |       -4.1347 | 0.000132886 | Yes                    | lower             | No               |                     40 |
| R160W     | deleterious       | pair ipTM (protein–RNA) |    0.7652 |         0.7558 |      -0.0095 |       -3.5495 | 0.000661196 | Yes                    | lower             | No               |                     40 |
| R160W     | deleterious       | PAE cross-mean (Å)      |    1.0755 |         1.3201 |       0.2446 |       12.6751 | 1.7358e-17  | Yes                    | higher            | No               |                     40 |
| R160W     | deleterious       | pLDDT                   |   96.9287 |        97.0451 |       0.1164 |        1.8442 | 0.0699998   | No                     | higher            | -                |                     40 |
| R160W     | deleterious       | pLDDT min (Cα)          |   76.2472 |        75.7165 |      -0.5307 |       -0.3859 | 0.700714    | No                     | lower             | -                |                     40 |
| R160E     | deleterious       | pair ipTM (protein–RNA) |    0.7652 |         0.7632 |      -0.002  |       -0.7072 | 0.481603    | No                     | lower             | -                |                     40 |
| R160E     | deleterious       | PAE cross-mean (Å)      |    1.0755 |         1.0994 |       0.0239 |        2.6795 | 0.00919245  | Yes                    | higher            | No               |                     40 |
| R160E     | deleterious       | pLDDT                   |   96.9287 |        96.8186 |      -0.1101 |       -1.6055 | 0.112764    | No                     | lower             | -                |                     40 |
| R160E     | deleterious       | pLDDT min (Cα)          |   76.2472 |        72.4473 |      -3.8    |       -3.2084 | 0.00233622  | Yes                    | lower             | No               |                     40 |
| R184C     | deleterious       | pair ipTM (protein–RNA) |    0.7652 |         0.7582 |      -0.007  |       -2.5923 | 0.0114065   | Yes                    | lower             | No               |                     40 |
| R184C     | deleterious       | PAE cross-mean (Å)      |    1.0755 |         1.103  |       0.0275 |        2.9041 | 0.00483158  | Yes                    | higher            | No               |                     40 |
| R184C     | deleterious       | pLDDT                   |   96.9287 |        96.6938 |      -0.2349 |       -3.8063 | 0.000343914 | Yes                    | lower             | No               |                     40 |
| R184C     | deleterious       | pLDDT min (Cα)          |   76.2472 |        70.8413 |      -5.406  |       -4.4844 | 3.98245e-05 | Yes                    | lower             | No               |                     40 |
| Y154H     | deleterious       | pair ipTM (protein–RNA) |    0.7652 |         0.7752 |       0.01   |        4.4173 | 3.38626e-05 | Yes                    | higher            | Yes *            |                     40 |
| Y154H     | deleterious       | PAE cross-mean (Å)      |    1.0755 |         1.0897 |       0.0142 |        1.5008 | 0.137579    | No                     | higher            | -                |                     40 |
| Y154H     | deleterious       | pLDDT                   |   96.9287 |        97.044  |       0.1153 |        1.7719 | 0.0810494   | No                     | higher            | -                |                     40 |
| Y154H     | deleterious       | pLDDT min (Cα)          |   76.2472 |        73.559  |      -2.6882 |       -1.894  | 0.0621094   | No                     | lower             | -                |                     40 |
| H212A     | deleterious       | pair ipTM (protein–RNA) |    0.7652 |         0.7703 |       0.005  |        1.892  | 0.062229    | No                     | higher            | -                |                     40 |
| H212A     | deleterious       | PAE cross-mean (Å)      |    1.0755 |         1.0473 |      -0.0283 |       -3.3553 | 0.00135696  | Yes                    | lower             | Yes *            |                     40 |
| H212A     | deleterious       | pLDDT                   |   96.9287 |        96.78   |      -0.1487 |       -1.8245 | 0.0719381   | No                     | lower             | -                |                     40 |
| H212A     | deleterious       | pLDDT min (Cα)          |   76.2472 |        78.776  |       2.5288 |        2.0904 | 0.0413622   | Yes                    | higher            | Yes *            |                     40 |
| D210E     | neutral           | pair ipTM (protein–RNA) |    0.7652 |         0.7693 |       0.004  |        1.5898 | 0.115917    | No                     | higher            | -                |                     40 |
| D210E     | neutral           | PAE cross-mean (Å)      |    1.0755 |         1.0825 |       0.007  |        0.7919 | 0.431138    | No                     | higher            | -                |                     40 |
| D210E     | neutral           | pLDDT                   |   96.9287 |        96.7838 |      -0.1449 |       -1.9498 | 0.0548313   | No                     | lower             | -                |                     40 |
| D210E     | neutral           | pLDDT min (Cα)          |   76.2472 |        74.547  |      -1.7002 |       -1.3457 | 0.183442    | No                     | lower             | -                |                     40 |
| S185N     | neutral           | pair ipTM (protein–RNA) |    0.7652 |         0.7708 |       0.0055 |        2.1642 | 0.0335068   | Yes                    | higher            | Yes *            |                     40 |
| S185N     | neutral           | PAE cross-mean (Å)      |    1.0755 |         1.1216 |       0.0461 |        3.4212 | 0.00106904  | Yes                    | higher            | No               |                     40 |
| S185N     | neutral           | pLDDT                   |   96.9287 |        97.0246 |       0.0959 |        1.4932 | 0.140306    | No                     | higher            | -                |                     40 |
| S185N     | neutral           | pLDDT min (Cα)          |   76.2472 |        72.9352 |      -3.312  |       -2.6397 | 0.0105892   | Yes                    | lower             | No               |                     40 |
| T187S     | neutral           | pair ipTM (protein–RNA) |    0.7652 |         0.7623 |      -0.003  |       -0.9519 | 0.344461    | No                     | lower             | -                |                     40 |
| T187S     | neutral           | PAE cross-mean (Å)      |    1.0755 |         1.1109 |       0.0354 |        3.1449 | 0.00236984  | Yes                    | higher            | No               |                     40 |
| T187S     | neutral           | pLDDT                   |   96.9287 |        96.7984 |      -0.1303 |       -1.8627 | 0.0664821   | No                     | lower             | -                |                     40 |
| T187S     | neutral           | pLDDT min (Cα)          |   76.2472 |        73.361  |      -2.8862 |       -2.3225 | 0.0237522   | Yes                    | lower             | No               |                     40 |
