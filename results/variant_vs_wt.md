# ZFP36L2 AF3 — Welch's t-test: variant vs WT (model_0 only; counts by variant: WT=8, C206F=8, C208A=8, G204D=8, R160W=8, R160E=5, R184C=8, Y154H=8, H212A=5, D210E=8, S185N=8, T187S=8)

## Methodology

For each metric, an independent two-sample Welch's t-test (unequal
variance) compares each variant's available values against WT's
available values for the same analysis set. Significance: α = 0.05.
For hypothesis-testing only, each group is normalized to assumed n=8
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
| C206F     | deleterious       | pair ipTM (protein–RNA) |    0.7775 |         0.7212 |      -0.0563 |      -16.6716 | 4.03099e-09 | Yes                    | lower             | No               |                      8 |
| C206F     | deleterious       | PAE cross-mean (Å)      |    1.0506 |         1.2219 |       0.1712 |        6.2592 | 7.84877e-05 | Yes                    | higher            | No               |                      8 |
| C206F     | deleterious       | pLDDT                   |   96.9882 |        95.5677 |      -1.4205 |      -11.7744 | 5.08674e-08 | Yes                    | lower             | No               |                      8 |
| C206F     | deleterious       | pLDDT min (Cα)          |   72.47   |        66.1012 |      -6.3688 |       -3.8857 | 0.00291496  | Yes                    | lower             | No               |                      8 |
| C208A     | deleterious       | pair ipTM (protein–RNA) |    0.7775 |         0.7588 |      -0.0188 |       -5.5572 | 0.000174634 | Yes                    | lower             | No               |                      8 |
| C208A     | deleterious       | PAE cross-mean (Å)      |    1.0506 |         1.1612 |       0.1106 |        5.3482 | 0.000132111 | Yes                    | higher            | No               |                      8 |
| C208A     | deleterious       | pLDDT                   |   96.9882 |        97.1424 |       0.1542 |        1.8954 | 0.0814924   | No                     | higher            | -                |                      8 |
| C208A     | deleterious       | pLDDT min (Cα)          |   72.47   |        73.2013 |       0.7313 |        0.3555 | 0.727517    | No                     | higher            | -                |                      8 |
| G204D     | deleterious       | pair ipTM (protein–RNA) |    0.7775 |         0.7788 |       0.0012 |        0.4472 | 0.662229    | No                     | higher            | -                |                      8 |
| G204D     | deleterious       | PAE cross-mean (Å)      |    1.0506 |         1.0956 |       0.045  |        1.7849 | 0.101631    | No                     | higher            | -                |                      8 |
| G204D     | deleterious       | pLDDT                   |   96.9882 |        96.9344 |      -0.0538 |       -0.7298 | 0.482463    | No                     | lower             | -                |                      8 |
| G204D     | deleterious       | pLDDT min (Cα)          |   72.47   |        70.1012 |      -2.3688 |       -1.5283 | 0.162594    | No                     | lower             | -                |                      8 |
| R160W     | deleterious       | pair ipTM (protein–RNA) |    0.7775 |         0.7637 |      -0.0138 |       -2.1282 | 0.0661639   | No                     | lower             | -                |                      8 |
| R160W     | deleterious       | PAE cross-mean (Å)      |    1.0506 |         1.2663 |       0.2156 |        4.3131 | 0.00263025  | Yes                    | higher            | No               |                      8 |
| R160W     | deleterious       | pLDDT                   |   96.9882 |        97.0879 |       0.0996 |        1.1133 | 0.284599    | No                     | higher            | -                |                      8 |
| R160W     | deleterious       | pLDDT min (Cα)          |   72.47   |        74.2262 |       1.7562 |        0.8544 | 0.407274    | No                     | higher            | -                |                      8 |
| R160E     | deleterious       | pair ipTM (protein–RNA) |    0.7775 |         0.775  |      -0.0025 |       -1      | 0.334622    | No                     | lower             | -                |                      8 |
| R160E     | deleterious       | PAE cross-mean (Å)      |    1.0506 |         1.0719 |       0.0212 |        1.4475 | 0.173894    | No                     | higher            | -                |                      8 |
| R160E     | deleterious       | pLDDT                   |   96.9882 |        97.0399 |       0.0517 |        0.7259 | 0.486711    | No                     | higher            | -                |                      8 |
| R160E     | deleterious       | pLDDT min (Cα)          |   72.47   |        72.4875 |       0.0175 |        0.0111 | 0.991372    | No                     | higher            | -                |                      8 |
| R184C     | deleterious       | pair ipTM (protein–RNA) |    0.7775 |         0.77   |      -0.0075 |       -2.3932 | 0.0345847   | Yes                    | lower             | No               |                      8 |
| R184C     | deleterious       | PAE cross-mean (Å)      |    1.0506 |         1.0681 |       0.0175 |        0.9937 | 0.337242    | No                     | higher            | -                |                      8 |
| R184C     | deleterious       | pLDDT                   |   96.9882 |        96.7951 |      -0.1931 |       -1.9497 | 0.0716888   | No                     | lower             | -                |                      8 |
| R184C     | deleterious       | pLDDT min (Cα)          |   72.47   |        70.5763 |      -1.8937 |       -1.0624 | 0.308328    | No                     | lower             | -                |                      8 |
| Y154H     | deleterious       | pair ipTM (protein–RNA) |    0.7775 |         0.785  |       0.0075 |        2.3932 | 0.0345847   | Yes                    | higher            | Yes *            |                      8 |
| Y154H     | deleterious       | PAE cross-mean (Å)      |    1.0506 |         1.0656 |       0.015  |        0.954  | 0.35727     | No                     | higher            | -                |                      8 |
| Y154H     | deleterious       | pLDDT                   |   96.9882 |        97.0031 |       0.0149 |        0.1916 | 0.851459    | No                     | higher            | -                |                      8 |
| Y154H     | deleterious       | pLDDT min (Cα)          |   72.47   |        70.6037 |      -1.8663 |       -1.075  | 0.303906    | No                     | lower             | -                |                      8 |
| H212A     | deleterious       | pair ipTM (protein–RNA) |    0.7775 |         0.7775 |       0      |        0      | 1           | No                     | equal             | -                |                      8 |
| H212A     | deleterious       | PAE cross-mean (Å)      |    1.0506 |         1.0375 |      -0.0131 |       -0.6907 | 0.501278    | No                     | lower             | -                |                      8 |
| H212A     | deleterious       | pLDDT                   |   96.9882 |        96.8816 |      -0.1066 |       -1.0395 | 0.316598    | No                     | lower             | -                |                      8 |
| H212A     | deleterious       | pLDDT min (Cα)          |   72.47   |        82.65   |      10.18   |        4.8668 | 0.00024941  | Yes                    | higher            | Yes *            |                      8 |
| D210E     | neutral           | pair ipTM (protein–RNA) |    0.7775 |         0.7763 |      -0.0012 |       -0.3444 | 0.737403    | No                     | lower             | -                |                      8 |
| D210E     | neutral           | PAE cross-mean (Å)      |    1.0506 |         1.0806 |       0.03   |        1.5121 | 0.153751    | No                     | higher            | -                |                      8 |
| D210E     | neutral           | pLDDT                   |   96.9882 |        96.9536 |      -0.0347 |       -0.2194 | 0.830767    | No                     | lower             | -                |                      8 |
| D210E     | neutral           | pLDDT min (Cα)          |   72.47   |        77.1325 |       4.6625 |        2.127  | 0.0518674   | No                     | higher            | -                |                      8 |
| S185N     | neutral           | pair ipTM (protein–RNA) |    0.7775 |         0.7812 |       0.0037 |        0.9697 | 0.355259    | No                     | higher            | -                |                      8 |
| S185N     | neutral           | PAE cross-mean (Å)      |    1.0506 |         1.0788 |       0.0281 |        1.283  | 0.222892    | No                     | higher            | -                |                      8 |
| S185N     | neutral           | pLDDT                   |   96.9882 |        97.1186 |       0.1304 |        1.5073 | 0.154801    | No                     | higher            | -                |                      8 |
| S185N     | neutral           | pLDDT min (Cα)          |   72.47   |        73.605  |       1.135  |        0.5209 | 0.610615    | No                     | higher            | -                |                      8 |
| T187S     | neutral           | pair ipTM (protein–RNA) |    0.7775 |         0.77   |      -0.0075 |       -1.8209 | 0.100084    | No                     | lower             | -                |                      8 |
| T187S     | neutral           | PAE cross-mean (Å)      |    1.0506 |         1.0931 |       0.0425 |        2.223  | 0.0436121   | Yes                    | higher            | No               |                      8 |
| T187S     | neutral           | pLDDT                   |   96.9882 |        96.9991 |       0.0108 |        0.1168 | 0.908718    | No                     | higher            | -                |                      8 |
| T187S     | neutral           | pLDDT min (Cα)          |   72.47   |        74.5188 |       2.0488 |        1.2375 | 0.24286     | No                     | higher            | -                |                      8 |
