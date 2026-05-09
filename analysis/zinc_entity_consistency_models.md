# Zinc runs: entity consistency per model

- `rna_seq` is extracted from the model residue order; compare across models for 5'→3' consistency.
- `n_zn_atoms` should be stable across models if Zn is defined consistently.

| ComplexBase                         |   Model |   n_protein_res |   n_rna_res |   n_zn_atoms | rna_seq   | rna_seq_rev   | protein_chain_ids   | rna_chain_ids   | zn_chain_ids   |
|:------------------------------------|--------:|----------------:|------------:|-------------:|:----------|:--------------|:--------------------|:----------------|:---------------|
| fold_zfp36l2_c206f_nonamer_wt_2zn   |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_c206f_nonamer_wt_2zn   |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_c206f_nonamer_wt_2zn   |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_c206f_nonamer_wt_2zn   |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_c206f_nonamer_wt_2zn   |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_d210e_nonamer_wt_2zn   |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_d210e_nonamer_wt_2zn   |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_d210e_nonamer_wt_2zn   |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_d210e_nonamer_wt_2zn   |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_d210e_nonamer_wt_2zn   |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_g204d_nonamer_wt_2zn   |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_g204d_nonamer_wt_2zn   |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_g204d_nonamer_wt_2zn   |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_g204d_nonamer_wt_2zn   |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_g204d_nonamer_wt_2zn   |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r160w_nonamer_wt_2zn   |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r160w_nonamer_wt_2zn   |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r160w_nonamer_wt_2zn   |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r160w_nonamer_wt_2zn   |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r160w_nonamer_wt_2zn   |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r184c_nonamer_wt_2zn   |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r184c_nonamer_wt_2zn   |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r184c_nonamer_wt_2zn   |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r184c_nonamer_wt_2zn   |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_r184c_nonamer_wt_2zn   |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_s185n_nonamer_wt_2zn   |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_s185n_nonamer_wt_2zn   |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_s185n_nonamer_wt_2zn   |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_s185n_nonamer_wt_2zn   |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_s185n_nonamer_wt_2zn   |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_t187s_nonamer_wt_2zn   |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_t187s_nonamer_wt_2zn   |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_t187s_nonamer_wt_2zn   |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_t187s_nonamer_wt_2zn   |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_t187s_nonamer_wt_2zn   |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_nonbind_2zn |       0 |              70 |           9 |            2 | UACAAACAU | UACAAACAU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_nonbind_2zn |       1 |              70 |           9 |            2 | UACAAACAU | UACAAACAU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_nonbind_2zn |       2 |              70 |           9 |            2 | UACAAACAU | UACAAACAU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_nonbind_2zn |       3 |              70 |           9 |            2 | UACAAACAU | UACAAACAU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_nonbind_2zn |       4 |              70 |           9 |            2 | UACAAACAU | UACAAACAU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_wt_2zn      |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_wt_2zn      |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_wt_2zn      |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_wt_2zn      |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_wt_nonamer_wt_2zn      |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_y154h_nonamer_wt_2zn   |       0 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_y154h_nonamer_wt_2zn   |       1 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_y154h_nonamer_wt_2zn   |       2 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_y154h_nonamer_wt_2zn   |       3 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |
| fold_zfp36l2_y154h_nonamer_wt_2zn   |       4 |              70 |           9 |            2 | UUAUUUAUU | UUAUUUAUU     | A                   | D               | B,C            |