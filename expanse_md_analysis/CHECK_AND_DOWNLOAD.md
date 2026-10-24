# Expanse: check completion + download everything back locally

Host / account details (from `anvil_md/expanse/EXPANSE_README.md`):
- **SSH:** `YOUR_USER@login.expanse.sdsc.edu`
- **Remote project dir:** `~/zfp36l2_md`
- **Remote outputs:** `~/zfp36l2_md/outputs/<VARIANT>/`
- **Account / partition:** `bio260224` / `gpu-shared`

The 10 variant folders use these names on Expanse (note the original
`C208A` name — it was *renamed to C212A only locally*, so on Expanse it is
still `C208A`):
`WT  Y154H  R160W  R184C  G204D  C206F  S185N  T187S  D210E  C208A`

---

## STEP 1 — SSH in

```bash
ssh YOUR_USER@login.expanse.sdsc.edu
```

## STEP 2 — Is anything still running in the queue?

```bash
squeue -u $USER -o "%.10i %.12j %.10P %.8T %.10M %.10l %.6D %R"
```
- If this prints **no job rows**, nothing is queued/running → all jobs have
  exited (either finished or hit wall-time and checkpointed).
- If you see rows in state `R` (running) or `PD` (pending), they are not done yet.

## STEP 3 — One-line completion check across all 10 variants

Paste this whole block on Expanse — it prints ns done / target / status for
every variant and a final "X/10 completed" tally:

```bash
cd ~/zfp36l2_md/outputs
printf "%-8s %-13s %10s %10s\n" VARIANT STATUS ns_done target; \
echo "------------------------------------------------"; \
done_ct=0; \
for V in WT Y154H R160W R184C G204D C206F S185N T187S D210E C208A; do \
  P="$V/progress.txt"; \
  if [ -f "$P" ]; then \
    st=$(grep '^status=' "$P" | cut -d= -f2); \
    ns=$(grep '^ns_completed=' "$P" | cut -d= -f2); \
    tg=$(grep '^target_ns=' "$P" | cut -d= -f2); \
    printf "%-8s %-13s %10s %10s\n" "$V" "$st" "$ns" "$tg"; \
    [ "$st" = "completed" ] && done_ct=$((done_ct+1)); \
  else \
    printf "%-8s %-13s %10s %10s\n" "$V" "NO_OUTPUT" "-" "-"; \
  fi; \
done; \
echo "------------------------------------------------"; \
echo "completed: $done_ct / 10"
```

**Interpretation**
- `status=completed` and `ns_done == target_ns` (400) → that variant is done.
- `status=checkpointed` → it hit wall-time before 400 ns; it needs to be
  re-submitted with `--resume` to finish (it will continue from its
  checkpoint).
- `completed: 10 / 10` at the bottom → everything is done; go to STEP 5.

## STEP 4 — (only if some are still `checkpointed`) resume them

```bash
cd ~/zfp36l2_md
# Resume just the ones that aren't at 400 ns yet. Example:
bash scripts/submit_all_expanse.sh 400        # resubmits all; --resume continues from checkpoints
# (or submit individually)
# bash scripts/submit_variant_expanse.sh R184C 400 --resume
squeue -u $USER
```
Then re-run STEP 3 later to confirm they reached `completed`.

## STEP 5 — Check total size before downloading (optional)

```bash
du -sh ~/zfp36l2_md/outputs
du -sh ~/zfp36l2_md/outputs/* | sort -h
```
Each variant's `production.dcd` at 400 ns is ~1 GB, so the full set is
~8–11 GB. Make sure you have the disk space locally.

---

## STEP 6 — Download EVERYTHING back to your Mac

Run this **from your Mac terminal** (not from inside the Expanse SSH
session). `rsync` is preferred over `scp` because it resumes partial
transfers and skips files you already have:

```bash
rsync -avz --progress \
  YOUR_USER@login.expanse.sdsc.edu:'~/zfp36l2_md/outputs/' \
  "/Users/eric/Desktop/ZFP36L2 AF3/expanse_md_results/"
```

Notes:
- The trailing slash on `outputs/` means "copy the *contents* of outputs/
  into expanse_md_results/" so you get `expanse_md_results/WT/…`,
  `expanse_md_results/C206F/…`, etc.
- `-a` preserves timestamps, `-z` compresses on the wire, `--progress`
  shows per-file progress, and rsync will **resume** if the connection
  drops — just re-run the same command.
- The variant that Expanse calls `C208A` will land in
  `expanse_md_results/C208A/`. Since you renamed it to `C212A` locally,
  after the download run:
  ```bash
  cd "/Users/eric/Desktop/ZFP36L2 AF3/expanse_md_results"
  # merge the freshly-downloaded C208A into your renamed C212A folder:
  rsync -av C208A/ C212A/ && rm -rf C208A
  # and fix the variant= line in its progress.txt:
  sed -i '' 's/^variant=C208A/variant=C212A/' C212A/progress.txt
  ```

### If you only want the small files first (fast sanity check)
To grab just the logs and progress files (a few MB) before the big DCDs:

```bash
rsync -avz --progress \
  --include='*/' --include='md.log' --include='progress.txt' \
  --include='final_frame.pdb' --include='equilibrated.pdb' --exclude='*' \
  YOUR_USER@login.expanse.sdsc.edu:'~/zfp36l2_md/outputs/' \
  "/Users/eric/Desktop/ZFP36L2 AF3/expanse_md_results/"
```

### Fallback with plain scp (no resume)
```bash
scp -r YOUR_USER@login.expanse.sdsc.edu:'~/zfp36l2_md/outputs/' \
  "/Users/eric/Desktop/ZFP36L2 AF3/expanse_md_results_new/"
```

---

## STEP 7 — Re-run local analysis with the fresh 400 ns data

```bash
cd "/Users/eric/Desktop/ZFP36L2 AF3"
# Widen the common window to 400 ns now that everything is 400 ns:
sed -i '' 's/^COMMON_WINDOW_NS = .*/COMMON_WINDOW_NS = 400.0  # all variants finished/' \
  expanse_md_analysis/plot_drift_and_hbonds.py
/Users/eric/miniconda3/bin/python3 expanse_md_analysis/run_full_analysis.py
/Users/eric/miniconda3/bin/python3 expanse_md_analysis/make_auroc_figure.py
```

That regenerates every CSV, the AUROC chart, and the H-bond / drift figures
on the full 400 ns common window.
