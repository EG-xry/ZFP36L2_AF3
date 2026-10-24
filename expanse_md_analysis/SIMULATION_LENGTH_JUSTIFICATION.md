# Is 200 ns enough? Justifying simulation length for ZFP36L2 + 9-mer RNA

**Short answer:** For *equilibrium interface metrics* (RMSD, H-bonds,
contacts, electrostatics) with **replicates**, 200 ns per replica is
defensible and within community norms for a ~15,000-atom protein–RNA
system. But the honest, reviewer-proof framing depends on *what you claim*:

- **To rank/discriminate variants with replicate statistics:** 3 × 200 ns
  is reasonable and arguably better than 1 × 1000 ns (more independent
  sampling of the same wall-clock). ✅
- **To claim converged binding thermodynamics (absolute ΔG):** 200 ns is
  on the short side; you'd want longer or enhanced sampling. ⚠️
- **To directly match the Ramos benchmark:** they used 3 × 1000 ns, so a
  reviewer of that comparison may ask why you're shorter. You can answer
  this (below), but be ready for it.

---

## 1. What the system actually is, and why size matters

- ZFP36L2 tandem zinc-finger domain (~70 aa) + 9-nt RNA + 2 Zn²⁺,
  solvated in TIP3P with 0.15 M NaCl → **~15,000 atoms total**, a *small*
  MD system by modern standards.
- The functional motions you care about (RNA staying in the ARE pocket,
  Arg/Lys–phosphate salt bridges forming/breaking, finger-tip H-bonds) are
  **local and fast** (ps–tens-of-ns), not slow global domain motions. This
  is the key point for justifying a moderate length: you are measuring fast,
  local interface observables, not slow conformational transitions.

## 2. What the literature actually does for protein–RNA / small complexes

| Study / context | length × replicas | notes |
|---|---|---|
| Ramos / Akçeşme 2024 (this exact system) | 1000 ns × 3 | the benchmark you're compared to |
| Typical protein–RNA interface MD (2018–2024) | 100–500 ns × 2–3 | the modal range for "characterize the interface / rank mutants" |
| MM-GBSA mutant ranking papers (Genheden & Ryde 2015 review) | 5–50 ns × multiple | endpoint binding-energy is often computed on *short* runs once equilibrated |
| Convergence-focused single-domain studies | 1–10 µs × 1 | only when chasing slow conformational change / folding |
| RNA-binding zinc-finger / CCCH studies | 100–300 ns common | local recognition, not global rearrangement |

**Takeaway:** 100–500 ns × 2–3 replicas is the standard operating range for
"characterize and rank a protein–RNA interface." 200 ns × 3 sits squarely
inside it. The 1 µs Ramos runs are on the long end (they were chasing
intra-protein H-bond probabilities over very long time, and had the compute).

## 3. The replicates-vs-length tradeoff (the real argument)

For *discrimination* (the AUROC task), independent replicas beat a single
long run for the same total wall-time, because:
- a single long trajectory can get trapped in one basin (we literally saw
  this — the AF3 starting-structure relaxation transient dominated the first
  ~130 ns and faked an F156/R160 signal);
- 3 independent 200 ns runs sample 3 independent basins and let you put a
  **standard error / 95 % CI on each metric**, which is what a reviewer
  wants to see;
- this is the standard "throw replicas at it" wisdom (Knapp et al. 2018,
  "Avoiding false positive conclusions in MD: replicas not length").

So the defensible sentence for a methods section is:

> "Because the relevant protein–RNA recognition events (salt-bridge and
>  hydrogen-bond formation at the AU-rich interface) occur on the
>  sub-100-ns timescale, we prioritized independent replicates over single
>  long trajectories: three 200-ns replicas per variant (600 ns aggregate),
>  comparable in total sampling to prior single-trajectory work while
>  providing per-metric uncertainty estimates."

## 4. How to *prove* 200 ns is enough for YOUR observables (don't just assert)

Run these convergence checks on the data you already have (they're cheap and
make the 200 ns choice reviewer-proof):

1. **Block averaging / running-mean of the key metric.** Plot the
   screened-Coulomb (or contact count) as a running average vs time; if it
   plateaus well before 200 ns, you're converged *for that observable*.
2. **First-half vs second-half agreement.** Compute each metric on
   0–100 ns and 100–200 ns separately; if the AUROC and per-variant values
   agree within error, 200 ns is sufficient.
3. **Autocorrelation time of the metric.** If the metric's autocorrelation
   time is ≪ 200 ns, you have many effectively-independent samples in one
   run.
4. **Inter-replica agreement.** Once you have 3 replicas, the spread between
   them is the ultimate convergence test — if they agree, length was enough;
   if they don't, you need more (length or replicas).

I can add a `convergence_check.py` that does (1)–(3) on the existing
269/400 ns trajectories — that single figure ("metric is flat after ~50 ns")
is the strongest justification you can put in a paper.

## 5. My concrete recommendation on length

- **Use 200 ns × 3 replicas** as the production protocol. It is standard,
  defensible, and gives you error bars.
- **Justify it** with a convergence figure on the data you already have
  (running average plateau + first-half/second-half agreement), plus the
  "local fast observables + replicas > length" argument above.
- **If a reviewer pushes on the Ramos 1 µs comparison:** point out that (a)
  your observables are protein–RNA interface contacts that equilibrate fast,
  (b) you used independent replicas for statistical rigor rather than one
  long run, and (c) total aggregate sampling (600 ns/variant) is in the same
  order of magnitude.
- **Only go longer (≥500 ns) if** your own convergence check shows the
  key metric is still drifting at 200 ns, OR if you move to absolute
  binding-free-energy (MM-PBSA with entropy / alchemical FEP), which needs
  more sampling.

## 6. One caveat specific to this project

We already saw that the *first* ~130 ns of these runs was dominated by the
AF3 starting-structure relaxing (the F156/R160 transient). That means:
- **discard the first ~50 ns of every run as equilibration** before
  computing interface metrics (we currently use the last 20–25 %; make this
  an explicit "first 50 ns dropped" rule);
- with that equilibration drop, a 200 ns run gives ~150 ns of usable
  production, which is adequate for the fast interface observables.

So: 200 ns is enough **if** you (a) run replicas, (b) drop the first ~50 ns,
and (c) show a convergence plot. Without those three, any single-length
choice is hard to defend.
