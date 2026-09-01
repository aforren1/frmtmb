# Fit-time benchmarks (2026-08-31)

Single runs, wall-clock seconds, Windows / R 4.6.1. Script:
scratchpad bench.R (fit only, default controls each package).

| model                          | frm  | glmmTMB | lme4 |
|--------------------------------|------|---------|------|
| LMM n=180 (sleepstudy)         | 0.19 | 0.08    | 0.03 |
| LMM n=50k, 200 grp, (x\|g)     | 4.61 | 2.86    | 0.35 |
| Poisson GLMM n=100k, 500 grp   | 6.00 | 5.18    | 9.89 |
| NB2 + dispformula n=20k        | 1.31 | 0.71    | -    |
| ZI-poisson n=20k               | 0.47 | 0.58    | -    |

Reading:
- lme4 owns pure gaussian LMMs (profiled deviance + specialized sparse
  Cholesky); nothing Laplace-generic touches that.
- vs glmmTMB (the architectural peer): within ~2x everywhere, at parity
  or faster on GLMM/zi workloads. The gap on gaussian models is tape
  overhead plus no OpenMP parallel accumulation (RTMB limitation).
- frm beats lme4 on the large Poisson GLMM (glmer's AGQ/PIRLS scales
  worse than Laplace-with-AD here).
- Not measured: glmmTMB's parallel threads (off by default), warm-start
  refits, REML variants.

# Would optimParallel help on larger problems? (2026-09-01)

Verdict up front: **no, do not adopt optimParallel.** Best case measured
is 1.03-1.22x wall-clock, and that only when the PSOCK cluster and the
worker tapes already exist; paying for them inside the fit makes it
1.31x *slower* than the same optimizer run sequentially. The native
finite-difference mode is 2.3-2.5x slower than exact AD and does not
converge as well. The levers that do move InstEval are `profile = TRUE`
and `optimizer = "optim"`, both already in the package; everything else
is inside TMB's C++.

## Setup

Machine: Windows 11, R 4.6.1, 12 physical / 16 logical cores, frmtmb
0.24.0 (worktree wt-bench), optimParallel 1.0-3, lme4 2.0.6, RTMB 1.9,
TMB 1.9.25. Data: `lme4::InstEval`, n = 73421, s = 2972 students,
d = 1128 lecturers.

Models:
- `service`: `y ~ service + (1 | s) + (1 | d)` (the canonical model)
- `plain`:   `y ~ 1 + (1 | s) + (1 | d)`

Tape size for both: **5 outer parameters** (2 beta, 1 betad = log sigma,
2 theta = log sds; 4 for `plain`) and **4100 inner parameters**; inner
Hessian 4100 x 4100 with 77521 nonzeros.

Scripts (all under `dev/`, none of them touch package code):
`bench-op-baseline.R`, `bench-op-profile.R`, `bench-op-ceiling.R`,
`bench-op-extptr.R`, `bench-op-optimparallel.R`, `bench-op-tight.R`,
`bench-op-context.R`, plus `bench-setup.R` / `bench-env.R`.
Timing scripts were run one at a time. The package test suite was not
run: no package code changed.

## 1. Baseline

`frm(..., se = TRUE, control = frmtmb_control(verbose = 1))`, medians of
3 runs, seconds:

| stage      | service | plain |
|------------|---------|-------|
| parse      | 0.00    | 0.00  |
| frame      | 0.14    | 0.14  |
| tape       | 0.66    | 0.66  |
| optimize   | 22.78   | 11.78 |
| restart 1  | 0.44    | 3.34  |
| sdreport   | 3.84    | 3.29  |
| total      | 28.63   | 19.41 |

`service` costs materially more than `plain` in the optimizer (40 vs 26
nlminb iterations), so both were carried through the comparison.
Run-to-run spread on this machine is roughly +/-15% on the optimize
stage; treat every number below as two significant figures.

Everything outside `optimize` is fixed cost: parse + frame + tape is
~0.8 s and `sdreport` ~3.8 s. **`optimize` is 80% of the fit**, which is
what makes the optimParallel question worth asking at all.

### What is inside `optimize`

A custom optimizer that wraps `fn`/`gr` in timers and then calls
`nlminb` exactly as `run_optimizer()` does (`restarts = 0`, so nothing
is masked by the restart pass):

| quantity                          | service | plain |
|-----------------------------------|---------|-------|
| nlminb iterations                 | 40      | 26    |
| fn evaluations                    | 63      | 37    |
| gr evaluations                    | 41      | 27    |
| optimizer wall (s)                | 19.62   | 12.98 |
| time inside `obj$fn` (s)          | 11.99   | 7.72  |
| time inside `obj$gr` (s)          | 7.63    | 5.26  |
| nlminb's own overhead (s)         | 0.00    | 0.00  |
| median `obj$fn` call (s)          | 0.17    | 0.17  |
| median `obj$gr` call (s)          | 0.18    | 0.19  |

**100% of the optimize stage is objective evaluation.** nlminb itself
costs nothing measurable on a 5-parameter problem. The first `fn` call
of a fit is ~1.6 s (cold inner Newton solve); every later one is ~0.17 s
because the inner solve warm-starts from the previous point.

Standalone per-call costs at the optimum (medians of 20 calls at
jittered parameter vectors, `dev/bench-op-ceiling.R`):

| call                                  | seconds |
|---------------------------------------|---------|
| `obj$fn(x)` (inner Newton solve)       | 0.165   |
| `obj$gr(x)` right after `obj$fn(x)`    | 0.30    |
| `obj$gr(x)` at a fresh x               | 0.28    |

Note the third row: `obj$gr` costs the same whether or not `obj$fn` has
just run at the same point. TMB's gradient carries its own inner solve,
so splitting `fn` and `gr` across two processes does not split the work
in half; the `gr` worker still does everything the `fn` worker does,
plus the adjoint sweep.

## 2. Analysis before tooling

optimParallel is a parallel L-BFGS-B. Reading `optimParallel:::FGgenerator`:

- **With an exact `gr` supplied** it builds an `exprList` of length
  **2** and dispatches it with `parLapply`. Two tasks. Its headline
  "~p cores" speedup applies only to the `gr = NULL` path, where the
  gradient is finite-differenced into `2p+1` (central) or `p+1`
  (forward) tasks. We have exact AD gradients, so the concurrency
  ceiling is 2.
- The master evaluates nothing itself; `testFn` only inspects formals.
  So the honest wrapper must give the workers something to evaluate.
- `evalFG` caches on `identical(par, par_last)`, so an L-BFGS-B step
  that wants both `f` and `g` costs one parallel round.

Theoretical ceiling on the evaluation time, from the table above:
sequential pair = 0.165 + 0.30 = 0.465 s; parallel round =
max(0.165, 0.28) = 0.28 s; **ceiling 1.66x on `optimize` only**, i.e.
about 1.5x on the whole fit if communication were free. Nowhere near
"p cores", and that is before any of the practical problems below.

### RTMB tapes do not survive a PSOCK worker

Verified concretely (`dev/bench-op-extptr.R`). The tape object holds
external pointers:

```
external pointers reachable from obj (depth <= 4):
  obj$env$ADFun$ptr  -> <pointer: 0x00000244d59ce910>
  obj$env$ADGrad$ptr -> <pointer: 0x00000244d59ce670>
```

`clusterExport(cl, "obj")` transfers 13.9 MB and both pointers arrive
dead:

```
--- pointers as seen ON THE WORKER after clusterExport ---
  obj$env$ADFun$ptr  -> <pointer: (nil)>
  obj$env$ADGrad$ptr -> <pointer: (nil)>
```

Calling the transferred object does not raise an R error, it kills the
worker process:

```
--- obj$fn(obj$par) ON THE WORKER ---
worker process DIED. master-side error:
  error reading from connection
```

(The same call on the master returns 121579.0038.)

The only viable pattern is therefore: ship the **frame**, rebuild the
**tape** in each worker. That works and reproduces the master's value
exactly. Per-fit setup cost for 2 workers:

| step                            | seconds |
|---------------------------------|---------|
| `makePSOCKcluster(2)`           | 0.26    |
| `library(frmtmb)` on both       | 1.28    |
| ship the frame (5.3 MB) to both | 0.02    |
| rebuild the tape on both        | 0.83    |
| **total**                       | **2.39** |

Shipping the frame is free; the cost is R startup, package loading and
retaping. It is amortizable only across repeated fits that reuse the
same cluster and the same frame (a bootstrap or a profile-likelihood
sweep), not within one fit.

## 3. Measurement

`restarts = 0` everywhere so the optimizer counts are honest; medians of
3 runs; `service` model; end-to-end `frm()` wall clock including frame,
tape, cluster setup and tape rebuilds.

| mode                                         | wall (s) | optimize (s) | rounds / iters | max abs grad |
|----------------------------------------------|----------|--------------|----------------|--------------|
| nlminb, sequential, exact gr (package default) | 20.39  | 19.50 (fn 12.08 + gr 7.42) | 40 it, 64 fn, 41 gr | 0.0022 |
| optim L-BFGS-B, sequential, exact gr          | 9.26     | 8.27 (fn 4.38 + gr 3.89)   | 20 it, 20 fn, 20 gr | 0.165 |
| optimParallel, exact gr, 2 workers, cold cluster | 12.11 | 7.45 + 2.08 setup          | 20 rounds       | 0.165 |
| optimParallel, exact gr, 2 workers, warm cluster | 8.93  | 6.08                       | 20 rounds       | 0.165 |
| optimParallel, central FD, 11 workers (2p+1)  | 21.55    | 13.56 + 5.22 setup         | 20 rounds       | 0.164 |
| optimParallel, forward FD, 6 workers (p+1)    | 23.26    | 16.14 + 4.04 setup         | 34 rounds       | 68.8  |

All modes reach the same objective, 118865.3077.

Reading:

- **(a) exact gradient.** Against the same algorithm run sequentially
  (`optim` L-BFGS-B), optimParallel with a pre-built cluster is
  9.26 -> 8.93 s, **1.04x**. On the optimize stage alone it is 8.27 ->
  6.08 s, 1.36x, well short of the 1.66x ceiling because the `gr`
  worker repeats the inner solve and each round pays a `clusterExport`
  plus two `parLapply` round trips. Built honestly, cluster and tapes
  inside the fit, it is 12.11 s: **31% slower than sequential**. The
  optimize stage saves 2.19 s (8.27 -> 6.08) and the cold cluster costs
  2.08 s of measured setup, 3.18 s of extra wall once the first-round
  export to fresh workers is counted.
- **(b) finite differences (its native mode).** 11 workers, 20 rounds,
  21.55 s: **2.3x slower** than sequential exact AD, exactly as
  predicted - 11 evaluations per step across cores cannot beat 1 fn +
  1 gr on one core when the gradient is exact. Forward differences with
  6 workers is worse still (23.26 s) and lands at max abs gradient 68.8
  instead of 0.16, i.e. it does not actually solve the problem.
- **(c)** every mode ran; nothing was blocked outright, but only because
  the wrapper rebuilds the tapes. The naive route - export the tape -
  segfaults the worker (section 2).

### Controlling for convergence tolerance

L-BFGS-B stops at max abs gradient 0.165 where nlminb reaches 0.0022, so
part of the L-BFGS-B advantage is a looser stop. Re-running both at a
tight `factr` (`dev/bench-op-tight.R`, medians of 3):

| model   | optimizer                       | wall (s) | optimize (s) | rounds | max abs grad |
|---------|---------------------------------|----------|--------------|--------|--------------|
| service | nlminb (defaults)               | 20.82    | 19.66        | 40     | 0.0022 |
| service | L-BFGS-B seq, factr = 1e7       | 8.97     | 8.09         | 20     | 0.165  |
| service | optimParallel(2w), factr = 1e7  | 8.67     | 5.92         | 20     | 0.165  |
| service | L-BFGS-B seq, factr = 1e2       | 11.65    | 10.57        | 26     | 0.0116 |
| service | optimParallel(2w), factr = 1e2  | 9.97     | 7.56         | 26     | 0.0115 |
| plain   | nlminb (defaults)               | 13.16    | 12.09        | 26     | 0.0056 |
| plain   | L-BFGS-B seq, factr = 1e7       | 10.64    | 9.46         | 17     | 0.285  |
| plain   | optimParallel(2w), factr = 1e7  | 8.70     | 5.94         | 17     | 0.285  |
| plain   | L-BFGS-B seq, factr = 1e2       | 20.70    | 19.25        | 49     | 0.0034 |
| plain   | optimParallel(2w), factr = 1e2  | 18.06    | 15.40        | 45     | 0.0034 |

optimParallel's advantage over the identical sequential optimizer is
1.03x, 1.17x, 1.22x, 1.15x across these four pairs. The optimize stage
alone improves by 1.25-1.59x. Neither figure survives the 2.4 s of
per-fit cluster setup that a real `frm()` call would have to pay.

## 4. Where the time actually goes

`Rprof` over a whole `service` fit (`se = TRUE`):

```
                self.time self.pct
".Call"             7.315    96.38
".External"         0.090     1.19
```

and over 20 bare `obj$fn` calls: `.Call` 100%. After `MakeADFun` has
taped the R closure once, every evaluation is TMBad tape replay plus a
sparse Cholesky of the 4100 x 4100 / 77521-nonzero inner Hessian, all in
TMB's C++. There is no R-level time left to parallelize, and the inner
solve is a single sparse factorization that outer-loop parallelism
cannot touch.

Options that already exist, `service` / `plain`, medians of 3,
`se = FALSE`:

| configuration                        | service (s) | plain (s) |
|--------------------------------------|-------------|-----------|
| frm default (nlminb, restarts = 1)   | 21.50       | 17.22     |
| frm default + `se = TRUE`            | 23.35       | 24.22     |
| `optimizer = "optim"` (L-BFGS-B)     | 10.44       | 9.95      |
| `profile = TRUE`                     | 13.54       | 15.94     |
| `sparse_x = TRUE`                    | 20.95       | 20.18     |
| `REML = TRUE`                        | 12.06       | 15.22     |
| `lme4::lmer(REML = FALSE)`           | 4.42        | 4.61      |
| `lme4::lmer(REML = TRUE)`            | 5.20        | 5.10      |

frm and lmer agree on the objective to 4 decimals (118865.3077 for
`service`, 118888.8630 for `plain`), so this is a like-for-like
comparison of the same optimum.

- **lmer is 4.9x faster** on the canonical model (4.42 vs 21.50 s), as
  expected and as the
  2026-08-31 table already says: it profiles both beta and sigma out of
  the deviance, leaving a 2-parameter problem over a purpose-built
  sparse Cholesky, while frm hands nlminb 5 parameters over a generic
  Laplace tape. The point here is context, not a race.
- `sparse_x` does nothing: X is 73421 x 2.
- `profile = TRUE` is worth 1.6x on `service` with only two fixed
  effects, and should scale with the number of fixed effects.
- `optimizer = "optim"` is worth 2.1x here, but the tight-tolerance
  table shows it is model-dependent: on `plain` at `factr = 1e2` it is
  slower than nlminb (20.70 vs 13.16). Not a safe default change on this
  evidence.

## 5. Verdict

**Do not adopt optimParallel.**

1. With exact AD gradients its concurrency is 2, not p. Measured
   ceiling 1.66x on the optimize stage; measured reality 1.25-1.59x on
   the optimize stage and 1.03-1.22x end to end with a warm cluster.
2. RTMB tapes hold external pointers that deserialize to `(nil)` and
   crash the worker, so every worker must rebuild its own tape. That
   costs 2.4 s per fit on InstEval (mostly R startup and
   `library(frmtmb)`), which wipes out the gain: cold-cluster
   optimParallel is 12.11 s against 9.26 s sequential.
3. Its native finite-difference mode is 2.3-2.5x slower than exact AD
   and converges worse. Exact gradients are the whole point of the
   RTMB design; giving them up to parallelize is a strict loss.
4. Adopting it would also mean owning a `parallel` dependency, cluster
   lifecycle, and a second optimizer path, for a gain that is inside the
   machine's timing noise.

What would help on large problems, in order of measured effect:

- **Nothing in the outer loop.** 100% of `optimize` is inside
  `obj$fn`/`obj$gr`, and 96-99% of that is `.Call` into TMB's C++. The
  real lever is inside the inner sparse Cholesky, which is TMB's
  domain; RTMB does not expose TMB's OpenMP `parallel_accumulator`, and
  that is the same limitation the 2026-08-31 notes already record.
- **`profile = TRUE`** for models with many fixed effects: 1.6x here
  with p = 2, and it removes fixed effects from the outer problem
  entirely, so the benefit grows with p. This is the one that scales.
- **Fewer outer iterations**: better starting values, or `autoscale`,
  attack the 40-iteration count directly. Each iteration is ~0.5 s of
  irreducible C++, so halving iterations halves the fit.
- **Not `sparse_x`** unless the design matrix is genuinely sparse and
  wide; it is a no-op here.
- **Not a different optimizer as a default.** `optimizer = "optim"` is
  2.1x on the canonical model but slower on `plain` at matched
  tolerance. Worth documenting as a per-model option, not as a change
  of default.
- For repeated fits over one dataset (bootstrap, profile likelihood,
  `allfit`), a persistent cluster where each worker rebuilds the tape
  once and then runs **whole fits** amortizes the 2.4 s properly. That
  is parallelism at the fit level, which is where the concurrency
  actually is, and it does not need optimParallel.

## Note: RTMBp (2026-09-01, unmeasured)

kaskr packages an experimental PARALLEL build of RTMB as 'RTMBp'
(r-universe only, not CRAN: install.packages('RTMBp', repos =
c('https://kaskr.r-universe.dev', ...))). This targets exactly the
gap the optimParallel benchmark identified: InstEval's optimize
stage is 100% taped-objective evaluation, and TMB's OpenMP
parallelizes the objective accumulation (the inner Cholesky stays
serial). A follow-up measurement would swap RTMB -> RTMBp on the
InstEval protocol above and compare fn/gr per-call costs across
thread counts. Caveats before adopting anything: RTMBp is
experimental and off-CRAN, so it could at most be an optional
enhancement (namespace-conditional), never a dependency; and the
"only real TMB-vs-RTMB gap is OpenMP" line in SPEC.md predates it.
