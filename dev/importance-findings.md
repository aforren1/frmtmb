# wt-islap findings (rewritten after a scratchpad wipe; keep appending)

## ENVIRONMENT HAZARD (happened twice)
1. The SHARED user library C:/Users/adf44/AppData/Local/R/win-library/4.6
   was emptied by something outside this worktree: 350 package
   directories still present, 212 with no DESCRIPTION and no contents
   (RTMB, Matrix, testthat, pkgload, lme4, GLMMadaptive all gone). No R
   process was running.
2. The SCRATCHPAD itself was then partly wiped, taking is-probe1.R,
   is-proto.R, is-measure*.R, is-smoke.R, is-dbg*.R, suite-49*.R and
   truncating this file.
Mitigation: a PRIVATE library at <scratchpad>/is-lib (228 packages,
CRAN binaries, is-install.R). Every command sets R_LIBS_USER to it.
The system library still supplies Matrix, mgcv, MASS, nlme, survival.
Numbers below are transcribed from the runs, not re-derivable from
files that no longer exist; re-measure before quoting anything.

## RTMB capability probe
- `S %*% RTMB::matrix(ll, n, N)` with S a Matrix sparse indicator works
  and gives a G x N advector matrix. No kron fallback needed.
- base `rowSums` DROPS the advector class (tape ends in a bare complex
  vector, "Invalid argument to 'advector'"). RTMB exports rowSums as an
  S4 generic: `RTMB::rowSums` is required. Same trap for colSums.
- The sparse product returns a Matrix object in the NUMERIC path, and a
  Matrix loses its dimensions through exp(); convert with as.matrix()
  when the result is not an advector.
- `rep(advector, N)` and `v[rep_idx]` both work.
- `obj$env$spHess(par, random = TRUE)` on a dim-2 block is exactly
  block diagonal, 2x2 per level, level-major. Confirms frame ordering.
- No R-side tape node count exists (env$ADFun is ptr+DLL, the ptr an
  externalptr with no methods). Tape size measured as process RSS delta
  via ps::ps_memory_info().

## Standalone prototype on the probe design (60 groups of 8 Bernoulli,
## correlated intercept and slope, seed 11)
| quantity | Laplace | IS round 1 | IS round 2 | GLMMadaptive nAGQ 15/25 |
|---|---|---|---|---|
| nll | 299.3017 | 297.9585 | 297.9488 | 297.9185 |
| intercept | -0.6915 | -0.6989 | -0.7026 | -0.7047 |
| slope | 0.3592 | 0.3557 | 0.3578 | 0.3617 |
N = 2000 antithetic, seed 7.

## Antithetic pairs: KEEP (sd of the objective over 8 seeds)
| N | plain | antithetic | variance ratio |
|---|---|---|---|
| 50 | 0.2669 | 0.1557 | 2.9x |
| 200 | 0.1007 | 0.0711 | 2.0x |
| 1000 | 0.0385 | 0.0250 | 2.4x |
| 2000 | 0.0330 | 0.0184 | 3.2x |

## ESS at each design's own Laplace optimum (N=1000 antithetic)
| design | q | G | n | min ESS/N | median | MCSE(nll) |
|---|---|---|---|---|---|---|
| probe (bin, 8/grp, corr slope) | 2 | 60 | 480 | 0.935 | 0.984 | 0.0351 |
| scalar intercept (bin, 4/grp) | 1 | 100 | 400 | 0.973 | 0.994 | 0.0275 |
| tiny+wide (bin, 3/grp, sd 2.0) | 1 | 40 | 120 | 0.429 | 0.868 | 0.0911 |
| gaussian | 1 | 50 | 300 | 1.000 | 1.000 | 0.0000 |

## ESS with the proposal frozen away from the optimum (probe, N=2000)
| log-sd shift | min ESS/N | median ESS/N | MCSE |
|---|---|---|---|
| 0.0 | 0.9334 | 0.9824 | 0.025 |
| 0.5 | 0.0089 | 0.2327 | 0.49 |
| 1.0 | 0.0027 | 0.0401 | 1.26 |
| 2.0 | 0.0022 | 0.0259 | 1.62 |

### THRESHOLD: warn when min ESS/N < 0.25 (imp_ess_floor)
The hardest LEGITIMATE design measured (40 groups of 3 binary rows, sd
2.0) sits at 0.429 at its own optimum, so 0.25 clears every valid fit
with a 1.7x margin. A proposal displaced by half a log-SD, already
worth 0.49 in nll, drops to 0.0089, thirty times below the threshold.
Two orders of magnitude separate the regimes.

## Tape cost vs Laplace (probe, n = 480, RSS delta)
| tape | RSS | gr | build |
|---|---|---|---|
| Laplace | +2.6 MB | <1 ms | - |
| IS N=50 | +14.6 MB | 4 ms | 0.09 s |
| IS N=200 | +47.9 MB | 2 ms | 0.42 s |
| IS N=1000 | +153.8 MB | 44 ms | 2.14 s |
Linear in N, about 320 bytes per (row x draw). fn timings were below
the clock resolution; re-measure with more repetitions.

## Per-level prior: PRECISION EXTRACTION (imp_prior_terms)
covstruct_registry[[k]]$nll returns the log density summed over ALL
levels; the per-group integrand needs one level's worth. Rather than
reimplement each structure, evaluate the registry's OWN density on a
one-level block at zero and at the canonical unit patterns and
difference out the normalizer and the precision matrix:
  dg[j]  = l(e_j) - l(0)          = -P[j,j]/2
  c[j,k] = l(e_j+e_k) - l(0) - dg[j] - dg[k] = -P[j,k]
  log p(u) = l(0) + sum_j dg[j] u_j^2 + sum_{j<k} c[j,k] u_j u_k
q + q(q+1)/2 + 1 small evaluations, taped once. Needs the density to be
Gaussian in the level's coefficients and independent between levels,
which is exactly what imp_covstructs whitelists.

## THE ANCHOR PROPERTY (measured in the package, decisive)
The IS objective evaluated AT ITS OWN ANCHOR equals the Laplace
objective there to 0.000e+00 for a gaussian response, at N = 2, 10,
200, 1000, over 12 seeds, with zero variance. Away from the anchor the
estimator is unbiased and its sd falls as N^-0.5:
| theta shift | N=10 sd | N=100 sd | N=1000 sd |
|---|---|---|---|
| 0.00 | 0 | 0 | 0 |
| 0.05 | 4.8e-2 | 2.0e-2 | 5.8e-3 |
| 0.20 | 2.2e-1 | 9.1e-2 | 2.5e-2 |
Bias at the shifted points is within about 2.5 standard errors of zero
over 12 seeds.

### CONSEQUENCE, and the design decision it forced
The optimizer minimizes truth(theta) + error(theta) where error is
EXACTLY zero at the anchor and grows away from it. So it is pushed
outward and stops at an optimistic value: on the gaussian design the
single-round optimum sat 1.97e-1 (N=2), 1.32e-2 (N=10), 6.98e-3
(N=200) BELOW the exact Laplace optimum, while the exact nll AT that
point was the same amount ABOVE it. No draw count removes this, and
refreezing alone does not either, because the displacement recurs each
round.
FIX (implemented): after the round loop converges, freeze ONE more
proposal at the final estimate and report the objective and the ESS
from that tape. At the anchor the weight variance is minimal (zero for
a gaussian response at any N), so the reported log-likelihood is the
honest one rather than the optimizer's. Costs one extra tape.
Residual: the PARAMETERS still drift by the noise-exploitation amount,
which is O(N^-1/2) and was about 0.1 of a standard error at N=200 on
the gaussian design. Report it; do not hide it.

### What test (a) can therefore demand
(a1) the identity, exact for ANY N: the corrected objective at its own
     anchor equals the Laplace objective there to machine precision,
     and every weight is equal. This is the free pin of the plumbing.
(a2) the FITTED gaussian model agrees with the Laplace fit to a
     tolerance that shrinks with N. Not machine precision.

## VALIDATION (b) IN THE PACKAGE: the probe design through frm()
`frm(bf(y ~ x + (x | g)) + binomial(), importance = N)`, seed 11,
60 groups of 8, correlated intercept and slope.

| fit | nll | intercept | slope | SE(int) | MCSE | min/med ESS | rounds |
|---|---|---|---|---|---|---|---|
| Laplace | 299.3017 | -0.6915 | 0.3592 | 0.1697 | - | - | - |
| importance 500 | 297.9598 | -0.7030 | 0.3611 | 0.1777 | 0.0557 | 0.930/0.978 | 3 (capped) |
| importance 2000 | 297.9147 | -0.7073 | 0.3648 | 0.1792 | 0.0308 | 0.931/0.970 | 3 |
| GLMMadaptive 15/25 | 297.9185 | -0.7047 | 0.3617 | - | - | - | - |

At N=2000 the nll misses AGHQ by 0.0038 against its own MCSE of
0.0308, i.e. 0.12 MCSE; the intercept misses by 0.0026 against an SE
of 0.1792, i.e. 0.015 SE. At N=500 the nll misses by 0.041 against
MCSE 0.0557, i.e. 0.74 MCSE. Laplace misses the nll by 1.38 and the
intercept by 0.013.
Cost: 14 s at N=500, 51 s at N=2000 (3 rounds, the report freeze and
sdreport included).

## STOPPING RULE (revised after measurement)
The round loop's fixed point is a theta at which the objective ANCHORED
at theta is stationary, so the re-anchored gradient is the criterion
that says so directly; the parameter move is the cruder backstop. The
loop stops on either. Measured need for this: at N=2000 the probe's
last move was 0.0013, just above the 1e-3 move tolerance, while the
re-anchored gradient was 0.0008, below grad_tol. Without the gradient
criterion that correct fit was flagged as capped.

## check_convergence and the Monte Carlo gradient
A corrected objective's gradient carries an O(N^-1/2) error that is
exactly zero only where the estimator is exact (a gaussian response at
its own anchor), so grad_tol cannot be met in general and judging an
importance fit by it warns on every correct fit. check_convergence
therefore judges an importance fit by ITS OWN criterion (the loop hit
its round cap with the estimates still moving) and leaves the gradient
visible in fit$importance$grad. Measured gradients at the reported
optimum: 0.0008 (probe, N=2000), 0.0269 (probe, N=500).

## Gaussian fitted agreement (validation a2), same design as a1
| N | nll - Laplace | max abs par diff | rounds | min ESS |
|---|---|---|---|---|
| 2 | +2.270e-01 | 5.7e-2 | 3 | 1.000000 |
| 10 | +1.395e-02 | 1.4e-2 | 3 | 1.000000 |
| 200 | +6.910e-03 | 9.6e-3 | 2 | 1.000000 |
| 1000 | +1.138e-07 | 4.0e-5 | 1 | 1.000000 |
Every weight is equal (ESS/N = 1.000000 exactly) at every draw count,
which is the (a) claim; the residual in the FITTED value is parameter
drift from the optimizer, and it vanishes by N=1000 where the loop
converges in one round.

## Determinism and RNG hygiene: both hold
Two fits at the same seed give identical objective and identical par
(`identical()`, not `all.equal`). A `set.seed(99); rnorm(3)` before and
after a fit gives identical draws, so the private stream neither reads
nor disturbs the session state.

## trunc() and cens(): SUPPORTED, measured not assumed
The Gauss-Kronrod refusal does not carry over. Its nodes map the whole
real line, so the truncation normalizer log(F(ub) - F(lb)) underflows
to exactly zero out there and the objective goes to -Inf; importance
draws sit within about four conditional standard deviations of the
mode, the same region the Laplace approximation itself lives in.

Left-truncated Poisson, 40 groups of 6, trunc(lb = 1):
  Laplace nll 448.8567 fixef 1.1815 0.3987
  importance 1000 nll 448.7131 fixef 1.1819 0.3986, min ESS 0.935,
  MCSE 0.0329, objective finite. `quadrature = TRUE` refuses the same
  model.
Right-censored lognormal, 40 groups of 6, 52 of 240 censored
(re-transcribed from the TEST's own construction; the original row came
from a scratch script lost in the wipe and quoted a different design):
  Laplace nll 311.1358
  importance 1000 nll 311.0847, min ESS 0.9753, median 1.0000,
  MCSE 0.0103.
weights() also works.

## Tape cost vs Laplace (probe design, n = 480, q = 2, G = 60)
| tape | RSS delta | fn | gr | build |
|---|---|---|---|---|
| Laplace | +2.3 MB | 1.5 ms | 2.5 ms | - |
| IS N=50 | +32.7 MB | 4.0 ms | 10.0 ms | 0.08 s |
| IS N=200 | +23.9 MB | 15.5 ms | 40.5 ms | 0.33 s |
| IS N=1000 | +177.3 MB | 73.5 ms | 284.5 ms | 1.95 s |
Times are means of 20 evaluations at perturbed parameters. They scale
linearly in N: the gradient costs about 0.28 ms per draw against a
2.5 ms Laplace gradient, so N=1000 is roughly 114x. The RSS column is
noisy at small N because the allocator reuses freed blocks; only the
N=1000 row should be quoted, and it says the tape is about 320 bytes
per (row x draw) as the standalone measurement did.

## tests/testthat/test-importance.R: 76 pass, 0 fail, 0 warn, 0 skip
Covers (a) through (h) plus trunc/cens, the fit record and the
per-group/joint pin. Two fixes it forced:
- ranef() keys by TERM LABEL ("1 | g"), not by grouping variable, so
  the test reads ranef(fit)[[1L]].
- importance_rounds DEFAULT RAISED 3 -> 5. Measured move sequences fall
  about tenfold per round and needed 2 rounds at N=2000, 3 at N=200 and
  4 at N=500; a cap of 3 produced spurious "still moving by 0.00105"
  warnings on correct fits. The loop exits on whichever criterion is
  met first, so an unused round costs nothing.

## A DEFECT FOUND BY THE VIGNETTE, and the guard it forced
Writing the vignette example turned up a design where the correction
returned a CONFIDENTLY WRONG answer: both variance components collapsed
to exactly 0.000, the log-likelihood came out WORSE than Laplace's, and
the reported diagnostics said min ESS 1.000 and MCSE 0.0000, because at
a collapsed variance component the conditional is a point mass and
every weight is trivially equal. The diagnostics were blind to it.

Root cause: the model is unidentified (4 to 6 binary rows per group
cannot identify a 2 x 2 covariance). The Laplace fit itself stops at
`singular convergence (7)` with the correlation parameter at -67.3. The
ESTIMATOR is fine there: the corrected objective at the Laplace
estimates is 152.829 against GLMMadaptive nAGQ=25's 152.851, better
than the Laplace value of 153.085. It is the fixed-point ITERATION that
diverges, wandering along the ridge to 156.056.

An anchor-ESS guard does NOT catch this: at a degenerate anchor the
proposal matches its own needle-thin conditional, so ESS at the anchor
is about 1. That guard was written and removed.

### The invariant that does catch it
The correction must not make the fit worse than the correction at the
Laplace estimates. Both values estimate the same exact marginal
likelihood, each at its own theta from a proposal frozen there, so they
compare directly. Measured rise (corrected at final minus corrected at
the Laplace estimates) across nine designs:

| design | rise | verdict |
|---|---|---|
| 40x4 seed11 (unidentified) | +0.352 | refused |
| 40x6 seed11 (unidentified) | +3.227 | refused |
| 50x6 seed3 | -0.009 | fine |
| 60x8 seed11 | -0.071 | fine |
| 40x5 seed7 wide | -0.028 | fine |
| 50x5 seed2 wide | -0.128 | fine |
| gaussian N=50 | +0.010 | fine |
| gaussian N=200 | +0.007 | fine |
| gaussian N=1000 | +0.000 | fine |

Tolerance `imp_worse_tol = max(3 * mcse_final, 3 * mcse_start, 0.1)`.
The 0.1 floor sits ten times above the worst legitimate rise (gaussian
optimizer drift, +0.010) and three and a half times below the smallest
divergence (+0.352); it is also a tenth of the log-likelihood a single
AIC unit is worth, so nothing that could change a model comparison
passes it.

## Vignette section: vignettes/diagnostics.Rmd ONLY
Per the coordinator, the section went in diagnostics.Rmd and
brms-migration.Rmd was NOT touched (confirmed by git status). New
section "When the Laplace approximation is the problem" with
subsections on reading the ESS, cost, and scope, placed after "Are the
standard errors trustworthy?" (whose list of remedies went from three
to four). Rendered output, 40 groups of 5 binary rows, seed 7,
importance = 500:
  laplace -119.5584, importance -118.4392
  "importance-corrected, 500 draws per group in 4 rounds
   (MCSE 0.11, min ESS 0.64 of 1)"
  ESS summary 0.6426 / 0.8973 median / 0.9924, MCSE 0.1085
Render time 21.5 s.

## tests/testthat/test-importance.R: 84 pass, 0 fail, 0 warn, 0 skip

## A REGRESSION THE SUITE CAUGHT (and the fix)
The `lp_eta_fixed()` refactor in R/objective.R moved the mo()/mi()/
offset arithmetic OUT of the objective closure. That closure sets
`"c" <- RTMB::ADoverload("c")` at its top; a separate function does
not inherit it, and the mo() simplex is built as `c(0, pars[[...]])`.
Base `c()` on (numeric, advector) silently DROPS the advector class,
so the tape ended in a bare complex vector and MakeADFun raised
"Invalid argument to 'advector' (lost class attribute?)".
Symptom: test-autoscale.R f=2 e=1 (the failing assertions were about
autoscale_plan() returning NULL, three steps downstream, because
fit_catching() had swallowed the real error).
Confirmed as MINE by running the same file against pristine HEAD
sources in a scratch copy: pristine 46 pass / 0 fail, mine 2 fail.
Fix: re-establish both overloads inside lp_eta_fixed(). After it the
autoscale objective is 239.5804661, bit-identical to pristine, and the
file passes 46/46.
Lesson worth keeping: factoring code out of an RTMB closure carries the
overloads with it, and the failure surfaces far from the cause.

## Refit paths: the correction is THREADED, not dropped
All four `fit_assembled()` refit sites already forwarded `quadrature`
by name (allfit.R:93, confint.R:1189, influence.R:88, sugar.R:239).
Without the same treatment, `confint(method = "profile")`,
`influence()`, `frm_allfit()` and `refit()` would silently rebuild an
importance fit as a LAPLACE fit and report intervals for a different
objective. `importance = <fit>$importance$draws %||% 0L` is now
forwarded at all four. Verified: refit() of a 100-draw fit comes back
with draws = 100 and the same objective to 1e-6.
`vcov_cluster()` REFUSES instead, matching its existing quadrature and
profile guards: the corrected objective is a sum over GROUPS of
reweighted integrals, not the per-row sum its cluster scores are read
off, so a cluster carrying no rows of its own still moves it.

## The several-blocks-one-factor follow-up: what it needs
Designed for, not built. What generalizes UNCHANGED:
- `build_importance_objective()`'s `zu` already loops over linear
  predictors and multiplies each `lp[["Z"]]` by the whole draw matrix,
  so two dpars each carrying `(1 | g)` need no change there.
- `imp_group_map()` already reads the row-to-level map off the
  sparsity of Z across ALL linear predictors, not off the formula.
- The Hessian block-diagonality gate is written against a per-position
  level vector, so it generalizes by concatenation.
- The precision extraction is per block already; it becomes a sum.
What has to change:
1. A group's coefficients stop being CONTIGUOUS. `b` is level-major
   within a block, but blocks are concatenated, so group k owns
   scattered positions. The `(k - 1) * q + seq_len(q)` spelling has to
   become an index vector per group, and it now appears in FIVE places,
   not the two the first draft of this note named: the scalar branch of
   imp_plan() (which assumes the layout implicitly), its general
   branch, the zsq assembly, the u_slices split in
   build_importance_objective(), and the group swap in imp_verify().
   The reviewer caught the zsq site; miss it and the half-norms would
   be read off the wrong rows with no error raised.
2. `imp_prior_terms()` becomes a sum over blocks, each contributing
   its own per-level Gaussian on its own slice of the group's
   coefficients. The blocks are independent given theta, so the
   per-group precision is block diagonal and no cross terms appear.
3. A check that the blocks share one grouping factor with identical
   levels, and that each row maps to the same level under every one of
   them. Blocks over DIFFERENT factors stay refused: that is the
   nested integral, which does not factorize.
None of this touches the estimator, the stacking, the reporting or the
diagnostics.

## THE API AS SHIPPED
`frm(formula, data, ..., quadrature = FALSE, importance = 0L, ...)`
- `importance` sits directly after `quadrature`, its sibling. Every one
  of the five `fit_assembled()` call sites passes by NAME after the
  first four positional arguments, so inserting a parameter there is
  safe (checked, not assumed).
- A positive integer is the draws per group; 0 is off. Validated with
  the package's own `check_count(importance, "importance", min = 0L)`,
  so a negative or fractional value is refused in the package's own
  words. An odd count is rounded UP by `imp_n_draw()` so the
  antithetic draws pair.
- `frm(importance =, quadrature = TRUE)` is refused.

`frmtmb_control(..., importance_seed = 1L, importance_rounds = 5L,
                importance_ess = 0.25, ...)`
- `importance_seed` feeds a PRIVATE RNG stream (save, set.seed, restore
  on exit), so the fit neither reads nor disturbs the session state.
- `importance_rounds` caps the refreeze loop (default raised from 3
  after measurement).
- `importance_ess` is the warning threshold as a fraction; > 1 is
  refused by name.

Recorded on the fit as `fit$importance`:
  draws, seed, rounds, moved, moves, grads, capped, grad,
  start_value, start_ess_min, ess_min, ess_median,
  ess (named by grouping level), mcse
`print()` and `summary()` both emit one line:
  "Marginal likelihood: importance-corrected, N draws per group in K
   rounds (MCSE m, min ESS e of 1)"
`logLik()`, `AIC()`, `BIC()`, `vcov()`, `confint()` all come from the
corrected objective through the ordinary sdreport path.
`ranef()`, `fitted()`, `predict()` work because the conditional modes
come from `solved_par_list(lap_obj, opt$par)`, the same recovery the
quadrature path already used.

### On the spelling
The brief said to implement `importance` and argue for a better one if
found. I did not find one. `importance = 2000` reads as "how many
importance draws", parallels `quadrature = TRUE` without colliding with
it, and leaves room for a future `importance = "auto"`. The
alternatives considered and rejected: `is =` (collides with the base
function and says nothing), `nIS =` (camel case, unlike every other
argument in the package), `draws =` (that word means posterior draws
everywhere else in this project, and `frmtmb.sample` owns it).

## DELIBERATE OMISSIONS
Refused by name, each with the remedy in the message:
- more than one random-effect block (nested integral)
- a block with no grouping factor: s(), t2(), gp(), hsgp()
- covariance structures that correlate the LEVELS: gr_cov, gr_prec,
  equalto, car, spde
- Student-t latent blocks (us_t, diag_t). They factorize over levels
  but are NOT gaussian in a level's coefficients, so the precision
  extraction in imp_prior_terms() would silently return a gaussian
  approximation to a t density. They fall out of the imp_covstructs
  whitelist and are refused with the reason.
- multi-membership: a row that reaches several levels, detected off
  the sparsity of Z rather than off the formula
- a second response, and rescor
- a family that supplies its own loglik (mixtures, HMM, latent class)
- residual correlation terms (ar, ma, arma, cosy, unstr)
- nonlinear predictors, cs(), mi()
- REML, frmtmb_control(profile = TRUE), quadrature
- dry_run = "objective" (no unfitted form: the proposal needs a mode)
- a map on b

Not refused, but not done:
- vcov_cluster() refuses a corrected fit rather than supporting it.
- rr blocks are outside the whitelist, so the b_idx/c_idx divergence
  that rr introduces is never exercised.
- autoscale and sparse_x are marked untested in the compat registry.
  Nothing suggests they break; nothing checks them either.
- No adaptive draw count. The user reads the MCSE and the effective
  sample sizes and decides.
- The parameter drift from optimizing a Monte Carlo surface is
  reported and bounded, not removed. It is O(N^-1/2) and was about a
  tenth of a standard error at 200 draws on the gaussian design.
- DESCRIPTION Version is untouched at 0.49.1 because DESCRIPTION is
  outside the lane's declared surface; NEWS.md opens a 0.50.0 section.
  Whoever integrates should bump it.

## FULL CORE SUITE (one file per process, load_all, NOT_CRAN=true)
102 / 102 files: failed=0 errors=0 warnings=0 skipped=4 passed=5347.
CORRECTED: the first run of this suite reported 43 skipped and 4903
passed, which was an artifact of running it BEFORE the Suggests were
installed into the private library. With all 60 Suggests present the
figure is 4 skips, and all four are opt-in tiers (three brms fit tests
behind FRMTMB_BRMS_FIT_TESTS, one Stan tier), not capability gaps.
Includes the new test-importance.R at 89 passed.
One regression was found and fixed during this run (the lp_eta_fixed
overload bug above); the run above is post-fix.

## PROBE RE-MEASURED WITH FINAL SETTINGS (round cap 5, gradient stop,
## divergence guard, re-anchored reporting)
| fit | nll | intercept | slope | MCSE | min/med ESS | rounds | grad |
|---|---|---|---|---|---|---|---|
| Laplace | 299.3017 | -0.6915 | 0.3592 | - | - | - | - |
| importance 500 | 297.9592 | -0.7033 | 0.3613 | 0.0559 | 0.930/0.978 | 5 | 0.0017 |
| importance 2000 | 297.9147 | -0.7073 | 0.3648 | 0.0308 | 0.931/0.970 | 3 | 0.0008 |
| GLMMadaptive 15 and 25 | 297.9185 | -0.7047 | 0.3617 | - | - | - | - |
N=2000: |297.9147 - 297.9185| = 0.0038 = 0.12 MCSE. Intercept off by
0.0026 against SE 0.1792 = 0.015 SE. Neither capped.
N=500 moves: 0.0725, 0.0301, 0.0088, 0.0023, 0.00061 (five rounds,
falling about threefold each, stopping on the move criterion).
N=2000 moves: 0.1494, 0.0191, 0.0013 (three rounds, stopping on the
gradient criterion at 0.0008 < grad_tol).
Wall clock: 41 s at 500 draws, 88 s at 2000, including sdreport.

## Per-round reporting (verbose = TRUE), as the brief asked
Every round prints its own min and median effective sample size beside
the move and the re-anchored gradient, so a proposal that stops
covering the integrand shows up round by round and not only at the end.
Sample trace, 60 groups of 4 binary rows, 200 draws:

  importance round [0.52s]: objective 147.11585, moved 0.127,
    re-anchored grad 0.396, ESS min 0.943 median 0.989
  importance round [0.45s]: objective 147.24297, moved 0.018,
    re-anchored grad 0.0291, ESS min 0.939 median 0.988
  importance round [0.50s]: objective 147.26185, moved 0.00133,
    re-anchored grad 0.00187, ESS min 0.938 median 0.988
  importance round [0.48s]: objective 147.26325, moved 8.49e-05,
    re-anchored grad 0.000119, ESS min 0.938 median 0.988
  importance report [0.00s]: objective 147.26334, min ESS 0.938,
    MCSE 0.071
  done [2.74s]: objective 147.26334, max|grad| 0.000119, 0 warnings

Note the objective RISES across rounds (147.116 -> 147.263). That is
the re-anchoring working as designed: each round's optimizer stops at
an optimistic value, and the fresh anchor prices it honestly.

## R CMD check --as-cran --no-manual  ->  Status: OK
_R_CHECK_CRAN_INCOMING_=false, all 60 Suggests installed so nothing was
skipped for a missing package. Zero NOTEs, zero WARNINGs, zero ERRORs.
  examples [22s] OK
  examples with --run-donttest [25s] OK
  tests [130s] OK
  re-building of vignette outputs [138s] OK
(Run on the tarball built before the per-round verbose ESS line; a
second build and check on the final source follows.)

## FINAL CHECK ON THE SHIPPED SOURCE -> Status: OK
Rebuilt and re-checked after the per-round verbose line.
  R CMD check --as-cran --no-manual, _R_CHECK_CRAN_INCOMING_=false
  0 NOTEs, 0 WARNINGs, 0 ERRORs. Status: OK
  examples --run-donttest [39s] OK, tests OK,
  re-building of vignette outputs [218s] OK
NOTE: R CMD check runs the suite WITHOUT NOT_CRAN, and
test-importance.R opens with skip_on_cran(), so the check never
exercises it. The per-file suite under NOT_CRAN=true is what tests the
correction; the check tests that everything else still builds and runs
around it.

## FINAL SURFACE (git status, nothing outside the brief)
 M NEWS.md
 M R/allfit.R R/compat.R R/confint.R R/fit.R R/influence.R
 M R/methods-fit.R R/objective.R R/sandwich.R R/sugar.R
 M man/frm.Rd man/frmtmb_control.Rd
 M vignettes/diagnostics.Rmd
 ?? R/importance.R
 ?? tests/testthat/test-importance.R
vignettes/brms-migration.Rmd is clean, as the coordinator asked.
dev/laplace-bias-probe.R was already untracked at handoff and was not
touched. No commits. Main checkout still at a233c3c.

# ===================================================================
# REVIEW ROUND: punch list P1 to P7 (reviewer verdict GO WITH FIXES)
# ===================================================================

## P3 (warm start) - APPLIED, and it did more than save time
`optimize_obj()` defaults `start_par` to `obj$par`, which MakeADFun
takes from the START TEMPLATE, so every round re-optimized from cold
and the Laplace warm start was used for the anchor only. Now
`start_par = par`.

Measured, cold against warm, same machine, scratch copy `is-cold` for
the cold build:

| design | cold | warm | rounds | result |
|---|---|---|---|---|
| probe, N=2000 | 46.6 s | 39.2 s | 3 both | nll 297.9147, fixef -0.7073 / 0.3648, min ESS 0.9312, IDENTICAL |
| gaussian, N=1000 | 5.41 s | 4.59 s | 1 both | diff from Laplace +1.14e-07 vs +8.08e-08 |

(gaussian is the median of four runs each; a single first run gave a
misleading 11.1 s, a cold-JIT outlier.)

About 15 percent faster on both, rounds unchanged, reported numbers
unchanged. The larger effect is on STABILITY: the 40x6 seed11 design,
whose cold-start iteration diverged to 156.056 and tripped the
divergence guard, now converges to 152.8382 against GLMMadaptive
nAGQ=25's 152.8506, a miss of 0.013. Its correlation parameter is
still unidentified and the round cap still warns, which is honest.
The guard is NOT dead: 40x4 seed11 still refuses.

Nine designs re-measured with the warm start:

| design | outcome |
|---|---|
| 40x4 seed11 | REFUSED (divergence guard) |
| 40x6 seed11 | nll 152.8382, 5 rounds, min ESS 0.9878 |
| 50x6 seed3 | nll 190.3768, 5 rounds, min ESS 0.9772 |
| 60x8 seed11 | nll 291.5411, 5 rounds, min ESS 0.8926 |
| 40x5 seed7 wide | nll 118.4392, 4 rounds, min ESS 0.6426 |
| 50x5 seed2 wide | nll 149.0300, 5 rounds, min ESS 0.7886 |
| gaussian N=50/200/1000 | 2 / 2 / 1 rounds, min ESS 1.0000 |

## P2 (profile intervals) - MEASURED, then guarded
The reviewer is right that R/confint.R:1190 is `anova_refit_ml()` and
unreachable for a corrected fit. `confint(method = "profile")` calls
`TMB::tmbprofile()` on the fit's own FROZEN tape.

Cost of the two options, measured on 40 groups of 4 binary rows, 500
draws:
- one importance fit: 4.3 s
- tmbprofile on the frozen tape: 20.2 s over 75 grid points
- one refreeze (plan + tape): 0.54 s
  so OPTION A (refreeze per grid point) is at least 75 x 0.54 = 41 s
  of taping ALONE per parameter, before any inner optimization: three
  times the profile it replaces at this size, and far worse at
  N = 2000 where a single tape is about 4 s.
- one ESS read on the frozen proposal: 0.128 s, and it reproduces the
  fit's reported anchor ESS exactly (0.951912 vs 0.951912).
  so OPTION B costs 0.26 s on a 20 s profile, 1.3 percent.

CHOSE OPTION B, and not only for the cost. Refreezing at each grid
point would make the profiled objective a DIFFERENT RANDOM FUNCTION at
every point, so the curve would acquire Monte Carlo noise between
adjacent points and the interpolation tmbprofile uses to find its
crossing would no longer be valid. Freezing the draws is what buys the
smooth deterministic curve in the first place; refreezing per point
throws that away to fix a smaller problem.

Implemented as `imp_frozen_proposal()` (rebuilds the fit's own
proposal deterministically from frame, estimates, draws and seed),
`imp_ess_at()` and `imp_profile_ess_warn()`, hooked into the profile
loop in R/confint.R. Verified end to end:
- `confint(fi, parm = "x", method = "profile")` silent.
- `confint(fi, parm = "theta_1", method = "profile")` warns naming the
  bound 0.704221 and the effective sample size 0.16 against the 0.25
  threshold.
- a Laplace fit profiles with none of the machinery.

## P1 (compat rows) - FIXED, and made self-policing
Added group `importance_refused_blocks` = us_t, diag_t, rr, gp, hsgp,
smooth with status "refused" and the reason. The three groups now
partition `covstruct_registry` exactly (13 + 5 + 6 = 24), and a test
asserts that partition plus the status of every key, so a structure
added later cannot quietly report "untested".

## P4 (the 0.43 ESS figure) - COULD NOT REPRODUCE EITHER NUMBER
Sweeping the design as described (40 groups of 3, sd 2.0, N=1000,
plan seed 1) over seeds 1, 2, 3, 7, 11, 101, 601:

| seed | min ESS at the Laplace optimum | min ESS at the corrected anchor (reported) |
|---|---|---|
| 1 | 0.6027 | 0.4907 |
| 2 | 0.5116 | 0.3545 |
| 3 | 0.1263 | 0.0807 |
| 7 | 0.2291 | 0.1598 |
| 11 | 0.7522 | 0.6035 |
| 101 | 0.7201 | 0.6247 |
| 601 | 0.4881 | 0.3222 |

My original 0.429 and the reviewer's 0.4802 are both inside this
spread; neither is "the" number. More important: this design is NOT a
witness that the floor clears every legitimate design, because at
seeds 3 and 7 the fit's own reported ESS is 0.08 and 0.16 and the
warning correctly FIRES.

So the justification is restated on what reproduces, and against the
anchor the fit reports from:
- ordinary designs sit far above: probe 0.93, scalar intercept 0.97,
  gaussian exactly 1.00, the vignette's design 0.64;
- a genuinely broken proposal sits far below: half a log-SD of
  displacement gives 0.009;
- and in between the threshold is MEANT to fire. It is not calibrated
  to clear every design anyone can build, and saying so is more useful
  than a margin that does not survive a seed change.

Restated in R/importance.R (imp_ess_floor), the vignette and NEWS.

## P5 (per-group pin) - STRENGTHENED, with its limit stated
`imp_verify()` now checks each group on its own, by comparing the
per-group log-joint DIFFERENCE between two draw columns: setting one
group to its column-j draw while the others stay at column k leaves
every other group's term unchanged, so it cancels and no per-group
absolute reference is needed. The total is still checked at both
columns, because a constant shared by every group cancels from the
differences. All groups up to 64, an evenly spaced 64 beyond.

Test: +5 on group 1 and -5 on group 2, in one draw column, leaves both
totals untouched and is caught per group.

STATED LIMIT, in the comment: an offset identical in every draw AND
summing to zero across groups cancels from differences and total
alike. The plain objective only ever returns sums, so no per-group
absolute reference can be recovered from it; that residue is the
price, and it is not a shape a plausible bug takes.

## P6 (findings accuracy) - CORRECTED in dev/importance-findings.md
- cens() row re-transcribed from the test's own construction: 52 of
  240 censored, Laplace 311.1358, corrected 311.0847, min ESS 0.9753,
  median 1.0000, MCSE 0.0103. The old row (40 of 240, 347.8633) came
  from a scratch script lost in the wipe. The compat note is corrected
  to match.
- suite figure corrected to 4 skips / 5347 passed, and the reason for
  the old 43 / 4903 recorded: the first run predated installing the
  Suggests into the private library.
- noted that R CMD check runs without NOT_CRAN and therefore never
  exercises test-importance.R at all.
- the several-blocks follow-up now names FIVE contiguity sites, not
  two: both branches of imp_plan(), the zsq assembly, the u_slices
  split, and the group swap in imp_verify().

## P6 (truncation) - STRENGTHENED rather than warned about
Restated in the compat note and the vignette: truncation does not
degrade the correction. Walking the bound past the response median so
that 64 percent of rows are discarded still leaves the worst group
holding 0.60 of its draws, more than twice the threshold, and a
gaussian right truncation at the median reports 0.97.

## P7 - APPLIED
- autoscale, sparse_x and mo() upgraded from "untested" to "works",
  each citing what was measured.
- `lp_eta_fixed()` now installs the same THREE overloads
  `frmtmb_ad_overload()` does, with a comment saying so and saying why
  it is spelled out rather than wrapped.
- one sentence added to ?frm: the reported covariance is a Hessian at
  a point stationary only to `fit$importance$grad`, and
  confint(profile) walks the frozen proposal.

## An accident during this round, and its repair
A `sed -i` whose line-number variable came back empty replaced every
line of R/fit.R with one comment line. Restored from the scratch copy
`is-cold/R/fit.R` (taken after P3 and before the roxygen edit),
verified at 1955 lines, +142/-11 against HEAD, parsing, with the
importance wiring intact; the doc edit was then reapplied with a
matched-text edit rather than a line number. R/importance.R was never
touched by the accident. Lesson: never `sed -i "${L}s/.*/.../"` without
asserting that L is non-empty.
