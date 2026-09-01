# Grammar fuzz findings

Run of the pairwise covering-array fuzzer (`tests/testthat/helper-fuzz.R`,
`tests/testthat/test-fuzz.R`, `dev/fuzz/run-fuzz.R`) on branch
`wt-fuzzres`, at package version 0.27.0, after the fixes described
under "Closing the REAL-NEW findings" below.

| | |
|---|---|
| seed | `20260901` |
| specs | 300 (108 pair-cover rows, 181 3-way pad rows, 11 refusal probes) |
| feasible value pairs | 737, all covered |
| fits | 300 model fits, plus permutation and unit-weight refits |
| wall time | 159 s |
| command | `Rscript dev/fuzz/run-fuzz.R` |
| machine-readable | `dev/fuzz-findings.json` |

**The tier is green.** Classification counts: 0 REAL-NEW, 0
GENERATOR-BUG, 55 KNOWN-PENDING, 6 KNOWN-REFUSAL, 4 KNOWN-DIVERGENCE,
16 UNCONVERGED. Every finding is on one of the three lists in
`helper-fuzz.R` with the reason it stands, which is what
`test-fuzz.R` asserts.

The previous run (branch `wt-registry`, version 0.25.0, same seed and
plan size) reported 15 REAL-NEW, 55 KNOWN-PENDING, 4 KNOWN-DIVERGENCE
and 24 UNCONVERGED. On this checkout, before the fixes below, the same
plan reported 14 REAL-NEW (9 `fit_error`, 4 `row_permutation`, 1
`vcov_psd`) and 22 UNCONVERGED: one of the two R2 findings had already
been closed by the v0.26/0.27 merges.

Stability was checked twice: the full plan and `FRMTMB_FUZZ_N = 60`
(60 specs, 23 s, 10 KNOWN-PENDING, 2 KNOWN-REFUSAL, 2 UNCONVERGED, 0
REAL-NEW), through both `dev/fuzz/run-fuzz.R` and the test tier.

`size` is a hard cap on the plan, refusal probes included:
`fuzz_plan(size = 20)` returns 20 specs (9 cover rows plus the 11
probes), not 119. A short plan keeps the prefix of the greedy cover,
which is the part that buys the most pairs, and
`attr(plan, "n_pairs_covered")` reports what that prefix reaches.

UNCONVERGED is not a defect class. The estimate-dependent invariants
(`row_permutation`, `vcov_psd`, `confint_wald`, `simulate_mean`,
`simulate_mean_rows`, `unit_weights`, `loglik_identity`) only say
something about a fit that reached an optimum; when the fit itself
warned about convergence, a violation is downgraded so it does not
compete with real defects. Assembly identities
(`predict_eq_fitted`, `vcov_dim`), design agreement, crashes, refusals
and `finite_loglik` are judged regardless. All 16 downgrades in this
run carry a real verdict from `check_convergence()`: false convergence,
a maximum gradient far above `grad_tol`, or both.

---

## Closing the REAL-NEW findings

### R1. `quadrature = TRUE` breaks down on hard likelihoods

**7 findings; 1 fixed, 6 are now documented refusals.**

The recorded explanation - a singular variance component, or nested
blocks - was only half right, so the seven were re-measured one by one
against a ladder of calibration templates (the Laplace anchor, the cold
template, the anchor with the conditional modes zeroed, and the anchor
with `theta` displaced by -2 to +2 in steps of 0.5, each with three
recalibration rounds). Two clean classes came out of it, and only one
of the seven Laplace optima has a singular component at all.

**Nested blocks (4 findings, not recoverable).** `obj$fn(obj$par)` is
`NaN` at *every* calibration point tried, before the optimizer takes a
step. `(1 | ga/gb)` is an iterated integral: the outer integrand is
itself the output of the inner rescaling, and no calibration of the
outer one repairs that.

```
family=weibull    aterm=cens     re=nested             seed 20276533
family=binomial   aterm=trials   re=nested             seed 20312682
family=cumulative aterm=weights  re=nested             seed 20336130
family=lognormal  aterm=cens_chr re=nested dpar=dpar_x seed 20456301
```

The `cumulative` case is the one with a singular component as well
(`sd = 9.9e-05` on the outer block); the other three have healthy
components (0.11 to 0.47) and fail purely on the iterated integral.

**A single scalar block (3 findings, 1 fixed).** Here the tape is fine
- `obj$fn` is finite at every calibration point - and the *optimizer*
dies as the frozen rescaling expires under it. Displacing `theta` is
what changes the widths the transform froze, and for one of the three
it is enough:

```
family=weibull  aterm=cens_chr re=ri    dpar=dpar_x seed 20520783  FIXED
family=Beta     aterm=weights  re=ri    dpar=dpar_x seed 20270671  refused
family=Beta     aterm=weights  re=ri_fx             seed 20492450  refused
```

`quad_fit()` now tries, in order: the Laplace optimum, the best point
the optimizer reached (as before), the cold template, and the Laplace
optimum with `theta` displaced by -0.5, +0.5 and -1. Seed 20520783
reaches a stationary point at the first displacement (objective
152.6868, max|grad| 5.0e-04 against a `grad_tol` of 1e-3; the anchor's
own tape reads 152.6895 there). The two `Beta` cases resist the whole
sweep. Every candidate costs a tape, and none of them runs unless the
one before it failed, so a healthy quadrature fit pays nothing.

The refusal message now names the configuration instead of saying
"this likelihood":

```
quadrature = TRUE could not marginalize this model: the Gauss-Kronrod
objective broke down (a non-finite objective or gradient) at every
calibration point - the Laplace optimum, the best point the optimizer
reached, the cold start, and the optimum with the integrand widths
displaced. Here the model asks for an iterated integral over 2 nested
blocks ('1 | gb:ga' (30 levels), '1 | ga' (15 levels)) on 90
observations, and the outer integrand is itself the frozen rescaling of
the inner one. Refit with quadrature = FALSE: the Laplace approximation
fits this model (it is what the warm start above already did).
```

The last clause is a fact, not a hope: `quad_fit()` fits the plain
Laplace objective first as its warm start, so by the time this error is
raised the Laplace fit has already succeeded. The refusal is listed in
`FUZZ_KNOWN_REFUSAL` with that reasoning. An integrator that
recalibrates per evaluation would be the real fix; TMBad's
`adaptive = TRUE` is meant to be that and is measurably worse (it fails
on negbinomial cases the frozen path fits).

### R2. `vcov()` returned non-finite entries with no warning

**2 findings, both closed.**

```
family=negbinomial aterm=none re=nested special=mo_int mode=sparse_x
  seed 20323429   vcov() has non-finite entries (36 of them)
family=lognormal aterm=cens_chr re=dbar_f special=smooth dpar=dpar_x
  mode=profile  seed 20389865   vcov() has non-finite entries (25)
```

The `profile` case no longer reproduces on this checkout: it converges
with `max|grad| = 6.8e-07`, `pdHess = TRUE`, and a finite covariance.
The v0.26/0.27 merges closed it.

The `sparse_x` case reproduces, and the mechanism is not the one the
previous run guessed. The fit converges (code 0, `max|grad| = 9.5e-05`
against a `grad_tol` of 1e-3), `sdreport` reports **`pdHess = TRUE`**,
and all 100 entries of `cov.fixed` are `NaN`. Recomputing the outer
Hessian by hand shows why: it is finite, `chol()` succeeds on it - which
is exactly what `pdHess` reports - and its eigenvalues run down to
5.4e-14 with a reciprocal condition number of 5.8e-17, so `solve()`
refuses it as computationally singular and TMB fills `cov.fixed` with
`NaN` while leaving `pdHess` alone. The `mo()` simplex parameters are on
the boundary (`zeta1 = -17.9, -29.9`), which is what flattens it.

`diagnose()` already had the verdict (`bad_se` names all ten
parameters). Nothing else did, and `check_convergence()` could not: it
only looked at `pdHess`, and the two verdicts do not coincide.

Three changes, so the value can never be returned mute:

- `vcov.frmtmb_fit()`'s ML branch warns when the covariance it is about
  to return holds non-finite entries, pointing at `diagnose()` and
  `vignette('diagnostics')` (`warn_nonfinite_cov()` in `R/utils.R`).
- `solve_joint_precision()` raises the same warning when the joint
  precision inverts *successfully* into non-finite entries, not only
  when `Matrix::solve()` throws. The `profile` branch was as quiet as
  the ML one in that case.
- `check_convergence()` reports it at fit time as a verdict of its own
  when `se = TRUE` has already run `sdreport`, since `pdHess` alone
  does not imply usable standard errors.

The invariant was reading the wrong warnings: it judged a `vcov()`
result against the *fit's* warnings only, so a warning issued by
`vcov()` itself could not clear it. It now reads the warnings of the
call it already makes (harness repair H8 below). A NaN standard error on
a rank-deficient fit is a legitimate degradation; a silent one is not,
and that is what the invariant now says.

### R4. Bare `NA/NaN gradient evaluation` outside quadrature

**2 findings, both fixed; and no bare optimizer error escapes `frm()`
any more.**

```
family=Gamma      aterm=weights re=ar1 special=mo_int mode=profile   seed 20299981
family=cumulative aterm=none    re=rs  special=smooth mode=autoscale seed 20530553
```

Two different failures that looked identical from outside.

**The path crosses a hole (seed 20299981).** The start is fine and the
optimizer walks into a region where the profiled objective is undefined
- a `Gamma` shape intercept on its way to `exp(34)`. Restarting from
`obj$env$last.par.best` does not help here: that point evaluates to
`NaN` too. The plain Laplace objective over the same model - the one
every other estimation mode optimizes - has no such barrier and
optimizes to `logLik = 20.03`, so its optimum is a starting point on the
far side. Under `profile = TRUE` that is now the fallback, which is the
same recipe `quad_fit()` uses to calibrate the Gauss-Kronrod tape at the
Laplace optimum. The fit then completes and warns loudly (false
convergence (8), `max|grad| = 60`), so the tier files it as
UNCONVERGED rather than trusting it.

**The start itself is unusable (seed 20530553).** `fn` at the starting
values is `Inf` and one gradient entry is `NaN`, so `nlminb` dies before
its first step. The cause is not autoscale's column scaling - the plan
touches one column, `x`, with a factor of 1.083 - but `(1 + x | g)`
collapsing to a perfect correlation, which drives `theta` to +-300 (the
plain ML fit of the same data lands at 302 / -346 and reports "singular
convergence (7)"). The autoscale pre-fit finds that runaway point and
hands it over as the warm start. The cold start `make_start()` would
have built fits the model, so it is now offered as a fallback whenever
the fit did not begin there. Result: `logLik = -117.86`, matching the
plain ML fit, with the singular-convergence warning.

Both fallbacks live in `fit_recovery_starts()` and run only after a
failure, so a healthy fit pays nothing for them. `optimize_obj()`
separately retries once from `obj$env$last.par.best` when the optimizer
aborts, which is the cheap case of the same idea.

**The wrapper.** `nl_start_context()` became `fit_error_context()` and
now covers every optimizer error, not just the nonlinear-zero-start one:

```
The optimizer failed on this model (Gamma, ML, profile, nlminb): NA/NaN
gradient evaluation. The likelihood was undefined or unbounded somewhere
the optimizer stepped. Refit with verbose = TRUE to see which stage
broke, try another optimizer (frmtmb_control(optimizer = "optim")), or
start the fit nearer the optimum with the `start` argument of frm().
```

The condition carries a `frmtmb_fit_error` class so the autoscale
pre-fit - an inner `fit_assembled()` that has already been through the
wrapper - is rethrown rather than wrapped twice. A nonlinear model
without `start` still gets the sharper `start=` advice it had.

### N5. `logLik` moved about 1e-6 under a row permutation

**4 findings, all closed - three by giving the invariant the right
tolerance, one by moving where it probes.**

```
family=student     aterm=weights re=diagbar special=smooth mode=profile  seed 20381072
family=negbinomial aterm=weights re=ar1     special=smooth mode=sparse_x seed 20433830
family=student     aterm=se      re=ar1     special=mo     mode=profile  seed 20501243
family=cumulative  aterm=weights re=ar1     special=smooth mode=profile  seed 20524691
```

Tightening the generator was tried first and does not work.
Re-running all four with `rel.tol = x.tol = 1e-13`, `grad_tol = 1e-6`
and `restarts = 5` reproduces the same differences to twelve digits.
The optimizer is not what limits them.

What limits them is the shape of the likelihood. The two fits land on
genuinely different points of a flat ridge - `||dp||_1` between 0.89 and
2.44 - each with a gradient at its stopping tolerance. To first order
the objective can differ across them by the gradient times the step, and
Hoelder pairs the sup-norm gradient the optimizer judges itself by with
the 1-norm step: `|dlogLik| <= max|g| * ||dp||_1`. Measured:

| seed | max abs g1 | max abs g2 | dp 1-norm | d logLik | bound |
|---|---|---|---|---|---|
| 20381072 | 9.3e-06 | 1.19e-05 | 0.895 | 1.02e-06 | 1.90e-05 |
| 20433830 | 1.33e-05 | 4.33e-06 | 1.193 | 2.21e-06 | 2.10e-05 |
| 20501243 | 2.46e-04 | 1.07e-04 | 2.442 | 1.01e-04 | 8.63e-04 |
| 20524691 | 6.89e-08 | 1.76e-05 | 1.778 | 1.31e-06 | 3.14e-05 |

Every one clears its bound by about an order of magnitude.
`fuzz_permutation_tol()` computes `(g1 + g2) * ||dp||_1` in the fit's own
natural units (`par_units`, so the product is unit-free) and keeps 1e-6
as a floor. This is not a blanket bump: on a well-conditioned fit the
bound is far below the floor and the threshold is exactly what it was.
It rises only when the two optima are genuinely far apart, which is the
one circumstance under which a gap is owed.

The fourth finding (seed 20501243) tripped the *other* branch, the
shared-parameter objective comparison, at `2.4455e-07` against a 1e-8
threshold. That one is a harness defect, not a tolerance question - see
H7.

---

## Harness repairs

Eight now. H1 to H6 were made before the previous run and are kept here
as history; H7 and H8 are new.

### H1. `fuzz_inv_confint` could not fail

Both of its assertions reduced to `se < 0`. "Is the lower bound below
the upper bound" and "is the estimate inside its own interval" are the
same statement about `est +/- z * se`, and a Wald interval satisfies
both by construction whatever the covariance is.

It now checks the arithmetic against its own input: `confint(method =
"wald")` must equal `est +/- qnorm(0.975) * sqrt(diag(vcov(full =
TRUE)))` row by row, which catches a misaligned or misscaled row rather
than restating the formula. Under REML and `profile = TRUE` the
coefficients are inner parameters and `vcov()` cannot report the outer
block, so the identity is skipped there and the ordering check stands
alone.

The invariant also accumulates coverage of the one coefficient the
generator knows the truth of (`FUZZ_BETA_X = 0.4`, the slope every spec
is drawn with) and the plan is judged as a whole, since a single
interval says nothing. Specs whose mean structure the generator
distorts are left out: `cens()` clamps the response at its own
quantiles, `se()` and `weights()` change what the reported standard
error means, and an ordinal response goes through a latent scale whose
sign convention is the family's. Both the `confint` and the `vcov` op
contribute, because the Wald check reads the sdreport covariance
`vcov()` has already built and is nearly free there.

**Result: 35 of 36 intervals covered 0.4** (the 1e-4 binomial floor at
n = 36 is 28), and no row of any interval departed from the identity.

### H2. brms grammar divergences were recorded and never asserted

`fuzz_inv_brms()` filed them under a class of their own, which
`fuzz_triage()` passed through untouched and `test-fuzz.R` never looked
at. A divergence is now a finding like any other unless it is listed in
`FUZZ_KNOWN_DIVERGENCE` with the reason it is acceptable. The list is
seeded with the one divergence the migration vignette documents
(binomial without `trials()`), which is what the four
KNOWN-DIVERGENCE records below are.

### H3. The unconverged demotion fired on line-search noise, and missed the real verdict

Two defects, in opposite directions.

`FUZZ_NONCONVERGENCE` matched `"NA/NaN function evaluation"`, which
`nlminb` prints whenever a line-search step lands outside the support.
Healthy fits do that routinely, so a violation on a perfectly converged
fit was being filed as "unconverged". It now matches only the three
verdicts `check_convergence()` (`R/fit.R`) issues about the optimum
itself: false convergence, a large maximum absolute gradient, and a
non-positive-definite Hessian.

That fix alone made things worse, and revealed the second defect. The
runner recorded `head(res$warnings, 3)`, and `check_convergence()`
warns *last*, after however much chatter the optimizer produced - so on
exactly the fits that had trouble, the verdict was truncated away and
the finding looked clean. A run with the pattern fixed but the
truncation in place reported 29 REAL-NEW, of which 14 were fits that
had in fact reported `false convergence (8)` with a maximum gradient
around 1e10. The convergence warnings are now picked out before the
list is trimmed.

`predict_eq_fitted` and `vcov_dim` were also removed from the
convergence-sensitive set. They are assembly identities - they say two
code paths read the same design and the same bookkeeping - and hold at
any parameter vector at all, so a bad optimum is no excuse for breaking
them.

### H4. `simulate_mean` collapsed to two scalars

It compared `mean(colMeans(draws))` with `mean(predict())`. Two scalars
agree whenever the two vectors are permutations of each other, so the
check was blind to exactly the row-order class of defect the rest of
the harness hunts.

`simulate()` conditions on the fitted modes (`re.form = NULL`), so row
i's draw mean is `predict(type = "response")[i]` plus noise of a known
size. The new check standardizes each row by the draws' own standard
error and compares the sum of squares against a one-in-a-billion
chi-square tail. **Result: no findings**, over every `simulate` spec in
the plan; the threshold clears Monte Carlo noise at a 200-draw budget
with a wide margin.

### H5. Two refusal-triage rules matched on the spec alone

`bar-crossed-star` and `mo-factor-interaction` claimed every finding on
their refusal probe. A refusal probe can fail three ways - the
combination was accepted, the refusal came out as an internal crash, or
the message does not name the reason - so a genuine crash on those
specs would have been triaged away as already known. Each rule now
names the one symptom its pending fix addresses. Neither probe produced
a finding in this run.

### H6. The brms oracle judged `sparse_x` specs against a dense design

The claim in the file was that estimation modes never reach the data
layer. That is true of REML, profile, quadrature and autoscale, and
false of `sparse_x`: `frm()` passes `control$sparse_x` into
`assemble_frame()` *before* `dry_run = "frame"` returns (`R/fit.R:94`).
The oracle was calling `frm(..., dry_run = "frame")` without the
spec's control, so it compared brms against a dense design the spec
never fits. It now builds the frame under the spec's own control; the
population-level design comes back as a `dgCMatrix` with the same
values, and the comparison is clean.

### H7. The shared-parameter objective was compared at the one point where it cannot be

`fuzz_inv_permutation()` evaluated both objectives at `fit$opt$par` -
the *first* fit's own optimum. `obj$fn` re-solves the inner problem
starting from wherever the tape was last left, and a tape that has just
finished optimizing is left AT its solution, so its inner Newton takes
no step, while the other tape arrives from its own optimum and stops one
inner tolerance away. The invariant was measuring that asymmetry.

The evidence is unambiguous. For seed 20501243 the two objectives differ
by `2.4455e-07` at `p1`, and by **exactly zero** at five random
perturbations of `p1`, at `p2`, and at the midpoint of the two optima.
Warming the tapes through other points first does not move the number;
only the choice of evaluation point does. The two inner solutions at
`p1` differ by 1.4e-08 elementwise.

The probe point is now the midpoint of the two optima - a point neither
optimizer stopped at, so both tapes re-solve the inner problem the same
way, and one that probes the same claim: that the two objectives are the
same function of the parameters. Over eleven control specs drawn from
the plan the midpoint difference is at most 1.7e-13, so the 1e-8
threshold is if anything now tighter than it was.

### H8. `vcov_psd` judged `vcov()` by the fit's warnings, not its own

A NaN standard error is a legitimate degradation on a rank-deficient
fit; a silent one is not, and that distinction is what the invariant is
for. But it only had `sp$.warned`, collected when the fit was built,
so a warning raised by `vcov()` itself could never clear it - which
made the R2 fix unobservable to the tier that found R2. `fuzz_try()`
already collects the warnings of the call it wraps; the invariant now
reads them, and files a finding only when the non-finite result came
with nothing said.

---

## KNOWN-PENDING (rediscovered)

55 findings across 3 of the in-flight items.

| item | findings | how it surfaced |
|---|---|---|
| `(1 \| factor(x))` brms translation | 34 | `brms_translation`: brms rejects `(1 \| factor(gc))`, frmtmb fits it |
| `(f \|\| g)` builds a correlated block | 17 | `brms_re_blocks`: 4 random coefficients per level where brms has 3 |
| truncation ignored by the post-fit surface | 4 | `trunc_support` / `simulate_mean`: `simulate()` draws below the bound |

Note on the truncation item: comparing `simulate()` against
`predict(type = "response")` does **not** see it, because both ignore
the bound in the same direction and therefore agree - and that is true
of the per-row comparison added in H4 as well, for the same reason. The
detector that works is the support check: a simulator must not draw
outside the support the model was fitted under. That is
`fuzz_inv_trunc_support()`.

The sibling defect list in flight alongside this run - `cens()` crossed
with `trunc()` in the likelihood, `quadrature` with `se.fit`, the
`vreal()` newdata drop, and the `frm_sample()` jitter bounds - is
outside the generated grammar. No spec carries two addition terms on
one response, no op requests `se.fit`, `vreal()` is not a generated
value, and `frm_sample()` is not an op. None of them can be
rediscovered here. `rescor x cens`, non-estimable rank-deficient
prediction, `ar1` with gapped integer levels, and `anova` across
different `nobs` are outside the grammar for the same kind of reason
(single-response specs, full-rank designs, balanced factor time grids,
no model comparison).

---

## KNOWN-REFUSAL

6 findings, all one refusal: **`quadrature = TRUE` on a likelihood
whose Gauss-Kronrod calibration cannot be made to hold** (R1 above).

The rule in `FUZZ_KNOWN_REFUSAL` claims only the `fit_error` invariant,
only on a `mode=quadrature` spec, and only on the message
`quad_breakdown_message()` produces. A crash, a silent acceptance, or
any other error on the same spec is still a finding.

The refusal is safe rather than merely tidy: `quad_fit()` fits the plain
Laplace objective first as its warm start, so a model reaching this
error has already been fitted successfully under `quadrature = FALSE`,
which is what the message tells the user to do.

---

## KNOWN-DIVERGENCE

4 findings, all one divergence: **frmtmb accepts `binomial()` without
`trials()`, brms does not.**

```
brms::bf(y ~ 1 + x + (1 + x | g))  with binomial()
  -> "Specifying 'trials' is required for this model."
frm(bf(y ~ 1 + x + (1 + x | g)) + binomial(), data = d)
  -> fits, with trials defaulting to 1
```

frmtmb follows `stats::glm`, where a 0/1 response under `binomial()`
is Bernoulli; brms requires `trials()` or `bernoulli()`. The divergence
is in the permissive direction, so brms code ports to frmtmb
unchanged - the reverse does not.

**Documented.** `vignettes/brms-migration.Rmd`, "What changes". It is
listed in `FUZZ_KNOWN_DIVERGENCE`, which is what keeps the tier from
failing on it; any divergence not on that list fails the tier instead
of being filed away unread.

---

## Generator notes

Four generator defects were found and fixed during the first triage,
and the harness repairs above during the second and third. This run has
no GENERATOR-BUG findings.

1. The unit-weight metamorphic refit reused the `weights(w)` formula
   (random weights) instead of `weights(one)`, so it "failed"
   everywhere.
2. The brms smooth and gp comparisons sliced the first columns of
   `lp$Z`, which spans every random column of the model, instead of the
   special's own block columns. This produced a spurious HSGP
   mismatch.
3. The brms acceptance comparison used fit success, so a quadrature
   estimation failure was reported as a grammar divergence. It now
   compares `dry_run = "frame"`.
4. `length(unique(y)) < 3` rejected every Bernoulli dataset, and the
   `vcov()` dimension reference counted dpars held at a constant (which
   `se()` creates and `vcov()` documents that it excludes).

Structural lessons worth keeping:

- Comparing two fits' objectives at a shared parameter vector is only
  valid when the two designs are elementwise identical. mgcv rebuilds a
  smooth basis in a different rotation on reordered data - the same
  model, a different meaning for `beta` - so `fuzz_inv_permutation()`
  checks the design first and only compares objectives when the
  parameterization is shared.
- An invariant that summarizes before it compares can be unfalsifiable
  (H1) or blind to a whole class of defect (H4), and a triage rule that
  is too generous (H3, H5) hides findings just as effectively as a
  missing invariant. Both failure modes look exactly like a green run.
- Where an invariant evaluates matters as much as what it compares. H7
  is a check that was exactly right about a quantity measured at
  exactly the wrong point - the one point in the parameter space where
  a re-solved inner problem is guaranteed to be asymmetric between the
  two fits.
- An invariant must read the warnings of the call it makes, not only
  the ones the fit made (H8). "Non-finite, and nothing said" is a
  different claim from "non-finite", and only the first is a defect.
