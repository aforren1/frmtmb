# Grammar fuzz findings

Run of the pairwise covering-array fuzzer (`tests/testthat/helper-fuzz.R`,
`tests/testthat/test-fuzz.R`, `dev/fuzz/run-fuzz.R`) on branch
`wt-registry`, at package version 0.25.0, after the invariant repairs
described under "Harness repairs" below.

| | |
|---|---|
| seed | `20260901` |
| specs | 300 (108 pair-cover rows, 181 3-way pad rows, 11 refusal probes) |
| feasible value pairs | 737, all covered |
| fits | 300 model fits, plus permutation and unit-weight refits |
| wall time | 134 s |
| command | `Rscript dev/fuzz/run-fuzz.R` |
| machine-readable | `dev/fuzz-findings.json` |

Classification counts: 15 REAL-NEW, 55 KNOWN-PENDING, 4
KNOWN-DIVERGENCE, 24 UNCONVERGED, 0 GENERATOR-BUG.

The previous run (branch `wt-quad`, version 0.23.0, same seed and plan
size) reported 16 REAL-NEW, 57 KNOWN-PENDING, 4 GRAMMAR-DIVERGENCE and
27 UNCONVERGED. The two runs are on different checkouts, so the
comparison is indicative rather than exact; the one REAL-NEW that
disappeared is R5, which the v0.26 backlog merges fixed.

`size` is now a hard cap on the plan, refusal probes included:
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
and `finite_loglik` are judged regardless. All 24 downgrades in this
run carry a real verdict from `check_convergence()`: 19 report false
convergence *and* a large gradient, 5 report false convergence alone.

---

## Harness repairs

The invariants themselves were audited before this run. Six of them
could not have failed, could not have been read, or muted findings they
had no business muting.

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

**Result: 34 of 35 intervals covered 0.4** (the 1e-4 binomial floor at
n = 35 is 27), and no row of any interval departed from the identity.

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

---

## REAL-NEW

### R1. `quadrature = TRUE` still breaks down on hard likelihoods

7 findings (`fit_error`), all the same informative refusal naming
`quadrature` rather than a bare RTMB string:

```
family=Beta        aterm=weights  re=ri     dpar=dpar_x  seed 20270671
family=weibull     aterm=cens     re=nested              seed 20276533
family=binomial    aterm=trials   re=nested              seed 20312682
family=cumulative  aterm=weights  re=nested              seed 20336130
family=lognormal   aterm=cens_chr re=nested dpar=dpar_x  seed 20456301
family=Beta        aterm=weights  re=ri_fx               seed 20492450
family=weibull     aterm=cens_chr re=ri     dpar=dpar_x  seed 20520783
```

`quad_fit()` tries three calibration points - the Laplace optimum, the
best point a broken tape reached, and the cold template - and reports
the failure only when none of them yields a finite objective and
gradient. Two shapes dominate:

- **A singular variance component.** The `cumulative` case above has
  `sd = 8.8e-05` on its inner nested block at the Laplace optimum.
  `logIntegrate_t::rescale_integrand` measures curvature by finite
  differences with `dx = 1`, which says nothing useful about an
  integrand 1e-4 wide, so `sigma = 1/sqrt(-h(mu))` comes back
  meaningless.
- **Nested blocks on a small fuzz dataset.** `(1 | ga/gb)` is an
  iterated integral, and the outer integrand is itself the output of a
  frozen inner rescaling. Clean data fits (the regression tests cover
  poisson, Gamma and Beta over `(1 | ga/gb)`); 15 inner levels over 90
  observations does not always.

Both would want an integrator that recalibrates per evaluation.
TMBad's `adaptive = TRUE` is meant to be that and is measurably worse:
it fails on negbinomial cases the frozen path fits.

### R2. `vcov()` returns non-finite entries with no warning

2 findings (`vcov_psd`), and the point is that neither fit warned about
anything:

```
family=negbinomial aterm=none re=nested special=mo_int mode=sparse_x
  seed 20323429   vcov() has non-finite entries (36 of them)
family=lognormal aterm=cens_chr re=dbar_f special=smooth dpar=dpar_x
  mode=profile  seed 20389865   vcov() has non-finite entries (25)
```

The `profile` case takes the `solve_joint_precision()` branch of
`vcov.frmtmb_fit`, which now degrades to `NaN` standard errors plus one
warning naming `diagnose()` and `vignette('diagnostics')` instead of
letting `solve()` throw. The ML branch reads an already-inverted
`cov.fixed` and returns `NaN` entries with no warning at all. The two
paths now agree on the value and still disagree on whether the user is
told. Six further `vcov_psd` findings in this run are on fits that did
warn about convergence and are filed as UNCONVERGED.

### R4. Bare `NA/NaN gradient evaluation` outside quadrature

2 findings, `ar1` or a wide smooth crossed with a special:

```
family=Gamma      aterm=weights re=ar1 special=mo_int mode=profile   seed 20299981
family=cumulative aterm=none    re=rs  special=smooth mode=autoscale seed 20530553
```

`ar1(tt + 0 | g) + s(xs)` on clean data under ML fits without
complaint, so these look data-conditional rather than structural. They
are listed because the error surface is uninformative, not because the
combination is known to be unsupported.

### N5. `logLik` moves about 1e-6 under a row permutation (marginal)

4 findings (`row_permutation`), three under `profile` and one under
`sparse_x`:

```
family=student     aterm=weights re=diagbar special=smooth mode=profile  seed 20381072  diff 1.0e-06
family=negbinomial aterm=weights re=ar1     special=smooth mode=sparse_x seed 20433830  diff 2.2e-06
family=student     aterm=se      re=ar1     special=mo     mode=profile  seed 20501243  diff 2.4e-07
family=cumulative  aterm=weights re=ar1     special=smooth mode=profile  seed 20524691  diff 1.3e-06
```

The designs are elementwise identical after undoing the permutation, so
this is not a design-assembly problem; it is an optimizer landing a
hair away from the same optimum, and the invariant's thresholds (1e-6
on `logLik`, 1e-8 on the objective at a shared parameter vector) are
tighter than `nlminb` converges on a flat surface. Investigated in the
previous run on clean `negbinomial + weights() + ar1()` data and closed
there: the standardization constants are bit-identical under the
permutation, dropping `autoscale` makes the gap larger rather than
smaller, and tightening `rel.tol`/`x.tol` to 1e-12 does not move it.
Recorded because the differences clear the stated thresholds, not
because they look like defects.

### R5 (gone). A model with no free parameters dies inside nlminb

`family=poisson re=none mode=profile`, seed 20477795, used to fail with
`'d' must be a nonempty numeric (double) vector` - `profile = TRUE`
moves every coefficient inside, leaving nothing for the outer
optimizer. The spec fits on this checkout, so the v0.26 backlog merges
closed it.

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
rediscovered here, and none of the 15 REAL-NEW findings is one of them.
`rescor x cens`, non-estimable rank-deficient prediction, `ar1` with
gapped integer levels, and `anova` across different `nobs` are outside
the grammar for the same kind of reason (single-response specs,
full-rank designs, balanced factor time grids, no model comparison).

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
failing on it; any divergence not on that list now fails the tier
instead of being filed away unread.

---

## Generator notes

Four generator defects were found and fixed during the first triage,
and the six harness repairs above during the second. This run has no
GENERATOR-BUG findings.

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

Two structural lessons worth keeping:

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
