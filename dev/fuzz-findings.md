# Grammar fuzz findings

Run of the pairwise covering-array fuzzer (`tests/testthat/helper-fuzz.R`,
`tests/testthat/test-fuzz.R`, `dev/fuzz/run-fuzz.R`) on branch `wt-quad`
after the quadrature defect wave, at package version 0.23.0.

| | |
|---|---|
| seed | `20260901` |
| specs | 311 (108 pair-cover rows, 192 3-way pad rows, 11 refusal probes) |
| feasible value pairs | 737, all covered |
| fits | 311 model fits, plus permutation and unit-weight refits |
| wall time | 120 s |
| command | `Rscript dev/fuzz/run-fuzz.R` |
| machine-readable | `dev/fuzz-findings.json` |

Classification counts: 16 REAL-NEW, 57 KNOWN-PENDING, 4
GRAMMAR-DIVERGENCE, 27 UNCONVERGED, 0 GENERATOR-BUG.

Against the same plan on the pre-fix checkout the counts were 28
REAL-NEW and 31 UNCONVERGED, so the quadrature work below cleared 12
REAL-NEW findings. (The plan itself changed by one refusal probe:
`quadrature` crossed with `trunc()` is now refused by design, so
`fuzz_ok()` no longer generates the pair and `fuzz_refusal_cases`
probes the refusal instead. On the previous plan the same comparison
was 28 -> 14.)

UNCONVERGED is not a defect class. The estimate-dependent invariants
(`row_permutation`, `vcov_psd`, `confint_wald`, `simulate_mean`,
`predict_eq_fitted`, `unit_weights`, `loglik_identity`) only say
something about a fit that reached an optimum; when the fit itself
warned about convergence, a violation is downgraded so it does not
compete with real defects. Assembly identities, design agreement,
crashes, refusals, and `finite_loglik` are judged regardless.

---

## Fixed since the previous run

The four items below were the quadrature cluster. All of them trace
back to one property of TMBad's `marginal_gk` transform: it rescales
each integrand **once**, finding the mode and curvature of the
log-integrand by finite differences and baking that `(mu, sigma)` pair
into the tape as a constant, at whatever parameter values the template
held when `MakeADFun` taped it. `adaptive = TRUE` retapes but is
strictly worse here - it fails on cases the frozen path survives.

`R/fit.R` now fits the plain Laplace objective first and tapes the
marginalized one at that optimum (`quad_fit()`), which is both the
calibration the transform wants and the only remaining source of the
conditional modes.

### N1 (fixed). `quadrature = TRUE` left every conditional mode but one at `NA`

`fit$estimates$b` kept its `parList()` value, and under `integrate =`
the random effects are gone from the tape: `parList()` left every
element after the first `NA` and slid an outer value into the first.
`ranef()`, `fitted()` and `predict(newdata =)` all read `b`, so they
returned `NA` for every group but one, with no warning, next to a
correct `logLik`.

The modes now come from the inner Newton solve of the Laplace
objective at the quadrature optimum (`solved_par_list()`), so they are
the modes of the fitted model rather than of the pre-fit. On the
`glmer(nAGQ = 25)` benchmark in `test-v11.R` they match `lme4`'s
`ranef()` to 3e-05 while differing from the plain-Laplace modes by
0.13. Regression tests: `test-quadrature-defects.R`.

### N2 (fixed). Bare `NA/NaN gradient evaluation` under quadrature

poisson, Gamma and Beta over `(1 | ga/gb)`, and Beta even over a
single `(1 | g)`, all satisfied the documented scalar-intercept guard
and died inside RTMB. Taped at the cold start the frozen rescaling
sits far from the real conditional mode, and every family whose
inverse link exponentiates the linear predictor then overflows;
gaussian survived only because its integrand is quadratic wherever it
is sampled. Warm-starting the tape fixes all of them, and also lifts
two negbinomial fits that used to "converge" with a gradient of 8.8
and `NaN`.

What is left of the family surface is a runtime limitation, not a
silent one: see the residual item below.

### N3 (fixed). `quadrature = TRUE` with `trunc()` returned `logLik = +Inf`

Refused at fit time. The truncation normalizer is
`log(F(ub) - F(lb))` over plain CDFs, and the Gauss-Kronrod nodes
reach random-effect values where that difference underflows to exactly
zero while the density itself is still representable. The integrand is
then `+Inf` and the marginalized objective `-Inf` - at the Laplace
optimum as well as at the starting values, so no amount of warm
starting helps. Laplace stays near the mode and never sees it, which
is why `quadrature = FALSE`, REML and profile all fit the same model.

A numerically stable log-difference-of-CDFs would make the pair
correct, but it needs a log-survival and a log-CDF for each of the six
families that carry one, which nothing else in the package needs. The
refusal names the remedy instead.

### N4-adjacent (fixed). `mixture()` under REML and `profile = TRUE`

Both are refused now. A mixture likelihood is invariant to permuting
its components, so it is multimodal in the fixed effects that REML and
`profile = TRUE` integrate out with a Laplace approximation about a
single inner mode. The fits used to stop at `NA/NaN gradient
evaluation` or report a gradient near 1e9 with no guard at all.

Quadrature stays allowed: it marginalizes the random effects, not the
coefficients, and `test-v19.R` pins down that it is exact when the
per-group integrand is univariate. A mixture whose fit collapses a
mixing weight to about `exp(-35)` still degenerates the rescaled
integrand and leaves a gradient near 1e14; the large-gradient and
non-convergence warnings both fire, so it is loud rather than silent,
and the registry records it as conditional.

---

## REAL-NEW

### R1. `quadrature = TRUE` still breaks down on hard likelihoods

7 findings (`fit_error`), all now the same informative refusal naming
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

### R2. `vcov()` is singular in two different ways

2 findings (`vcov_psd`):

```
family=lognormal aterm=cens_chr re=dbar_f special=smooth dpar=dpar_x
  mode=profile  seed 20389865
  vcov() errored: system is computationally singular: rcond = 2.8e-22
family=negbinomial aterm=none re=nested special=mo_int mode=sparse_x
  seed 20323429
  vcov() has non-finite entries
```

Unchanged from the previous run, and the pair is the point. The
profile case takes the `solve(sdr_of(object)$jointPrecision)` branch
of `vcov.frmtmb_fit` (`R/methods-fit.R:216`), which inverts by hand
and lets `solve()` throw; the ML branch reads an already-inverted
`cov.fixed` and returns `NaN` entries instead. Low severity, but the
two paths behave differently for the same condition and neither
message names the model or `diagnose()`.

### R3. `logLik` moves under a row permutation (marginal)

4 findings (`row_permutation`), under `profile` and `sparse_x`. The
designs are elementwise identical after undoing the permutation, so
this is not a design-assembly problem; it is an optimizer landing a
hair away from the same optimum. Recorded because the differences
clear the stated 1e-6 threshold, not because they look like defects.

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

### R5. A model with no free parameters dies inside nlminb

1 finding: `family=poisson re=none mode=profile`, seed 20477795, fails
with `'d' must be a nonempty numeric (double) vector`. This is the
open backlog item glmmTMB#1325/#1317, rediscovered by the generator:
`profile = TRUE` moves every coefficient inside, leaving nothing for
the outer optimizer.

---

## KNOWN-PENDING (rediscovered)

57 findings across 3 of the in-flight items.

| item | findings | how it surfaced |
|---|---|---|
| `(1 \| factor(x))` brms translation | 35 | `brms_translation`: brms rejects `(1 \| factor(gc))`, frmtmb fits it |
| `(f \|\| g)` builds a correlated block | 18 | `brms_re_blocks`: 4 random coefficients per level where brms has 3 |
| truncation ignored by the post-fit surface | 4 | `trunc_support` / `simulate_mean`: `simulate()` draws below the bound |

Note on the truncation item: comparing `simulate()` against
`predict(type = "response")` does **not** see it, because both ignore
the bound in the same direction and therefore agree. The detector that
works is the support check - a simulator must not draw outside the
support the model was fitted under. That is `fuzz_inv_trunc_support()`.

`rescor x cens`, non-estimable rank-deficient prediction, `ar1` with
gapped integer levels, and `anova` across different `nobs` are outside
the generated grammar (single-response specs, full-rank designs,
balanced factor time grids, no model comparison) and were not probed.

---

## GRAMMAR-DIVERGENCE

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
unchanged - the reverse does not. Worth a line in the migration
vignette rather than a code change.

---

## Generator notes

Four generator defects were found and fixed during the first triage;
the run has none.

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

One structural lesson worth keeping: comparing two fits' objectives at
a shared parameter vector is only valid when the two designs are
elementwise identical. mgcv rebuilds a smooth basis in a different
rotation on reordered data - the same model, a different meaning for
`beta` - so `fuzz_inv_permutation()` now checks the design first and
only compares objectives when the parameterization is shared.
