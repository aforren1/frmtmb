# Grammar fuzz findings

Run of the pairwise covering-array fuzzer (`tests/testthat/helper-fuzz.R`,
`tests/testthat/test-fuzz.R`, `dev/fuzz/run-fuzz.R`) on branch `wt-fuzz`
at package version 0.22.0.

| | |
|---|---|
| seed | `20260901` |
| specs | 310 (110 pair-cover rows, 190 3-way pad rows, 10 refusal probes) |
| feasible value pairs | 738, all covered |
| fits | 310 model fits, 310 permutation refits, 75 unit-weight refits |
| wall time | 133 s |
| command | `Rscript dev/fuzz/run-fuzz.R` |
| machine-readable | `dev/fuzz-findings.json` |

Classification counts: 19 REAL-NEW, 84 KNOWN-PENDING, 5
GRAMMAR-DIVERGENCE, 36 UNCONVERGED, 0 GENERATOR-BUG.

UNCONVERGED is not a defect class. The estimate-dependent invariants
(`row_permutation`, `vcov_psd`, `confint_wald`, `simulate_mean`,
`predict_eq_fitted`, `unit_weights`, `loglik_identity`) only say
something about a fit that reached an optimum; when the fit itself
warned about convergence, a violation is downgraded so it does not
compete with real defects. Assembly identities, design agreement,
crashes, refusals, and `finite_loglik` are judged regardless.

---

## REAL-NEW

### N1. `quadrature = TRUE` leaves every conditional mode but one at `NA`

7 findings (`predict_eq_fitted`), across student/gaussian/negbinomial/
binomial/lognormal/zero_inflated_poisson, with `(1 | g)` and
`(1 | ga/gb)`. Highest severity of the run: the numbers are silently
missing, not visibly wrong.

`fit$estimates$b` keeps its `parList()` value for every element after
the first, and under `integrate =` that value is `NA`. `ranef()`,
`fitted()`, and `predict(newdata =)` all read `b`, so they return `NA`
for every group except the first. No warning fires.

```r
set.seed(4)
ng <- 20; nt <- 5; n <- ng * nt
d <- data.frame(g = factor(rep(1:ng, each = nt)), x = rnorm(n))
d$y <- 1 + 0.4 * d$x + rnorm(ng, 0, 0.5)[d$g] + rnorm(n)
f <- frm(bf(y ~ 1 + x + (1 | g)) + gaussian(), data = d, quadrature = TRUE)

sum(is.na(f$estimates$b))   # 19 of 20
sum(is.na(fitted(f)))       # 95 of 100
head(ranef(f)[[1]], 3)      # (Intercept): -0.547, NA, NA
as.numeric(logLik(f))       # -144.69, i.e. the FIT itself is fine
head(predict(f, re.form = NA, type = "response"), 3)   # population level is fine
```

`R/fit.R:303-311` patches the outer components out of `opt$par`
because `parList()` leaves them `NA` under `integrate =`; the inner
`b` gets no equivalent treatment. Either the conditional modes have to
be recovered (a Laplace inner solve at the optimum would do it) or the
post-fit surface has to refuse `re.form = NULL` for a quadrature fit.

### N2. `quadrature = TRUE` dies with a bare `NA/NaN gradient evaluation`

5 findings (`fit_error`) plus 2 more from `ar1` under `profile`/`REML`
(below). Reproducible independently of the fuzz data:

```r
set.seed(4)
ng <- 20; nt <- 5; n <- ng * nt
d <- data.frame(g = factor(rep(1:ng, each = nt)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                x = rnorm(n))
d$ga <- d$g
eta <- 0.5 + 0.4 * d$x + rnorm(ng, 0, 0.4)[d$g] +
  rnorm(ng * 2, 0, 0.3)[as.integer(d$ga) * 2 + as.integer(d$gb) - 2]

for (fam in list(gaussian(), poisson(), negbinomial(),
                 Gamma(link = "log"), Beta())) {
  # ... y simulated from eta for this family ...
  frm(bf(y ~ 1 + x + (1 | ga/gb)) + fam, data = d, quadrature = TRUE)
}
```

| family | `(1 \| g)` | `(1 \| ga/gb)` |
|---|---|---|
| gaussian | ok | ok |
| poisson | ok | `NA/NaN gradient evaluation` |
| negbinomial | ok | ok |
| Gamma | ok | `NA/NaN gradient evaluation` |
| Beta | `NA/NaN gradient evaluation` | `NA/NaN gradient evaluation` |

Every one of these satisfies the documented restriction ("every block
must be a dim-1 us/diag term"), so the guard in `R/fit.R:232-247`
passes and the failure lands inside RTMB instead. Two problems in one:
the Gauss-Kronrod path does not survive several supported families,
and when it fails the user gets an RTMB string with no mention of the
model, the parameter, or `quadrature =`.

The same bare message ends two non-quadrature fits, both `ar1()`
crossed with a special:

```
family=gaussian aterm=weights re=ar1 special=mo_int mode=profile   seed 20483657
family=Gamma    aterm=weights re=ar1 special=smooth mode=reml      seed 20493427
```

`ar1(tt + 0 | g) + s(xs)` on clean data under ML fits without
complaint, so these two look data-conditional rather than structural;
they are listed because the error surface is the same and equally
uninformative.

### N3. `quadrature = TRUE` with `trunc()` returns `logLik = +Inf`

2 findings (`finite_loglik`), weibull and poisson. Reproducible on
clean gaussian data, so it is the pair and not the family:

```r
d3 <- d[d$y > -0.6, ]
f <- frm(bf(y | trunc(lb = -0.6) ~ 1 + x + (1 | g)) + gaussian(),
         data = d3, quadrature = TRUE)
f$opt$objective        # -Inf
as.numeric(logLik(f))  #  Inf
```

The same model under Laplace (`quadrature = FALSE`) has a finite
objective. The fit object is returned with nothing but a gradient
warning, so `AIC`, `anova`, and `confint` all inherit an infinite
likelihood. An unbounded objective should stop the fit, not be
reported as a number.

### N4. `vcov()` under REML or `profile = TRUE` raises a raw LAPACK error

2 findings (`vcov_psd`):

```
family=cumulative aterm=weights re=dbar special=mo_int mode=profile  seed 20296073
  vcov() errored: system is computationally singular: rcond = 4.2e-18
family=student aterm=none re=nested special=gp dpar=dpar_x mode=reml seed 20426991
  vcov() errored: system is computationally singular: rcond = 5.0e-25
```

Both take the `solve(sdr_of(object)$jointPrecision)` branch of
`vcov.frmtmb_fit` (`R/methods-fit.R:216`). The ML branch reads an
already-inverted `cov.fixed` and returns `NaN` entries for a singular
fit; the REML/profile branch inverts by hand and lets `solve()` throw.
Low severity, but the two paths behave differently for the same
condition and the message names neither the model nor `diagnose()`.
The second case fitted without any convergence warning.

**Fixed.** `solve_joint_precision()` (`R/utils.R`) wraps the inversion:
a singular precision now yields `NaN` standard errors plus one warning
naming `diagnose()` and `vignette('diagnostics')`, matching the ML
branch. `hyp_par_cov()` (`R/confint.R`) carried the same raw solve and
uses the wrapper too, so `confint()` and `hypothesis()` degrade the
same way. Regression test in `tests/testthat/test-osa-inference.R`.

### N5. `logLik` moves 2.1e-06 under a row permutation (marginal)

1 finding (`row_permutation`), `family=negbinomial aterm=weights
re=ar1 mode=autoscale`, seed 20292165. The designs are elementwise
identical after undoing the permutation, and neither fit warned, so
this is not a design-assembly problem; it is the autoscale pre-fit
landing a hair away from the same optimum. Recorded because it clears
the stated 1e-6 threshold, not because it looks like a defect.

**Not a defect; closed.** Reproduced on clean `negbinomial + weights()
+ ar1()` data (gap 2.136e-06, the same magnitude). Three observations
settle it:

- The standardization constants are bit-identical under the
  permutation (`sd(x)` and `mean(x)` differ by exactly 0), so
  `autoscale_plan()` is not the source.
- Dropping `autoscale` entirely leaves a gap of -5.26e-06 on the same
  pair, i.e. **larger**. The finding is a property of the ar1 x
  negbinomial surface, not of the two-stage warm start.
- Tightening `nlminb`'s `rel.tol`/`x.tol` to 1e-12 does not move the
  gap. The maximum absolute gradient at the two reported optima is
  1.9e-04 and 1.3e-04; a quadratic with that gradient sits about 1e-06
  above the true optimum, which is exactly the observed spread.

The invariant's 1e-6 threshold is simply tighter than the optimizer's
own convergence on a flat correlated-count surface. No code change.

---

## KNOWN-PENDING (rediscovered)

84 findings across 4 of the in-flight items. Each was found by the
generator, not by a targeted test, which is the evidence the harness
works.

| item | findings | how it surfaced |
|---|---|---|
| `(1 \| factor(x))` errors | 34 | `fit_error`: "couldn't evaluate grouping factor factor(gc) within model frame ... unique() applies only to vectors" |
| truncation ignored by the post-fit surface | 29 | `trunc_support`: `simulate()` draws below the bound, e.g. lb = 1 with 76 of 720 poisson draws at 0; lb = 0.5 with 197 of 1000 weibull draws down to 0.03 |
| `(f \|\| g)` builds a correlated block | 20 | `brms_re_blocks`: 4 random coefficients per level where brms has 3 (lme4 `expandDoubleVerts` gives `(1\|g) + (0+f\|g)`, brms gives one uncorrelated `model.matrix(~f)` block) |
| bar terms crossed with `*` accepted | 1 | `refusal_is_error`: `y ~ 1 + x + x * (1 \| g)` fits instead of failing |

Note on the truncation item: comparing `simulate()` against
`predict(type = "response")` does **not** see it, because both ignore
the bound in the same direction and therefore agree. The detector that
works is the support check - a simulator must not draw outside the
support the model was fitted under. That is `fuzz_inv_trunc_support()`.

Three in-flight items were probed and did **not** reproduce on this
checkout:

- **HSGP basis vs brms** - `gp(xs, k = 10)` matches
  `Xgp_1[Jgp_1, ]` and `Lgp_1` exactly on all 59 gp specs. (An earlier
  iteration of this harness reported a mismatch; that was the fuzzer
  comparing the whole `lp$Z` instead of the gp block's own columns.)
- **`mo()` crossed with a factor** - `mo(xo) * fac` is refused cleanly
  with "mo() interactions support numeric multipliers only: fac", not
  an NA-gradient crash. Kept as a refusal probe.
- **character `cens()` codes** - all 32 `cens(ccl)` specs with
  "left"/"none"/"right" assemble and fit.

`rescor x cens`, non-estimable rank-deficient prediction, `ar1` with
gapped integer levels, and `anova` across different `nobs` are outside
the generated grammar (single-response specs, full-rank designs,
balanced factor time grids, no model comparison) and were not probed.

---

## GRAMMAR-DIVERGENCE

5 findings, all one divergence: **frmtmb accepts `binomial()` without
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

**Documented.** `vignettes/brms-migration.Rmd`, "What changes".

---

## Generator notes

Four generator defects were found and fixed during triage; the final
run has none.

1. The unit-weight metamorphic refit reused the `weights(w)` formula
   (random weights) instead of `weights(one)`, so it "failed"
   everywhere.
2. The brms smooth and gp comparisons sliced the first columns of
   `lp$Z`, which spans every random column of the model, instead of the
   special's own block columns. This produced the spurious HSGP
   mismatch noted above.
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
