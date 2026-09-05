# brms post-fit method tests: plan, measurements and findings

Status: implemented on branch `wt-brms-methods`. This document is both
the specification of the tier and its ledger: what was measured, on
which side, and what it means for a user porting a brms script.

## Claim under test

The log-density tier (`tests/testthat/test-brms-likelihood.R`,
`dev/brms-likelihood-tests.md`) proves that frmtmb's objective IS the
Stan program's log density at a point. That is a claim about the
OBJECTIVE. It says nothing about what the two packages return from
`fitted()`, `predict()`, `ranef()`, `conditional_effects()` and the
rest, and before this file nothing in the repository compared any of
them against a brms fit. Every other test that names brms beside a
method calls frmtmb's method next to a brms FORMULA object, and the
only value comparisons against a fitted brms model are
`frmtmb.sample`'s posterior-level ones, which carry an estimand
mismatch and a tolerance.

The claim here is narrower and exact: at the same parameter vector, the
two packages' post-fit methods return the same numbers, and where they
do not, the difference is a property of the METHOD and is written down.

## The mechanism: a brms fit whose draws are frmtmb's estimates

Stan's `Fixed_param` algorithm returns its initial values as the draws.
The log-density tier's translator, `stan_pars_from_fit()`, already
produces the constrained, Stan-named parameter list that `init` takes,
including the non-centered `z` blocks for random-effect models. Feeding
frmtmb's estimates in as `init` therefore produces a `brmsfit` whose
every draw IS frmtmb's estimate, so brms's posterior mean is a point
estimate, its posterior standard deviation is zero, and an exact
comparison is legitimate.

Every draw being identical is the tier's own precondition, not an
assumption: `expect_draws_degenerate()` asserts it per shape.

### The two routes, and the measurement that chose between them

Both were built and measured on row 1 (`y ~ x + z, sigma ~ x`,
gaussian), on Windows 11, R 4.6.1, rstan 2.32.7, brms 2.23.0, with the
C++17 repair the log-density tier documents.

- **(a) `brms::brm(..., algorithm = "fixed_param", init = list(pars),
  chains = 1, iter = 10, warmup = 0, backend = "rstan")`.** Public API
  end to end.
- **(b) `brms::brm(..., empty = TRUE)` for the scaffold, then
  `rstan::sampling(algorithm = "Fixed_param", init = list(pars))` on
  the log-density tier's CACHED `stanmodel`, the result assigned into
  the scaffold's `fit` slot, then `brms::rename_pars()`.**
  `rename_pars()` is exactly what `brm()` itself calls after
  `fit_model()` (`brms:::.fit_model_rstan`, then `x <- rename_pars(x)`),
  and it is exported.

| route | first call on a shape | later calls, same session | needs |
| --- | --- | --- | --- |
| (a) `brm(algorithm = "fixed_param")` | 153.6 s | 0.7 s via `brm(fit = )` | nothing beyond public API |
| (b) scaffold + cached program | 1.04 s | 1.04 s | `brms::rename_pars()`, exported, plus `exclude_pars()` |

Route (b) reaches for exactly one brms internal, `exclude_pars()`, and
only for fidelity: `brm()` passes its result to `rstan::sampling()` as
`pars = , include = FALSE`, and without it the object keeps both
spellings of every group-level block, `r_1[i,j]` beside
`r_Subject[level,coef]`, which is not what a `brm()` fit looks like.
The call is guarded with `exists()`, so a brms that drops the function
loses the tidying and not the tier.

Both produce byte-identical results: `variables()` identical,
`posterior_epred()` identical, and both equal to frmtmb's `fitted()`
exactly.

**Route (b) is the one implemented**, because the 153.6 s is not paid
once. brms compiles through its own machinery, which cannot see the
log-density tier's content-addressed cache, so route (a) recompiles a
program that cache already holds, in every session, for every shape.
Across the 23 shapes of this tier that is roughly an hour of compiling
per cold run, DUPLICATING what the log-density tier compiled from the
same formulas and the same flat priors. Route (b) shares that cache, so
a repository running both tiers compiles each program once. The
`brm(fit = )` reuse in route (a) is a within-session effect only and
does not survive to the next run.

### The mechanism proved before it is used

On row 1, and then on every shape:

- `posterior_epred(x)[1, ]` equals `fitted(fit)`. On row 1 the two are
  `identical()`, bit for bit. That is NOT true of every shape, and an
  earlier draft of this line claimed it was: the zero-inflated row is
  4.44e-16 apart with `identical()` FALSE, and sleepstudy is 5.68e-14
  apart on values between 250 and 450. Machine precision reached by a
  different order of operations, orders inside the 1e-8 this file
  asserts, but not zero.
- `log_lik(x)[1, ]` equals frmtmb's per-row log density, and the row
  sum equals `logLik(fit)` to 12 significant digits (-239.538977525 on
  row 1) for a model with no random effects. With random effects the
  sum is a different quantity: finding 16.

The per-row check is the one that matters: a row-aggregation error
hides inside the summed identity the log-density tier already checks,
so the vector is compared before the sum is.

## Ownership and what was not touched

New files: `tests/testthat/test-brms-methods.R`,
`tests/testthat/helper-brms-methods.R` and this document. Edited:
one line in `.gitignore` for the cache directory, one NEWS bullet, and
`.github/workflows/brms-likelihood.yaml` to run the tier.
`helper-brms.R` and `test-brms-likelihood.R` are CALLED and not
edited. Nothing under `R/` changed.

### CI: one job, one gate, one cache

The tier joins the existing `brms-likelihood.yaml` job rather than
getting a workflow of its own, because the two tiers share the
compiled-model cache: the methods tier builds its brms fits from the
SAME programs the identity tier compiles, through the same
`brms_stan_model()` in `helper-brms.R`. Splitting them would compile
every program twice and would let the two gates drift apart. They share
the gate too, `FRMTMB_BRMS_FIT_TESTS` through
`skip_unless_brms_fit()`, so no second gate was added.

Three things had to change for that to be correct rather than merely to
run:

1. **The cache key hashes all four tier files, not two.**
   `actions/cache` saves at post-job only when the primary key was not
   an exact hit. Hashing the identity tier's two files alone means the
   first run compiles the methods tier's one extra program, restores an
   exact key hit, saves nothing, and recompiles that program on every
   run afterwards. The `md5sum()` vector now carries
   `test-brms-methods.R` and `helper-brms-methods.R` as well.
2. **Two steps, so two PROCESSES.** `test_file()` re-sources the
   helpers, which resets the session-level `.brms_stan_models` cache. A
   program compiled earlier in the same session would then be re-read
   from its own RDS and come back with a DSO that will not initialize,
   which is the "NULL value passed for DllInfo" failure
   `dev/brms-likelihood-tests.md` records under "Cache design". A
   second `test_file()` in the same process would walk straight into
   it.
3. **`MASS` joins `extra-packages`**, because the merged `(1 | q | g)`
   shape draws its group-level values from `MASS::mvrnorm()`. It is a
   recommended package and would almost certainly be present anyway,
   but a missing one would surface as a skip, and the step's own
   `stopifnot` treats a skip with the gate set as a failure.

The workflow's `name:` changed from `brms-likelihood` to `brms-tiers`,
since it is no longer one tier. **A branch protection rule that names
the old check has to be updated with it**, which is the one operational
consequence of this lane outside the test suite.

The cost is one extra Stan compile on a cold cache, the `y ~ x * f`
shape that no identity row carries, plus the tier's own non-Stan work.
The 90 minute timeout is not at risk.

### The C++17 repair is a recipe, not a file

rstan 2.32.7 with StanHeaders 2.39.1 compiles nothing here unless
`CXX17FLAGS` ends in `-std=gnu++17`; `dev/brms-likelihood-tests.md`
diagnoses why under "Toolchain finding", and the workflow's
`Force C++17 for rstan's generated models` step is the same repair in
executable form.

An earlier version of this lane carried the one line as
`dev/makevars-stan-methods`. That file is gone, deliberately. A
checked-in `.mk` cannot be right for two platforms at once: a sibling
lane's copy carries `-mfpmath=sse -msse2 -mstackrealign`, which mean
nothing on the CI runner, so shipping either invites someone to point
`R_MAKEVARS_USER` at the wrong one. To run this tier locally, put

```
CXX17FLAGS = -O2 -Wall -std=gnu++17
```

in `~/.R/Makevars`, or in any file `R_MAKEVARS_USER` points at, and
keep it outside the repository.

## Missing accessors

Reached through the fit object in the test rather than added, because
sibling lanes own `R/` this round.

- **`log_lik()` on a `frmtmb_fit`.** The pointwise log density is
  `frmtmb.sample`'s internal `draws_row_loglik()`; core exports the
  pieces it is built from (`row_lpdf()`, `with_cs_offsets()`) but not
  the composition. `frm_row_loglik()` in `helper-brms-methods.R`
  reproduces it. A user comparing per-row densities with brms has no
  core route to the numbers.
- **`posterior_epred()`, `posterior_linpred()`.** Neither generic has a
  `frmtmb_fit` method. The values are reachable as
  `fitted()` / `predict(type = "response")` and
  `predict(type = "link")`, and this tier proves they are the same
  numbers, but the brms spellings do not resolve.
- **`fitted(newdata = )`.** `fitted.frmtmb_fit` takes `(object, ...)`
  only: no `newdata`, no `scale`, no `dpar`, no `re_formula`. Every one
  of those is reachable through `predict()` instead, so this is a
  spelling gap and not a capability gap.

## Results

23 shapes: the log-density tier's rows 1, 2, 5, 7, 12 (four ordinal
families and `cs()`), 13, 14 (`cens()` and `se()`), 15, 16, 17, 20, 21
(four families), the random-effect variants of rows 7, 16 and the
sleepstudy anchor, and one shape this tier adds, `y ~ x * f`, because
every row of that matrix has numeric covariates only and half of
`conditional_effects()`'s held-value semantics is unreachable without a
factor.

"exact" means a relative difference below 1e-8, the tolerance this file
carries and never widens. Nothing marked exact needed it: the
comparisons land at machine precision, and the tolerance is there for
the last bits of a link round trip, not for a disagreement.

The file as a whole, one file per process, under `pkgload::load_all()`
with `FRMTMB_BRMS_FIT_TESTS=true` and `NOT_CRAN=true`:

| run | tests | assertions | failed | warning | skipped | errors | wall | programs |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| cold, empty cache | 43 | 814 | 0 | 0 | 0 | 0 | 2573.2 s | 23 compiled |
| warm | 46 | 872 | 0 | 0 | 0 | 0 | 97.8 s | 23 restored |

The cold row is from before the review round, which added three blocks
and 58 assertions: the exclusion guard below, the hurdle reproduction
(finding 1b) and the ordinal silence (finding 3). None of the three
compiles anything, so the cold program count is unchanged and the cold
wall time moves by the warm delta.

The cold figure is contended: sibling lanes were compiling Stan
programs on the same machine throughout, and an earlier cold attempt
was killed at 18 of 23 programs after 50 minutes without reaching the
end. The 2573.2 s above is a completed run under lighter load. Read it
as the order of magnitude, not as a benchmark: what it is measuring is
23 rstan compiles, and those are the same 23 the log-density tier needs
for the same 23 shapes.

Warm is what a repeat run costs, and 145.4 s of it is frmtmb fits,
`brm(empty = TRUE)` scaffolds and the methods themselves; the Stan side
of a warm run is 23 RDS reads and 23 `Fixed_param` samplings of ten
identical draws.

| method | result | scope and note |
| --- | --- | --- |
| `posterior_epred()` vs `fitted()` | exact | all 23 shapes, response scale, including the 3-d ordinal and categorical returns |
| brms `fitted()` Estimate vs `fitted()` | exact | all 23 |
| `posterior_linpred()` vs `predict(type = "link")` | exact | 20 shapes; the 3 with a multi-column mu predictor are a shape divergence, finding 13 |
| `posterior_linpred(transform = TRUE)` vs `predict(type = "response")` | exact | 13 shapes; the rest differ by the family's own definition of the mean, finding 14 |
| `posterior_epred(dpar = )` vs `predict(type = "response", dpar = )` | exact | every dpar except a mixture's `theta` (finding 1c) and `se()`'s unused `sigma` (finding 15) |
| `posterior_linpred(dpar = )` vs `predict(type = "link", dpar = )` | exact | every dpar that HAS a linear predictor; a scalar dpar is finding 11 |
| `log_lik()` per row | exact | all 23. The per-row vector, not the sum |
| `sum(log_lik())` vs `logLik()` | divergence | exact without random effects; 75.72 nats apart with them, finding 16 |
| `residuals(type = "ordinary")` vs `type = "response"` | exact | 17 shapes; the nominal and ordinal families have no shared definition |
| `residuals(type = "pearson")` | divergence | same estimand, brms's is a Monte Carlo of it, finding 8 |
| `fixef()` Estimate | exact | every name the two spell in common, all 23 |
| `ranef()` values | exact | 3 random-effect shapes, including the merged `(1 \| q \| g)` block |
| `coef()` values | exact | where both define it; the contract differs, finding 12 |
| `VarCorr()` standard deviations | exact | 3 shapes |
| `hypothesis()` point estimates | exact | 7 expressions on row 1, plus `sd` and `cor` by both spellings, finding 7 for the container |
| `predict()` point column | divergence | different estimand, finding 1 |
| `conditional_effects()`, one-way | exact | 12 shapes, grid and estimate |
| `conditional_effects()`, `conditions = ` | exact | values agree; the conditioning column is absent, finding 4 |
| `conditional_effects()`, `re_formula = NA` | exact | |
| `conditional_effects()`, two-way | divergence | order and held value, findings 2 and 5 |
| `conditional_effects()`, `int_conditions = ` | divergence | ignored, finding 3 |
| `conditional_effects()`, `dpar = ` | divergence | different effect list, finding 6 |
| `conditional_effects()`, `re_formula = NULL` | divergence | different group, finding 6d |
| `conditional_effects()`, ordinal and categorical | divergence | different layout and key, finding 6b |
| `conditional_effects()`, zero-inflated | divergence | the conditional mean, finding 1b |
| `conditional_effects()`, `mo()` | divergence | a continuous grid, finding 17 |
| `conditional_effects()`, nonlinear | divergence | refused, finding 18 |
| `conditional_effects()`, mixture | divergence | refused, finding 1d |
| `pp_check()` | structural | both return a ggplot |
| `loo()`, `bayes_R2()` | structural | frmtmb refuses by design, finding 10 |
| interval and SE columns | structural | different in kind, finding 9 |
| `newdata` | exact | plain, permuted, extra column, dropped factor level, character column; both refuse a missing column and an unseen level |

## The exclusion table, and why it is guarded

Each agreement loop runs over a list of shapes, and a shape that
diverges is kept OUT of that list so the loop asserts agreement where
agreement is the claim. Every divergence is separately pinned by a
block that asserts the CURRENT wrong output against a named
expectation, so a fix produces a labeled failure.

That leaves one hole, and it is not in the assertions. Fix the
zero-inflated defect, leave `r16` and `rC16` out of
`brms_ce_shapes()`, and every assertion still passes while the tier has
quietly stopped covering the shape it was built for. A comment above
the list cannot catch that; only a check can.

So the five lists are DERIVED from one table, `brms_exclusions()` in
`helper-brms-methods.R`. Each row names the finding that put it there
and that finding's class, and each row whose class is **defect**
carries a live-defect probe: `brms_ce_agrees()` for a
`conditional_effects` row, `brms_dpar_epred_agrees()` for a `dpar` row.
The guard block in `test-brms-methods.R` runs every probe and fails the
moment one starts agreeing, with the message

```
finding 1b looks fixed: r16 now agrees with brms, so drop its row from
brms_exclusions() and let brms_ce_shapes() cover it again
```

Paradigm and design-choice rows carry no probe, because there is
nothing to wait for and a probe that can never flip is noise. The
guard also checks that every key names a registered shape, so a typo
excludes nothing silently, and that the defect list is not empty: if it
ever empties, the lists themselves should go, not merely their contents.

The table holds 32 rows: 11 defect, 17 paradigm, 4 design choice. The
probes discriminate, which is the property that matters and is worth
stating as a measurement rather than an intention: all 11 defect rows
return FALSE today, while `brms_ce_agrees()` on `r1`, `r14a` and `rC0`
and `brms_dpar_epred_agrees(r1, "sigma")` all return TRUE.

## Findings

Every finding below ends with a **Class** verdict, because the class is
what the next person needs and prose that leans one way without saying
so is not a verdict.

- **defect** a frmtmb bug. The exclusion that keeps its shape out of an
  agreement loop is temporary, and `brms_exclusions()` in the helper
  carries a probe that fails once the defect is fixed.
- **paradigm difference** right for a maximum-likelihood fit. Nothing
  to change in either package.
- **design choice** neither is wrong; someone has to decide.

26 findings: 10 defect, 8 paradigm, 8 design choice. Twenty-two of them
are the divergences the review classified, and its calls are adopted
unchanged: I found nothing to refute and the mixture evidence
(`theta1` reaching 1.111956599700, above 1 on 1.25% of rows) settles
finding 1c more firmly than my own argument did. The four outside that
list are findings 6 and 6e, both design choices, and findings 9 and 10,
both paradigm differences that were already labeled "NOT a divergence
to fix".

### 1. `predict()` is not the same estimand in the two packages

brms's `predict.brmsfit` summarizes `posterior_predict()`: it DRAWS
from the response distribution and returns
`Estimate, Est.Error, Q2.5, Q97.5`. Its point column is therefore a
Monte Carlo mean and two consecutive calls differ.

frmtmb's `predict(type = "response")` is the conditional mean itself,
returned as a bare numeric vector. It is exactly brms's `fitted()`:

- `predict(fit, type = "response")` equals `fitted(fit)` to 0.
- `fitted(brmsfit)[, "Estimate"]` equals `fitted(fit)` to 0.
- `colMeans(posterior_predict(brmsfit))` against `fitted(fit)`:
  relative 1.09 at 10 draws, 0.0852 at 2000. The gap is Monte Carlo
  error and shrinks with the draw count, so the estimand underneath is
  the same one.

**Which is right: both, for their own paradigm.** brms is Bayesian and
a predictive distribution is a first-class object there; frmtmb is a
maximum-likelihood fit and its `predict()` follows the `stats` and
`lme4` convention, where `predict()` is the linear predictor or its
inverse link and simulation is `simulate()`.

**What a porting user experiences:** `predict(fit)` returns a vector
where brms returned a four-column matrix, so `predict(fit)[, "Estimate"]`
fails loudly. That is the good case. The quiet case is a script that
used `predict()` for predictive UNCERTAINTY: the numbers it gets are the
mean and are correct as a mean, and the spread it wanted is simply not
there.

**Class: paradigm difference.** Neither package is wrong. `predict()` means the predictive distribution in a Bayesian package and the conditional mean in the `stats` and `lme4` lineage frmtmb follows, where the draw is `simulate()`'s job. Nothing to change.

### 1b. `conditional_effects()` plots the wrong mean on a zero-inflated fit

The most consequential finding in this file, and the one the tier was
built to catch.

On `bf(y ~ x, zi ~ x) + zero_inflated_poisson()`, and again on the same
model with `(1 | g)`:

| quantity | value at the first grid point |
| --- | --- |
| brms `conditional_effects()` | 0.6085311262 |
| brms `posterior_epred()` | 0.6085311262 |
| brms `(1 - zi) * exp(eta)` | 0.6085311262 |
| frmtmb `fitted()` / `predict(type = "response")` | 0.6085311262 |
| frmtmb `conditional_effects(method = "predict")` | 0.6085311262 |
| **frmtmb `conditional_effects()`, the default** | **0.7012895769** |
| brms `posterior_linpred(transform = TRUE)`, i.e. `exp(eta)` | 0.7012895769 |
| frmtmb `predict(type = "conditional")` | 0.7012895769 |

frmtmb's DEFAULT `conditional_effects()` curve is the CONDITIONAL mean
of the count component, `exp(eta_mu)`, not the expected response
`(1 - zi) * exp(eta_mu)`. It is 15.2% high at the first grid point and
the excess grows with `zi` along the curve.

**The cause is in the estimate, not in the band.**
`R/conditional-effects.R` in the `method = "epred"`, `band = "wald"`
branch computes

```
p <- predict(x, newdata = nd, type = "link", dpar = dpar, ..., se.fit = TRUE)
df$estimate__ <- lp[["link"]]$linkinv(p$fit)
```

so the point estimate is the inverse link of the MU linear predictor,
taken from the same call that supplies the delta-method standard error.
For every family whose mean IS `linkinv(eta_mu)` that is the expected
response and the curve is right, which is why every other shape in this
tier agrees exactly. For a family whose mean is not the inverse link
of its mu predictor, which covers the zero-inflated families here and
by the same argument the hurdle ones, it is a different quantity.

**frmtmb contradicts itself, and its own documentation.** The `method`
argument is documented as "`\"epred\"` (default): Wald bands for the
expected response". `fitted()`, `predict(type = "response")` and
`conditional_effects(method = "predict")` all return the expected
response and all three match brms exactly. Only the default path does
not.

**Which is right: brms's**, and frmtmb's own three other answers agree
with brms.

**What a porting user experiences:** a plot. No error, no warning, a
smooth curve of the right shape 15% above the mean of the data, and
`fitted()` in the same session disagreeing with it.

**How wrong, across the curve.** 15.2% is the error at the FIRST grid
point and the smallest one on the curve. It grows monotonically:

| grid point | ce `estimate__` | `predict(type = "response")` | ratio |
| --- | --- | --- | --- |
| 1 | 0.701289576884 | 0.608531126228 | 1.152430 |
| 50 | 1.837330 | 1.127389 | 1.629721 |
| 100 | 4.909230 | 1.334783 | **3.677924** |

So the plotted curve is 268% above the expected response at the far end
of the data, not 15%.

**Scope: the hurdle families too, and it needs no Stan.** Every family
in this tier's roster whose mean IS the inverse link of its mu linear
predictor agrees to 1e-16, so the defect is confined to families whose
mean the mu link does not produce. The Stan-backed roster covers
`zero_inflated_poisson` in two shapes. The hurdle families are not in
the log-density tier's matrix, but they do not need to be: this is a
disagreement between two of frmtmb's OWN methods, so it reproduces from
the installed package alone.

Measured on `bf(y ~ x, hu ~ x) + hurdle_poisson()`, 300 rows, seed 13:

```
ce[1]     = 1.25014841313      predict(response)[1] = 1.45040279098
ratio[1]  = 0.86193188603      ratio[100] = 3.37134980095
max |ratio - 1| = 2.37134980095
conditional_effects() == predict(type = "conditional") : 0
```

Two things are worse here than in the zero-inflated case:

1. The error **changes sign**. The curve is 14% BELOW the expected
   response at one end of the grid and 237% above at the other, so it
   is not even monotonically biased, and no eyeball check of the plot
   would read as "too high".
2. The workaround does not exist. `conditional_effects(method =
   "predict")` gives brms's answer on a zero-inflated fit; on
   `hurdle_poisson` it signals "method = 'predict' needs a family with
   a simulator". There is no argument that produces the expected
   response for a hurdle fit.

Pinned by a block in `test-brms-methods.R` that is gated on
`skip_on_cran()` and nothing else, since it needs neither brms nor a
Stan compile.

**Class: defect.** frmtmb's own `fitted()`, `predict(type = "response")` and `conditional_effects(method = "predict")` all give brms's answer; only the default plotting path does not, and the argument documentation calls that path the expected response.

### 1c. A mixture's `theta` is not a probability on frmtmb's response scale

Same class as finding 1b, on a different surface.

On `bf(y ~ 1, theta1 ~ x) + mixture(gaussian(), gaussian())`:

- `family(fit)$links$theta1` is **identity**, so
  `predict(fit, type = "response", dpar = "theta1")` returns the linear
  predictor unchanged, equal to
  `predict(type = "link", dpar = "theta1")` to 0. At the first grid point that is 0.8268843027.
- brms applies the softmax over the component predictors, which for two
  components is `plogis(eta)`:
  `posterior_epred(x, dpar = "theta1")` is 0.6956957293, equal to
  `plogis(eta)` to 2.22e-16.

The two differ by a relative 0.323. frmtmb's number is not a
probability and is not bounded by 1. The likelihood is unaffected: the
log-density tier proves the objective is identical, so the softmax IS
applied inside the density, and it is the `predict(dpar = )` surface
that does not apply it.

`mixture_probs(fit)` is not the missing accessor: it returns the
POSTERIOR class responsibilities given `y` (0.9997218617 at row 1),
which is a different quantity again, 0.744 away from brms's theta.
**frmtmb has no accessor that returns the mixing weight on the
probability scale.**

brms also exposes `dpar = "theta2"` and returns its fixed zero, which
is the reference component's linear predictor rather than its
probability, so brms is not self-consistent here either. frmtmb refuses
the name and lists the dpars it has, which is the better refusal.

**Class: defect.** A scale convention could excuse a different number; it cannot excuse this one. On the fitted data `range(predict(type = "response", dpar = "theta1"))` is [0.606969187849, 1.111956599700] and 1.25% of rows exceed 1. A mixing weight above one is not a value any convention makes right.

### 1d. `conditional_effects()` refuses a mixture whose predictor is on theta

Same model. brms returns one effect, `"x"`, with 100 rows: `x` is in
the model, so it has a panel. frmtmb signals

```
No plottable predictors found for dpar 'mu1'
```

because it looks only at `mu1`'s predictors, and `y ~ 1` has none. With
an explicit `dpar = "theta1"` both packages return `"x"` and then
disagree by the 0.323 of finding 1c.

**Which is right: brms's.** A model with a covariate somewhere has
something to plot, and the default should find it rather than name the
one linear predictor it looked at.

**Class: design choice.** frmtmb chose to enumerate one linear predictor's plottable terms and to say which one it looked at. Widening that search to every dpar is a decision about what `conditional_effects()` means on a model whose mu is an intercept, not a repair of arithmetic.

### 2. `conditional_effects()` rounds the value it conditions on

For a two-way effect `"x:z"`, both packages hold the second predictor
at three values, `mean(z) + c(-1, 0, 1) * sd(z)`.

- brms evaluates at those values exactly
  (`-0.9362169409490, 0.0267302282881, 0.9896773975252` on row 1) and
  rounds only the LABEL it puts in `effect2__`, a factor with levels
  `"-0.94", "0.03", "0.99"`.
- frmtmb evaluates at `signif(mean +- sd, 3)`
  (`-0.936, 0.0267, 0.99`). `ce_second_values()` at
  `R/conditional-effects.R:87` applies `signif(, 3)` to the values
  themselves.

Aligning the two grids on `(x, signif(z, 3))` so that row order is not
what is being measured, the estimates differ by 1.64e-4 on row 1. The
curve frmtmb draws is the model at a slightly different covariate
value, and the size of the error is set by how big `signif(, 3)` is
relative to the coefficient, so it is not bounded by anything the user
can see.

**Which is right: brms's.** A display rounding belongs in the label,
not in the arithmetic. frmtmb has conflated the two.

**What a porting user experiences:** nothing visible. The plot looks
right and the numbers are close. It is only wrong to the precision
someone eventually needs.

**Class: defect.** Display rounding belongs in the label, which is exactly where brms puts it. `signif(, 3)` in `ce_second_values()` changes the point the model is evaluated at, and the size of the resulting error is set by the coefficient rather than by anything the user can see.

### 3. `conditional_effects(int_conditions = )` is accepted and ignored

`int_conditions` is not an argument of
`conditional_effects.frmtmb_fit` and the string does not occur anywhere
under `R/`. It lands in `...` and has no effect on the grid.

Measured on row 1 with `int_conditions = list(z = c(-1, 0, 1))`:

- brms conditions at exactly -1, 0, 1.
- frmtmb returns its default grid, `-0.936, 0.0267, 0.99`, with
  `identical()` values AND `identical()` estimates to the call without
  the argument.

It is not swallowed in silence, and the first draft of this section
said it was, on a measurement taken under `suppressWarnings()`. The
test that pinned it is what corrected the record.
`conditional_effects()` forwards its dots to `predict()`, and
`predict()` warns:

```
ignoring unknown arguments to predict(): int_conditions
```

**And the warning is not there on every shape.** That correction needed
correcting in turn. The warning comes from the `method = "epred"`,
`band = "wald"` branch forwarding its dots to `predict()`. The ordinal
and categorical branches do not forward, so on an ordinal fit:

```
conditional_effects(fo, categorical = TRUE)                -> NO WARNING
conditional_effects(fo, int_conditions = list(x = c(-1, 1))) -> NO WARNING
conditional_effects(fo, nosucharg = 1)                     -> NO WARNING
```

against the gaussian fit's `ignoring unknown arguments to predict():
nosucharg`. So on an ordinal or categorical fit an unknown argument is
discarded in complete silence, which makes finding 6b's ignored
`categorical =` genuinely silent and not merely misattributed.

**Which is right: brms's**, and the argument is not exotic; it is how
brms's own vignettes pick the levels of a moderator. **Class: defect.**

**What a porting user experiences:** a plot conditioned on the wrong
values. On a gaussian fit, a warning that names `predict()`, a function
they did not call, rather than the plot they did: enough to notice, not
enough to explain. On an ordinal fit, nothing at all.

### 4. `conditional_effects()` returns different columns

brms's frame carries the held-constant covariates, `cond__` and
`effect1__` (and `effect2__` for a two-way effect). frmtmb's carries
the varying predictor and the band only:

- brms: `x, y, z, cond__, effect1__, estimate__, se__, lower__, upper__`
- frmtmb: `x, estimate__, se__, lower__, upper__`

`cond__` appears in frmtmb's output only when `conditions = ` is passed
explicitly; brms always has it. This is the data-level root of the
faceting defect recorded in `dev/brms-vignette-audit.md`: brms's
`plot()` calls `facet_wrap(facets = "cond__")`, and on a frmtmb frame
there is no such column to facet on.

**What a porting user experiences:** any code that reads `effect1__`,
or facets on `cond__`, or pulls the held covariate value out of the
frame, gets an error or a missing column.

**Class: defect.** `cond__` is not decoration: brms's own `plot()` method facets on it, so its absence is what makes a ported faceting call impossible rather than merely different. The held covariate values are the other half of the frame's contract.

**Which is right: brms's.**

### 5. The two-way grid is in the opposite row order

brms varies the SECOND effect fastest (`x` repeats three times while
`z` moves); frmtmb varies the first (`z` sits in blocks of 100 while
`x` moves). Both frames hold the same 300 points. No elementwise
comparison of the two is meaningful, and a script that indexed rows
positionally reads a different point.

The order divergence is independent of the rounding one, and the
factor shape separates them: with a FACTOR moderator (`y ~ x * f`)
there is nothing to round, both packages use `sort(unique(f))`, and the
row order still differs, so an elementwise comparison of the `x` column
alone is out by 1.99.

**Class: design choice.** Both frames hold the same 300 points and neither order is wrong. What it costs is that positional indexing does not port, and a package that wanted to match brms would match the order.

**Which is right: neither, but matching brms costs nothing.**

### 6. `conditional_effects(dpar = )` enumerates different effects

Asked for `dpar = "sigma"` on `bf(y ~ x + z, sigma ~ x)`:

- brms returns panels for `x` AND `z`, because it enumerates the
  model's population-level predictors; the `z` panel is constant, since
  sigma does not depend on `z`.
- frmtmb returns a panel for `x` only, because it enumerates the
  predictors of the requested linear predictor.

Values agree exactly on the shared panel. The two also record the dpar
in different places: brms sets `attr(, "response")` to `"sigma"`;
frmtmb keeps `attr(, "response") = "y"` and adds `attr(, "dpar")`.

**Which is right: frmtmb's, on the merits.** A flat panel for a
predictor the dpar does not contain is not information. But
`names(conditional_effects(fit, dpar = ))` differs, so a loop over the
returned list runs a different number of times.

**Class: design choice.**

### 6b. `conditional_effects()` on an ordinal fit: three differences

On `cumulative(y ~ x)` with three categories:

| | brms, default | brms, `categorical = TRUE` | frmtmb, either |
| --- | --- | --- | --- |
| effect name | `"x"` | `"x:cats__"` | `"x"` |
| rows | 100 | 300 | 300 |
| `cats__` | absent | 1, 2, 3 | 1, 2, 3 |
| estimate | expected CATEGORY (1.2555 at the first grid point) | per-category probability | per-category probability (0.7984 at the first grid point) |

1. frmtmb's DEFAULT is brms's `categorical = TRUE` behavior. brms's
   default, the expected category number, has no frmtmb spelling at
   all.
2. `categorical = ` is accepted and ignored by frmtmb: passing `TRUE`
   returns the same frame as passing nothing.
3. The effect is named `"x"` by frmtmb and `"x:cats__"` by brms, so
   `ce[["x:cats__"]]` is NULL on a frmtmb result and `ce[["x"]]` is
   NULL on a brms one.

**Which is right: neither wholly.** brms's default is arguably the
wrong summary for an ordinal model, since an expected category number
treats an ordered factor as if it were a count; frmtmb's default is the
more honest one. But frmtmb reaches it by ignoring the argument that
selects it, so a user cannot ask for the other, and the effect name is
brms's documented key for the categorical layout.

**Class: defect**, for the ignored argument specifically. The default layout is a defensible choice and the effect key is a compatibility decision, but `categorical =` being accepted and doing nothing leaves the user no way to ask for the other layout at all.

### 6c. `conditional_effects(method = )` uses a different vocabulary

brms takes the `posterior_*` generic names
(`"posterior_epred"`, `"posterior_linpred"`, `"posterior_predict"`);
frmtmb takes `c("epred", "predict")`. A ported call with
`method = "posterior_epred"` fails at `match.arg()` with
"'arg' should be one of \"epred\", \"predict\"", which names the
choices but does not say the argument was renamed. Where both spellings
resolve the values agree exactly.

**Class: design choice.** frmtmb's names are shorter and its `match.arg()` failure lists the choices. Accepting brms's spellings as aliases would cost one line and remove a porting stumble.

**Which is right: frmtmb's names, but they should accept brms's.**

### 6d. `conditional_effects(re_formula = NULL)` conditions on different groups

The most consequential of the `conditional_effects()` findings, because
both packages return a plausible curve and neither says which group it
belongs to in a way the other's reader would recognize.

On `Reaction ~ Days + (Days | Subject), sigma ~ Days` (sleepstudy):

- `re_formula = NA`, the default on both sides: identical to 0. This is
  the population curve, and it is the case almost every user hits.
- `re_formula = NULL`: relative 0.181 apart.
  - brms sets the grouping column to **NA**, meaning a NEW, unobserved
    group, and draws that group's random effects from the fitted
    covariance per draw. Its curve is therefore stochastic: 253.28 at
    `Days = 0` with a slope of 7.40 per day on this run, scattered
    around the population 252.86 and 10.09.
  - frmtmb uses the **first observed level**, Subject 308. Its curve is
    249.44 at `Days = 0` with a slope of 20.08 per day, which is
    exactly `fixef()[1] + ranef()[["308", "(Intercept)"]]` and
    `fixef()[2] + ranef()[["308", "Days"]]`.

**Which is right: brms's**, and by a wide margin. A new group is the
only answer that does not privilege an arbitrary level, and brms's
frame carries the `Subject` column so a reader can see it is NA.
frmtmb's frame has no grouping column at all, so it presents one
subject's curve with nothing to say it is one subject's curve, and
which subject depends on factor level order.

The prediction methods do NOT share this problem:
`posterior_epred(re_formula = NA)` matches `predict(re.form = ~ 0)` to
0, and `posterior_epred(re_formula = NULL)` matches
`predict(re.form = NULL)` to 1.94e-16. The divergence is in
`conditional_effects()` alone.

**Class: defect.** Not for choosing a level, but for choosing one silently: the frame carries no grouping column, so the plot presents one subject's curve with nothing to say whose it is, and which subject depends on factor level order.

### 6e. `re_formula` is spelled `re.form` on `predict()`

`predict.frmtmb_fit` takes `re.form`, the lme4 spelling;
`conditional_effects.frmtmb_fit` takes `re_formula`, the brms spelling.
A ported `predict(fit, re_formula = NULL)` reaches `...`, and frmtmb
WARNS: "ignoring unknown arguments to predict(): re_formula". On a
direct call that is exactly right.

The same warning is what surfaces from `conditional_effects()` in
finding 3, because `conditional_effects()` forwards its dots to this
same `predict()`. So an argument dropped by the plotting function is
reported against a function the user never called.

**Class: design choice.** `re.form` is the lme4 spelling and frmtmb is an lme4-lineage package, so the name is defensible; accepting `re_formula` as an alias is the cheap repair.

### 7. `hypothesis()` returns a different object

Point estimates agree EXACTLY on every expression tried, including
linear combinations and dpar coefficients: `x`, `z`, `Intercept`,
`x - z`, `sigma_x`, `sigma_Intercept`, `2 * x + z` all match to 0.

The containers do not:

- brms returns a `brmshypothesis` LIST
  (`hypothesis, samples, prior_samples, class, alpha`) whose
  `$hypothesis` is a data frame with
  `Hypothesis, Estimate, Est.Error, CI.Lower, CI.Upper, Evid.Ratio, Post.Prob, Star`.
- frmtmb returns the table itself, a `frmtmb_hypothesis` data frame
  with `hypothesis, estimate, se, lwr, upr, z, p`.

brms also rewrites the expression (`"x = 0"` becomes `"(x) = 0"`);
frmtmb keeps it verbatim.

**What a porting user experiences:** the documented brms idiom
`hypothesis(fit)$hypothesis$Estimate` reaches frmtmb's `hypothesis`
COLUMN, a character vector, and then errors on the second `$` with
"`$` operator is invalid for atomic vectors". Loud, which is the good
case, but the error names nothing about hypothesis tables.

**Class: design choice.** frmtmb's flat table is the better return for a frequentist fit, and it fails loudly rather than quietly on the brms idiom. A `$hypothesis` alias would remove the stumble.

**Which is right: frmtmb's shape, but the brms idiom deserves an alias.**

### 8. `residuals(type = "pearson")` divides by different quantities

frmtmb divides the response residual by the MODEL's sigma, exactly:
`residuals(fit, "pearson")` equals `(y - mu) / sigma` to 0.

brms divides by the standard deviation of its posterior predictive
DRAWS, a Monte Carlo estimate of the same number. The gap to frmtmb's
answer is 0.687 at 10 draws, 0.091 at 500 and 0.027 at 4000: it shrinks
like a Monte Carlo error rather than staying put, so the estimand is
the same and brms's realization of it is noisy. brms also deprecates
the type.

`type = "ordinary"` against frmtmb's `type = "response"` is exact.

**Class: paradigm difference.** frmtmb's denominator is exact and brms's is a Monte Carlo estimate of the same number, on a type brms itself deprecates. frmtmb is the more accurate of the two and there is nothing to change.

### 9. The interval and standard error columns differ in KIND

Recorded rather than forced into agreement, per the plan.

Under `Fixed_param` brms has nothing to summarize, so `Est.Error` is
exactly 0 and `lower__`/`upper__` collapse onto `estimate__`. frmtmb's
`se__` is a Wald standard error from the observed information, a
frequentist quantity with no draws behind it, and is strictly positive
with a band on both sides of the estimate. The two columns have the
same NAME and the same shape and are not the same thing. That is a
property of the paradigms, not a defect, and the tier asserts the
difference so it stays visible.

**Class: paradigm difference.** Not one of the 22: it is what a frequentist fit has instead of draws, and forcing agreement would mean inventing a posterior.

### 10. `loo()` and `bayes_R2()` refuse, by design

`loo(fit)` and `bayes_R2(fit)` signal errors naming the reason and the
route (`R/loo.R`). brms answers both. This is a documented design
choice and the tier pins that the refusal is a refusal and not a wrong
number, and that brms's answers carry the `SE` column and the
`Estimate/Est.Error/Q2.5/Q97.5` columns that frmtmb has no draws to
fill.

**Class: paradigm difference.** Not one of the 22. The refusal is correct and names the route; what this tier pins is that it stays a refusal rather than becoming a wrong number.

### 11. `posterior_linpred(dpar = )` for a dpar with no predictor

When a dpar HAS a linear predictor the two agree exactly: on row 1,
`posterior_linpred(x, dpar = "sigma")` equals
`predict(fit, type = "link", dpar = "sigma")` to 0.

When it does NOT, as in `y ~ x * f` with no `sigma` formula where brms
declares `sigma` as a scalar on its natural scale, they disagree by a
relative 1.09. brms returns the parameter as declared, on the NATURAL
scale, because there is no linear predictor to be on the link scale of.
frmtmb returns `log(sigma)`, the link scale, because that is what
`type = "link"` means to it.

Both are self-consistent. brms's answer is the one that makes
`posterior_linpred(dpar = )` return the same thing as
`posterior_epred(dpar = )` for a constant dpar, which is arguably the
point of a parameter that is not modeled; frmtmb's is the one that
makes `type = "link"` mean the same thing for every dpar. The
`posterior_epred(dpar = )` comparison agrees exactly on both shapes, so
only the link-scale spelling is affected.

**Class: paradigm difference.** Both are self-consistent, `posterior_epred(dpar = )` agrees exactly on both kinds of dpar, and `type = "link"` meaning the same thing for every dpar is the more defensible rule of the two.

### 11b. `ranef()` and `coef()` key their lists differently

Values agree exactly on sleepstudy (1.38e-16 for the intercepts,
1.29e-16 for the slopes, and 1.25e-16 / 1.99e-16 through `coef()`), and
so do the `VarCorr()` standard deviations (25.478191322 and
5.988822187 on both sides) and the correlation (-0.02805918521).

What differs is the list key and the coefficient set:

- brms keys both lists by the GROUPING FACTOR, `"Subject"`. frmtmb
  keys `ranef()` by the BLOCK, `"Days | Subject"`, and `coef()` by the
  grouping factor, `"Subject"`. So `ranef(fit)$Subject` is NULL on a
  frmtmb fit while `coef(fit)$Subject` is not, in the same model.
- brms's `coef()` broadcasts EVERY dpar's fixed effects over every
  grouping factor: on this model it returns four coefficients per
  subject, `Intercept, Days, sigma_Intercept, sigma_Days`, even though
  sigma has no random effect on Subject and the last two are constant
  across levels. frmtmb returns the two that vary.
- brms names the intercept `Intercept`; frmtmb keeps model.matrix's
  `(Intercept)`.

frmtmb's coefficient set is the more useful one and its `(Intercept)`
is the more standard R spelling. The inconsistent list key between its
own `ranef()` and `coef()` is the part with nothing to recommend it.

**Class: defect**, for the inconsistency inside frmtmb rather than for the difference from brms. Its own `ranef()` and `coef()` key the same model two different ways, so `ranef(fit)$Subject` is NULL while `coef(fit)$Subject` is not. The coefficient set and the `(Intercept)` spelling are frmtmb's to keep.

### 12. `coef()` means different things

brms's `coef()` is group-level: fixed effects broadcast over the levels
of each grouping factor, plus that level's random effect. On a fit with
no group-level effects it has nothing to return.

frmtmb's `coef()` returns brms's quantity when there are random
effects, and falls back to `stats::coef()` otherwise
(`R/methods-fit.R:475`): the bare `mu` coefficient vector when every
other dpar is constant or intercept-only, and the whole `fixef()` list
when one of them is modeled. So `coef()` on a frmtmb fit returns a
named numeric vector for `y ~ x * f`, a list of two vectors for
`bf(y ~ x + z, sigma ~ x)`, and a list of per-group data frames for
`y ~ x + (1 | g)`, and which of the three it is depends on the model.

frmtmb's fallback is the right choice for a package whose fits are
maximum-likelihood and whose users reach for `coef()` the way they do
on a `glm`. What it costs is that the return type is not a function of
the generic alone.

**Class: design choice.** The `stats::coef()` fallback is right for a maximum-likelihood package; the cost is that the return type is a function of the model rather than of the generic, which is worth documenting on the help page.

**Which is right: frmtmb's, with the return type documented.**

### 13. A multi-column linear predictor is one vector to frmtmb

Where the mu predictor has more than one column per observation,
`posterior_linpred()` returns `draws x N x columns` and
`predict(type = "link")` returns the N-vector of the mu predictor
alone:

- categorical `y ~ x` with three categories: brms 10 x 300 x 2, frmtmb
  300. The two columns ARE frmtmb's `predict(type = "link", dpar =
  "mu2")` and `"mu3"`, exactly, so nothing is missing: it is reached
  under a different name.
- `sratio(y ~ x + cs(z))`: brms 10 x 300 x 2 for the two thresholds,
  frmtmb 300. Here the extra columns are NOT reachable through
  `predict()`: `cs()` offsets have no `dpar` spelling.
- `mixture(gaussian, gaussian)`: brms 10 x 400 x 2 for the two
  components, frmtmb 400. Reachable as `dpar = "mu1"` and `"mu2"`.

**Class: design choice.** The categorical and mixture columns are reachable under `dpar = `; only `cs()` has no spelling, and giving it one is a decision about API surface.

**Which is right: frmtmb's, except that `cs()` needs a spelling.**

### 14. `transform = TRUE` is the inverse link, not always the mean

brms's `posterior_linpred(transform = TRUE)` applies the mu link's
inverse and stops. frmtmb's `predict(type = "response")` is the MEAN.
For 13 of the 23 shapes those are the same number to 1e-16. Where they
are not, the difference is the family's own definition:

- `y | trials(n)`: brms returns the probability, frmtmb the expected
  count. `p * n` equals frmtmb's answer exactly.
- zero-inflated: `(1 - zi)` multiplies, finding 1b.
- ordinal and categorical: there is no single mean, and the two return
  different-length objects.

This is a vocabulary difference, not a defect on either side, but it is
the one place where a brms idiom returns a number of the right type and
the wrong magnitude.

**Class: paradigm difference.** The two return different named quantities and both are correct; `p * n` converts exactly between them.

### 15. `se()` leaves a residual sigma each package reports differently

With `y | se(s) ~ x` and no `sigma = TRUE`, the residual standard
deviation beyond the known `s` is zero, and both packages fit that
model: the log-density tier's row 14c is exact and so is this tier's
per-row `log_lik()` comparison. What each REPORTS for the unused
parameter differs.

- brms declares `sigma` and holds it at **0**, so it enters in
  quadrature and switches off.
- frmtmb leaves its `sigma` dpar at the link-scale zero, which the log
  link turns into **1** on the response scale.

Neither number enters a density. A user reading
`predict(fit, dpar = "sigma")` gets 1 where brms gives 0, for a model
whose residual spread is entirely `s`.

**Class: defect**, mildly. Neither number enters a density, but reporting 1 for a residual standard deviation that is not there reads as an estimate, where brms's 0 reads as absent.

### 16. `sum(log_lik())` is not `logLik()` once there are random effects

The trap this tier's per-row comparison exists to avoid, and the reason
the plan asked for a per-row check rather than a summed one.

brms's `log_lik()` is CONDITIONAL on the group-level values in the
draw. Under this tier's mechanism those values ARE frmtmb's conditional
modes, so the row sum is the conditional log-likelihood at the modes.
frmtmb's `logLik()` is the MARGINAL likelihood, the Laplace
approximation with the modes integrated out.

| shape | sum of brms `log_lik()` | frmtmb `logLik()` | gap |
| --- | --- | --- | --- |
| `(1 \| q \| g)` across mu and sigma | -555.805856 | -631.527774 | 75.721918 |
| sleepstudy `(Days \| Subject)` | -814.813566 | -870.000348 | 55.186782 |
| zero-inflated poisson `(1 \| g)` | -418.097315 | -418.097316 | 0.000001 |

The per-row comparison is unaffected: frmtmb's row density is
conditional too, and it matches brms to 1e-16 on all three. The third
row's near-zero gap is a property of that fit, not of the method, which
is why the tier asserts the inequality and the two large gaps rather
than a number.

Nothing here is wrong in either package. What is wrong is porting
`sum(log_lik(fit))` as a stand-in for `logLik(fit)`, and it would be
wrong by 75 nats.

**Class: paradigm difference.** The Laplace correction is the frequentist quantity and `logLik()` is right to carry it. Nothing to change in either package; what is wrong is the port.

### 17. `conditional_effects()` gives a `mo()` predictor a continuous grid

A monotonic effect is defined at the ordered LEVELS of its variable and
nowhere between them: the simplex assigns one increment per step.

- brms plots the four levels: 4 rows, `inc` in 0, 1, 2, 3.
- frmtmb builds the same 100-point numeric grid it builds for any
  numeric predictor, and evaluates the monotonic effect at 0.0303,
  0.0606 and so on.

At the four points both define, the estimates agree exactly, so the 96
interpolated points are the whole of the difference. The other effect
in the same model, the plain numeric `z`, agrees exactly too.

**Which is right: brms's.** The interpolated values are not wrong
arithmetic, they are the model evaluated where it has no meaning.

**Class: defect.** The interpolated points are not wrong arithmetic, they are the model evaluated where a monotonic effect has no definition.

### 18. `conditional_effects()` refuses a nonlinear predictor

brms plots `bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE)` with no
special argument. frmtmb refuses:

```
conditional_effects() cannot put a wald band on a nonlinear predictor:
predict() has no standard error for it. Use band = "boot", ...,
method = "predict", ..., or display one nonlinear parameter with
dpar = "a"
```

The refusal is deliberate and the reason is real, since the default
band is a delta-method interval and there is no analytic standard error
here, and the message names three ways out. `method = "predict"` produces the
curve, and it is brms's curve to 1e-16.

What it costs is that the brms call does not port: the same expression
that draws a plot there signals an error here. This is the
`conditional_effects` nonlinear guard already recorded in
`dev/brms-vignette-audit.md`, now pinned against brms rather than
against a reading of the source.

**Class: paradigm difference.** No analytic delta-method standard error exists for a nonlinear predictor, the message names three ways out, and one of them reproduces brms to 1e-16. The refusal is right.

**Which is right: frmtmb's refusal, though the brms call does not port.**

## Deliberate omissions

- **Row 3, `mo()` interactions.** A sibling lane is changing frmtmb's
  monotonic semantics this round, so no shape here spells
  `mo(x) * z` or `mo(x):z`. Row 2, `mo(inc) + z`, is present and is the
  spelling both packages already agree on.
- **Rows 4, 6, 8, 9, 10, 11, 18, 19** of the log-density tier's matrix.
  The tier's own translator has no rule for them, so no fixed-parameter
  fit can be built for them here either. They join this tier the same
  day they join that one.
- **`residuals()` for the categorical and ordinal rows.** frmtmb
  refuses a nominal response outright and scores the ordinal
  categories; brms does neither. `brms_resid_shapes()` names the
  exclusion instead of skipping inside a loop.
