# Lane arcovsample: `log_lik()` and `loo()` for brms's `cov = FALSE` ARMA

Round of 2026-09-28. Gap closed: frmtmb.sample's `log_lik()`, `loo()`,
`waic()`, `psis()` and `loo_compare()` refused every fit carrying an
`ar()`, `ma()` or `arma()` term, including brms's DEFAULT `cov = FALSE`
form that frmtmb 0.64.0 added. That form has a per-observation column
and it is the column brms's own `log_lik()` returns
(`dev/arcov-findings.md`, section 6, which filed this as the most useful
of the remaining refusals). The covariance form stays refused.

Machine: R 4.6.1, brms 2.23.0, rstan 2.32.7, tmbstan 1.2.1, loo. Lane
library `C:/Users/adf44/source/r/wt-arcovsample-lib`; reference build of
the base commit (frmtmb 0.64.0, frmtmb.sample 0.12.0)
`C:/Users/adf44/source/r/rellib-r3`, read-only.

## 1. What the reference build did, measured

`dev/arcovsample-before.R`, log `dev/arcovsample-log/before.txt`. Seed
4021, N = 76 in 8 ragged groups with interior gaps, rows shuffled; NUTS
draws (1 chain, 200 kept) from a `cov = FALSE` `ar(t, g)` gaussian fit
with `(1 | g)` and from an `arma(t, g, p = 1, q = 1)` student fit with
`(1 | g)`. Both cases gave the same answers:

- `log_lik()`: refused. "log_lik() needs a likelihood that factors into
  one term per observation, and the response 'y' does not: its smallest
  independent unit is an R-side residual correlation
  (ar/ma/arma/cosy/unstr) block, so a column of the matrix would be a
  GROUP and leaving one out would drop a whole sequence. brms has no
  family in this position, so there is no leave-one-out convention to
  follow. Compare these models with AIC() on the ML fits, or with
  frm_bootstrap()"
- `loo()`, `waic()`, `psis()`, `loo_compare()`: the SAME message, which
  they reach through `log_lik()`.
- `loo_moment_match()`, `loo_subsample()`, `kfold()`: refused for their
  OWN reasons, which are about having no stored program to refit. Not
  about autocorrelation, and unchanged by this lane.
- `bayes_R2()`: WORKED, a 1 x 4 matrix. It reads `posterior_predict()`
  and `posterior_epred()`, not `log_lik()`.
- `pp_check(type = "dens_overlay")`: worked.
- `pp_check(type = "loo_pit_overlay")`: `rlang_error`, "One of 'lw' and
  'psis_object' must be specified." Section 6: PRE-EXISTING, and on
  every model.
- `posterior_predict()`, `posterior_epred()`: worked, 200 x 76.
- `loo_R2()`: does not exist. Not in frmtmb.sample's NAMESPACE, so there
  was nothing to refuse.

## 2. What changed

The shift `mu* - mu` has ONE definition, `autocor_cond_shift()` in
core's `R/autocor.R`, and the taped objective applies it to `mu` before
any density is taken. The change gives the sampling extension the same
mu through the same function rather than a second implementation.

### frmtmb (core)

- `R/autocor.R`: two new functions beside the existing post-fit helpers.
  `arma_cond_resp(fit)` names the responses whose autocorrelation block
  is brms's `cov = FALSE` form; `arma_cond_dpars(fit, dpv)` takes the
  whole `eval_dpars()` list and returns it with each such response's
  `mu` replaced by `autocor_cond_dpars()`'s one-step mean, leaving every
  other response untouched. Safe to call unconditionally, which is what
  keeps the extension from carrying a branch of its own. One added
  sentence in `?frmtmb-autocor`'s `cov = FALSE` bullet: `log_lik()` and
  `loo()` DO cover this form and refuse the covariance one.
- `R/sampling-api.R`: both names added to the objective-seam section of
  the contract, to `@aliases` and to the `@rawNamespace export()` list.
  The section says WHY they exist and in what order to call them
  (`arma_cond_dpars()` before `row_lpdf()` and before
  `rescor_row_loglik()`, which is the objective's order).
- Nothing else in core changed. The objective, the tape and every
  existing export are untouched, which is why no core behavior moved.

### frmtmb.sample

- `R/loo.R`, `draws_loglik_factors()`: a `cov = FALSE` response is no
  longer refused. The covariance form gets a refusal of its OWN, because
  it has a remedy the generic message could not name: it quotes the
  term's label, says its smallest unit is a whole group, and for
  `ar`/`ma`/`arma` adds that brms's default `cov = FALSE` form does
  factor. `cosy()` and `unstr()` have no such form and the sentence is
  not added for them (asserted both ways).
- `R/loo.R`, `draws_row_loglik()`: one line,
  `arma_cond_dpars(fit, with_cs_offsets(fit, NULL, eval_dpars(fit)))`.
  Everything downstream is unchanged, which is what carries `weights()`,
  `cens()`, `trunc()`, a distributional `sigma`, `student()`'s `nu` and
  `set_rescor(TRUE)`: they enter through `row_lpdf()` and
  `rescor_row_loglik()` exactly as they did.
- `R/loo.R`, `draws_chain_id()`: guard a `NULL` `stanfit`. `%||%` cannot
  do this, because `NULL@sim` is an error and not a `NULL`. See section
  6; this is why the new tests need no sampler.
- `R/loo.R` documentation: a new section on `?sample-log_lik`
  ("Autocorrelation, and what a row conditions on") and one on
  `?sample-loo` ("A time series, and what is left out"). Both say the
  same thing in the place a reader of each will be: a column is the
  density of one row GIVEN THE OBSERVED EARLIER ROWS of its group, so
  PSIS-LOO leaves out one CONDITIONAL density and not the series; every
  retained column still reads `y_t` through its own lagged residuals, so
  it is neither a held-out forecast nor leave-one-group-out. The
  "Likelihoods with no per-observation column" section now names the
  covariance MATRIX rather than autocorrelation in general.
- `R/methods-draws.R`: `posterior_predict()`'s `arma_cond` flag now
  calls `arma_cond_resp(fit)` instead of reading
  `fit$frame$autocor[[resp]]$cov`. Behavior identical; the point is that
  two files no longer answer "is this brms's cov = FALSE form" in two
  places.

- `R/predict-brms.R`, `rescor_row_loglik()`: the Student-t branch now
  goes through `lgamma_shift_diff()` and `log(nu) + log(pi)`. This is a
  BUG FIX with its own measurement, section 6; it is here because
  `set_rescor(TRUE)` with a `cov = FALSE` ARMA term on each response is
  one of the models this lane makes reachable, and it landed on the bug.

`mi()` on the response is still refused, and for its own reason: a
latent response is a parameter, so its column is not an observation's
likelihood. The refusal is `draws_loglik_factors()`'s `mi_map` branch,
untouched, and a test asserts that an `mi()` response WITH an `ar()`
term lands on that message rather than on the autocorrelation one.

## 3. Validation

### (a) The row sum against the taped objective

`dev/arcovsample-validate.R`, log `dev/arcovsample-log/validate.txt`.
Seed 4021, N = 95 in 10 ragged groups with 5 interior gaps, rows
shuffled; NUTS 1 chain, 250 kept draws per model, EVERY draw checked.

This is an IDENTITY, not a measurement: both sides compose `row_lpdf()`
over the same `arma_cond_dpars()` mu, and `build_objective()` is the
bare likelihood plus each random-effect block's own prior. What the run
reports is the residual at full precision, and the content is that the
random-effect correction is the only correction needed.

```
model                              max |rowsum - (-nll - re_prior)|   max rel
gaussian ar(1), no random effect   0                                  0
student arma(1,1) + (1 | g)        1.42108547152e-14                  1.61e-16
gaussian ar(1), sigma ~ x          0                                  0
gaussian ma(1), weights(w)         0                                  0
gaussian arma(1,1), cens(cc)       0                                  0
gaussian ma(1), trunc(lb = 0)      0                                  0
gaussian ar(2) + (1 | g)           1.42108547152e-14                  1.62e-16
```

Zero here is EXACT, not a rounded print: the residual is computed as a
difference of doubles and `max(abs(.))` of it is `0`, which `format()`
would have shown as `1e-16` had it been one. The two non-zero rows are
the two with a random effect, where the correction term is subtracted.

`loo()` ran on all seven and every `elpd_loo` is finite, with finite
pointwise values:

```
model                             elpd_loo      se   p_loo   max k
gaussian ar(1), no random effect  -103.639712   6.94    3.97  0.4925
student arma(1,1) + (1 | g)       -100.270530   6.19   12.70  0.9409
gaussian ar(1), sigma ~ x         -103.673708   6.88    4.55  0.5721
gaussian ma(1), weights(w)        -145.158093  10.47    4.79  0.4632
gaussian arma(1,1), cens(cc)       -89.453361   5.77    4.99  0.5930
gaussian ma(1), trunc(lb = 0)     -112.352130   6.67    3.39  0.4707
gaussian ar(2) + (1 | g)          -100.420372   6.23   12.41  0.9778
```

The `elpd_waic` values in the same order: -103.61967,
-100.12782, -103.61823, -145.10800, -89.40922, -112.31713,
-100.07823.

The two Pareto k above 0.9 are the two random-effect models, where
`p_loo` is about 12.5 on 10 groups: that is the ordinary
leave-one-observation-out-with-b-fixed situation the `?sample-log_lik`
help already warns about, not something the ARMA term causes.

`log_lik()`'s mu against `posterior_epred()`'s, same log: exactly 0 for
all six models where the two are the same quantity (`trunc()` excluded,
where epred is the truncated mean). A MEASUREMENT and not an identity:
`posterior_epred()` reaches the one-step mean through `frm_linpred()` /
`linpred_arma_cond()` and `log_lik()` through `arma_cond_dpars()`, two
call paths. They agree bitwise, which says the two paths perform the
same arithmetic in the same order, not merely to round-off.

### (a2) `set_rescor(TRUE)` with a `cov = FALSE` term on each response

`dev/arcovsample-rescor.R`, log `dev/arcovsample-log/rescor.txt`. Seed
4023, N = 69 in 8 ragged groups, rows shuffled, `ar(t, g)` on both
responses, NUTS 1 chain, 200 kept draws, every draw checked. Here the
row density is the JOINT one (`rescor_row_loglik()`) and the shift has
to reach it, which it does because `arma_cond_dpars()` runs first.

```
model                             max |row sum + nll|   max relative
rescor + ar(1) on both, gaussian  2.84217094304e-14     1.63254808534e-16
rescor + ar(1) on both, student   2.84217094304e-14     1.63190387166e-16
```

The shift is measurably in there: at draw 1 the shifted joint density is
-178.3768838 (gaussian) against -180.0273064 unshifted, and
-176.1550883 (student) against -176.820696. `elpd_loo` is finite in both
(-181.4805285 and -182.5533767).

The student row of that table read `max |row sum + nll| = 47859.8`,
relative 1.83, before the `rescor_row_loglik()` fix of section 6.

### (b) Against brms 2.23.0, row by row

`dev/arcovsample-brms.R`, log `dev/arcovsample-log/brms.txt`. Seed
4022, N = 69 in 8 ragged groups with 3 interior gaps, rows shuffled.

Construction: brms is fitted briefly (1 chain, 300 iterations, 100
kept; its posterior is irrelevant), then EVERY ONE of its draws is read
out with `as_draws_matrix()` and written into a frmtmb draws object by
parameter name. Both packages then report a row log-density at the same
parameter vector, so nothing rests on the two samplers agreeing. The
map is explicit and guarded: `b_*`, `sigma`, `nu` and `r_g[...]` carry
brms's own spelling and scale in frmtmb's stored matrix, `thetaac_k`
takes `ar[k]` then `ma[k]` and `theta_1` takes `log(sd_g__Intercept)`;
any frmtmb column left unfilled stops the run. The transplant is checked
BEFORE the claim rests on it, by `posterior_epred()`, which is a
function of every transplanted parameter.

brms's `log_lik()` undoes its own `order(gr, time)` sort with
`reorder_obs(old_order)` (read from `brms:::log_lik.brmsprep`), so both
matrices are in the user's row order and are compared without
reordering.

```
model                    epred max|d|     ll max|d|    rel      sd(ll)
gaussian ar(1)            0              0             0        0.70107
student arma(1,1)         4.218847e-15   4.440892e-15  2.6e-15  0.66861
gaussian arma(1,1)+(1|g)  1.110223e-15   1.776357e-15  1.1e-15  0.65001
```

Gaussian `ar(1)` agrees BITWISE, cell for cell, over 100 x 69 cells.
The other two agree to 4.4e-15 absolute against cells whose own spread
is 0.65 to 0.70, which is double-precision round-off in two different
orderings of the same sum. The groups' FIRST rows, where brms's
convention (`e_s = 0` before a group starts) could have differed
silently, agree to 1.8e-15, 0 and 0 respectively.

### (c) `loo()`

On the SAME (brms) draws, `loo()` from each package:

```
model                        brms elpd_loo   frmtmb elpd_loo
gaussian ar(1)               -103.10545      -103.10545
student arma(1,1)            -101.17574      -101.17574
gaussian arma(1,1) + (1 | g)  -92.92819       -92.92819
```

That is an IDENTITY given (b) and not a second result: `loo()` here is
`loo::loo.matrix()` in both packages, so equal log_lik matrices give
equal elpds, and the eight printed digits are the aggregated form of the
1e-15 cell agreement.

The LOOSE comparison, each package's `loo()` on its OWN 100 draws, which
is what the task asked for and which is not an identity:

```
model                        brms (se)         frmtmb (se)       diff
gaussian ar(1)               -103.10545 (5.91) -102.40761 (5.78) 0.70
student arma(1,1)            -101.17574 (5.52) -101.72196 (5.57) 0.55
gaussian arma(1,1) + (1 | g)  -92.92819 (5.20)  -97.84677 (5.26) 4.92
```

100 draws from two different samplers on two priors (frmtmb.sample
applies brms's defaults but tmbstan explores a differently parameterized
posterior), so these differences carry no information beyond "the same
order of magnitude": the third is 0.9 of its own reported SE. The
identity in the table above it is the evidence; this table is only
recorded so that nobody re-derives it and reads the 4.92 as a defect.

## 4. Tests

### Pins seen to fail on the reference build

`dev/arcovsample-log/pin-sample.txt`, `pin-core.txt`,
`gated-pin-ref.txt`.

- `extensions/frmtmb.sample/tests/testthat/test-loo.R` against
  rellib-r3: `pass=77 fail=0 err=4 skip=2`, the four new ungated blocks
  each erroring on the old refusal, quoted in the log.
- The same file GATED (`FRMTMB_BRMS_FIT_TESTS=true`) against rellib-r3:
  `pass=83 fail=0 err=5 skip=0`. All five new blocks fail, the brms
  row-by-row one included, each on the refusal.
- `tests/testthat/test-autocor-cond.R` against rellib-r3:
  `pass=46 fail=0 err=1`. This is the WEAK form of a pin: the block
  errors because `arma_cond_resp` does not exist there. The block's
  behavioral content is carried inside it instead, by two assertions
  that would fail if the exports were no-ops: the one-step mu must
  differ from the unshifted mu by more than 0.05 sd(y), and on a
  two-response model only the response with the term may move.
- The frmtmb.sample side has the same kind of guard: one draw sets
  `ar = 0`, where `log_lik()` must equal the plain density, and another
  with `ar = 0.55` must differ from it by more than 5 percent of the
  total log-likelihood. A change that dropped the `arma_cond_dpars()`
  call would pass the first and fail the second.

### Lane build, one file per process, `NOT_CRAN=true`

`dev/arcovsample-log/tests-lane.txt` and the `t-*.txt` files beside it.

```
RESULT lane frmtmb test-autocor-cond.R pass=56 fail=0 err=0 skip=0
RESULT lane frmtmb test-autocor.R pass=192 fail=0 err=0 skip=0
RESULT lane frmtmb test-conditions.R pass=150 fail=0 err=0 skip=0
RESULT lane frmtmb test-bracket-access.R pass=33 fail=0 err=0 skip=0
RESULT lane frmtmb test-portability.R pass=95 fail=0 err=0 skip=0
RESULT lane frmtmb test-adefects.R pass=85 fail=0 err=0 skip=0
RESULT lane frmtmb test-generic-collision.R pass=56 fail=0 err=0 skip=0
RESULT lane frmtmb test-priors-autocor-classes.R pass=62 fail=0 err=0 skip=0
RESULT lane frmtmb.sample test-loo.R pass=103 fail=0 err=0 skip=2
RESULT lane frmtmb.sample test-draws-methods.R pass=148 fail=0 err=0 skip=0
RESULT lane frmtmb.sample test-simulators.R pass=43 fail=0 err=0 skip=3
RESULT lane frmtmb.sample test-message-uniqueness.R pass=6 fail=0 err=0 skip=0
RESULT lane frmtmb.sample test-conditions-census.R pass=11 fail=0 err=0 skip=0
RESULT lane frmtmb.sample test-bracket-access.R pass=1 fail=0 err=0 skip=0
RESULT lane frmtmb.sample test-sampling-ported.R pass=200 fail=0 err=0 skip=4
RESULT lane frmtmb.sample test-brms-pins.R pass=334 fail=0 err=0 skip=0
RESULT lane frmtmb test-families.R pass=252 fail=0 err=0 skip=0
RESULT lane frmtmb test-multivariate.R pass=23 fail=0 err=0 skip=0
RESULT lane frmtmb test-mv-gaps.R pass=57 fail=0 err=0 skip=0
RESULT lane frmtmb test-prior-mv-resp.R pass=36 fail=0 err=0 skip=0
RESULT lane frmtmb.sample test-brms-shapes-draws.R pass=61 fail=0 err=0 skip=0
RESULT gated frmtmb.sample test-loo.R pass=111 fail=0 err=0 skip=0
```

Every skip was read rather than counted (`dev/arcovsample-skips.R`):
test-loo.R's 2 are the two `FRMTMB_BRMS_FIT_TESTS` blocks, which the
gated run then executes with skip=0; test-simulators.R's 3 are
`{frmtmb.latent} is not installed`, which is the lane library's
`dependencies = FALSE` and not a regression.

The whole frmtmb.sample suite and `R CMD check` for both packages are
recorded in section 7.

## 5. Decisions, and what is NOT done

- **brms's latent-residual AR for non-gaussian families is STILL OPEN.**
  Under `cov = FALSE` brms gives a non-gaussian family's AR part N
  latent residuals (`err = sderr * zerr`, added to `eta` with
  ar-weighted lags of itself) and refuses an MA part outright. frmtmb
  refuses both (`dev/arcov-findings.md` section 6), so there is no fit
  for a pointwise density to be taken of, and nothing in this lane
  changes that. It is a different MODEL, not a missing method, and it
  stays open for a later round. A `log_lik()` of it would be
  straightforward once the fit exists: the density is per-row given the
  latent residual, so the columns are observations there too.
- **The covariance form stays refused**, and the message is now true of
  it alone. No approximation was invented: the smallest independent unit
  of `ar(cov = TRUE)`, `cosy()` or `unstr()` is a group, and a
  per-observation column would have to condition on the rest of the
  group, which is not a quantity brms or loo defines.
- **`mi()` stays refused** even with a `cov = FALSE` term, for the
  mi_map reason, which this lane did not touch. The task listed `mi()`
  among the terms to carry; the measurement is that the objective DOES
  carry it (`dev/arcov-findings.md`, `mi() response, arma(1,1)` at
  1e-16) but `log_lik()` refuses the whole model one step earlier, and
  that refusal is about latent responses and not about ARMA. Removing it
  would be a different change with its own argument.
- **`newdata` and `re_formula` stay refused** on `log_lik()`, as they
  are for every model: the pointwise density runs the objective's
  composition over the assembled frame.
- **`loo_moment_match()`, `loo_subsample()` and `kfold()` stay refused**
  for their own reasons, which are about refitting and not about the
  likelihood. They now refuse a `cov = FALSE` fit with their own message
  rather than with the autocorrelation one, which is more informative,
  and no new work is needed for them.
- **No version numbers were changed** and none are proposed here. The
  floor frmtmb.sample needs RISES: it now imports `arma_cond_resp()` and
  `arma_cond_dpars()`, which do not exist in frmtmb 0.64.0. The NEWS
  heading says so in words.
- **`bayes_R2()` was already right and is untouched.** It reads
  `posterior_epred()`, which was already the one-step mean, and section
  3(a) measures that mean as bitwise equal to `log_lik()`'s and section
  3(b) as equal to brms's to 1.1e-15. So `bayes_R2()` on such a fit is
  brms's number by construction; no assertion was added for it.

## 6. Defects found

### Fixed: `rescor_row_loglik()` lost the Student-t density at a large nu

Found by construction, not by a test: `dev/arcovsample-rescor.R` was
written to check the `set_rescor(TRUE)` plus `cov = FALSE` ARMA case and
the student row came back at `max |row sum + nll| = 47859.8`, relative
1.83, against 2.8e-14 for the gaussian one.

Localized to something with NO autocorrelation in it.
`dev/arcovsample-rescor-t.R`, logs `rescor-t-lane.txt` and
`rescor-t-ref.txt`, ran three models on the LANE build and on the
REFERENCE build:

```
model                          lane (before the fix)   reference
gaussian rescor, no autocor    2.842170943e-14         2.842170943e-14
student rescor, NO autocor     NaN                     NaN
student rescor + ar(1) on both 47859.80992             (refused)
```

So it is PRE-EXISTING and has nothing to do with autocorrelation: a
Student-t `rescor` model with no ARMA term at all gave `NaN` on both
builds. The sampled `nu` ranged to 3.05e306, because `nu` is weakly
identified on data with no heavy tail, and `sum(abs(res) > 1e-8)` was
197 of 200 draws; the three that agreed had `nu` in 772 to 19985 and the
197 that did not had `nu` from 1.15e7 up. That is what named the cause.

Characterized numerically, away from any fit:
`dev/arcovsample-tnu.R`, log `dev/arcovsample-log/tnu.txt`. K = 2,
n = 5 rows, C with rho 0.5, seed 7. Three expressions of the same
density: the as-written one, the `lgamma_shift_diff()` one, and the
objective's own `mvt_std_loglik()`.

```
     nu     sum naive      sum stable      objective    naive-obj
  5e+00  -19.10766381  -19.10766381  -19.10766381     3.553e-15
  3e+01  -20.69271468  -20.69271468  -20.69271468     3.553e-15
  1e+03  -21.57991869  -21.57991869  -21.57991869    -4.334e-13
  1e+05  -21.61485249  -21.61485249  -21.61485249    -1.532e-11
  1e+07  -21.61520482  -21.61520481  -21.61520481    -7.970e-09
  1e+10  -21.61521393  -21.61520837  -21.61520837    -5.563e-06
  1e+20 -248.40798177  -21.61520837  -21.61520837    -2.268e+02
 1e+100 -1169.44201896  -21.61520837  -21.61520837   -1.148e+03
 1e+250 -2896.38083871  -21.61520837  -21.61520837   -2.875e+03
 1e+300 -3472.02711196  -21.61520837  -21.61520837   -3.450e+03
 1e+305 -3529.59173928  -21.61520837  -21.61520837   -3.508e+03
 1e+306           NaN  -21.61520837  -21.61520837          NaN
```

The stable form reaches the gaussian limit exactly and stays there; the
as-written one leaves it at `nu = 1e10` and diverges without bound. So
the as-written form is the wrong side, and the fix is one line using a
helper core already has and whose own comment says why it exists
(`lgamma_shift_diff()`, `R/families.R`, with the overflow arithmetic
worked out there). `log(nu * pi)` became `log(nu) + log(pi)` for the
same reason, which that comment also states.

After the fix, all three rows of the localization table read
2.842170943e-14 with 0 of 200 draws disagreeing.

Pinned by a new block in `tests/testthat/test-mv-gaps.R`, beside the
existing `sum(rescor_row_loglik()) == logLik()` block, which passes with
the old code because the ML optimum has a moderate `nu`. The new block
walks `nu` to 1e2, 1e10, 1e20 and 1e300 through its `logm1` link, gives
the objective and the pointwise density ONE parameter vector, and also
asserts convergence (1e20 and 1e300 are the same double) and that `nu`
matters at all (1e2 differs from the limit by more than 1 percent).
SEEN TO FAIL: `pass=53 fail=4` against rellib-r3
(`dev/arcovsample-log/pin-mv-gaps-ref.txt`), `pass=57 fail=0` on the
lane build.

`logLik()`, `frm()` and every maximum-likelihood quantity were never
affected: they read the objective, which was already stable.

### Not fixed

- **`pp_check()`'s `loo_*` types are broken for EVERY model**, not only
  for this one. `dev/arcovsample-ppcheck.R`, logs
  `dev/arcovsample-log/ppcheck-lane.txt` and `ppcheck-ref.txt`: on a
  plain `y ~ x` gaussian fit with no autocorrelation and no random
  effect, on the LANE build and on the REFERENCE build identically,

  ```
  [dens_overlay]    ggplot2::ggplot OK
  [loo_pit_overlay] rlang_error: One of 'lw' and 'psis_object' must
                    be specified.
  [loo_pit]         getvarError: argument "lw" is missing, with no
                    default
  [loo_intervals]   getvarError: argument "psis_object" is missing,
                    with no default
  [loo_ribbon]      getvarError: argument "psis_object" is missing,
                    with no default
  ```

  bayesplot's `ppc_loo_*` functions need the PSIS weights passed in, and
  brms's `pp_check.brmsfit` computes them from `log_lik()` and does so.
  `pp_check.frmtmb_draws` passes neither. Pre-existing, four types
  affected, not specific to autocorrelation, and outside this lane's
  gap. Before this lane it was also unreachable on an ARMA fit for a
  second reason (`log_lik()` refused); now it is reachable and fails the
  same way every other model does, which is why it is recorded here.

- **`nchains.frmtmb_draws()` and `draws_derived_matrix()` have the same
  `NULL`-`stanfit` bug** that `draws_chain_id()` had, and are NOT fixed.
  Seen while probing: `print(VarCorr(ds))` on a `stanfit = NULL` draws
  object dies in `nchains.frmtmb_draws -> %||%`. Nothing in this lane
  reaches them, and a sweep of every `@` read behind a `%||%` is its own
  change with its own test.

- **`loo_R2()` does not exist in frmtmb.sample.** brms has it beside
  `bayes_R2()`, and `bayes_R2()` is implemented here. Not filed as a
  regression, because nothing claims it: recorded because the task asked
  what `loo_R2()` does on such a fit and the answer is "there is no such
  function", which is a gap in the brms surface rather than a defect in
  this one. It would now be one line of `loo()` plus `posterior_epred()`
  on any model whose `log_lik()` works, ARMA included.

- **A shell script must not be edited while bash is reading it.** The
  first core `R CMD check` of this lane died with
  `dev/arcovsample-check.sh: line 30: syntax error near unexpected
  token ')'` and `` `e) PKG=...` ``, because `dev/arcovsample-check.sh`
  was edited (to add a library path) while that run was between
  commands: bash resumes from a byte offset, so it read the tail of the
  shifted file. This is `dev/lane-rules.md`'s "do not edit a runner
  while Rscript is reading it" and it applies to bash the same way. Cost
  one `R CMD build` of core. Recorded so the rule can name both.

- A `%.6f`-style caution that did not apply: the zeros in section 3 were
  checked to be exact zeros rather than small numbers rendered as zero,
  because `format(x, digits = 12)` prints `1e-16` as `1e-16`.

## 7. The suite and the checks

### The whole frmtmb.sample suite, once, one file per process

`bash dev/arcovsample-suite.sh`, log
`dev/arcovsample-log/suite-sample.txt`, per-file output in
`dev/arcovsample-log/suite/`. The totals line is computed from the
RESULT lines by the script, and the script also lists any file that
produced no RESULT line at all:

```
launched 36 files
files=36 pass=1903 fail=0 err=0 skip=15
--- files with no RESULT line (aborted):
(none)
--- BAD lines
  (none)
```

FIVE of those 36 reported `pass=0 fail=0 err=0 skip=0`, which is not a
count. They are the generated brms-suite ports, whose file-level
`skip_unless_brms_suite()` aborts the whole file and records
nothing. Run gated (`dev/arcovsample-log/suite-gated/`) they are:

```
gated test-brms-suite-brmsformula.R  pass=16 fail=0 err=0 skip=0
gated test-brms-suite-brmsterms.R    pass=4  fail=0 err=0 skip=0
gated test-brms-suite-families.R     pass=84 fail=0 err=0 skip=0
gated test-brms-suite-helper-copy.R  pass=33 fail=0 err=0 skip=0
gated test-brms-suite-methods.R      pass=61 fail=0 err=0 skip=0
gated test-loo.R                     pass=111 fail=0 err=0 skip=0
```

`test-brms-suite-methods.R` matters most of the five: brms fixture 1 is
`arma(visit, patient, cov = TRUE)` and feeds 94 assertions, so it is the
file that would notice the covariance-form refusal changing wording. It
passes.

One file is still unrun by design: `test-scale.R`, `pass=0 skip=1`,
"set FRMTMB_SCALE_TESTS=true to run the scale tier". That tier is a cost
benchmark with its own switch and its own runner
(`dev/release/run-scale.ps1`); this lane did not run it.

### `R CMD check --as-cran`

Once per changed package, built WITH vignettes and without
`--no-manual`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and
TinyTeX on PATH: `bash dev/arcovsample-check.sh core|sample`, logs
`dev/arcovsample-log/check-core.log` and `check-sample.log`.

The lane library holds frmtmb and frmtmb.sample alone
(`dependencies = FALSE`), and frmtmb.sample SUGGESTS frmtmb.eam,
frmtmb.latent and frmtmb.ode, so a first attempt stopped at
`checking package dependencies ... ERROR / Packages suggested but not
available`. The reference build supplies those three, READ-ONLY and
AFTER the lane library in `R_LIBS`, so the frmtmb and frmtmb.sample
under test stay this lane's. That is a harness fact, not a finding.

frmtmb.sample: **Status: 1 WARNING, 1 NOTE**, with
`[ FAIL 0 | WARN 15 | SKIP 12 | PASS 1932 ]` in `tests/testthat.Rout`.

Both lines are one environmental cause and NOT this lane's:

- WARNING `checking PDF version of manual`, "LaTeX errors when creating
  PDF version", followed by an EMPTY list of errors. The manual is
  produced: `frmtmb.sample-manual.pdf`, 40 pages, 180442 bytes, and
  `Rdlatex.log` ends in "Done".
- NOTE `frmtmb.sample-manual.tex` left in the check directory, which is
  what that step leaves behind when it is judged to have failed.

Discriminated rather than assumed. `dev/arcovsample-rd2pdf.sh` runs
`R CMD Rd2pdf` with the index on the BASE sources (the main checkout,
read-only) and on the LANE sources, same TinyTeX and same PATH:

```
base exit=0 pdf=187566  dest warnings: 114  LaTeX '! ' lines: 0
lane exit=0 pdf=189793  dest warnings: 114  LaTeX '! ' lines: 0
```

Identical, and both clean. Two further facts point the same way: the
release harness's own run of the base commit
(`dev/release/check.log`) has frmtmb.sample at `Status: OK` with
`checking PDF version of manual ... OK`; and between two runs of the
SAME sources in this lane the companion step
`checking PDF version of manual without index` went from ERROR to OK,
which is a flaky step and not a property of an Rd file. The first run's
log shows TinyTeX generating a font (`mktextfm pcrr8t`) during the
step, which is the usual source of that flakiness. So this WARNING is
the bash-launched check on this box, not an Rd defect, and the
consolidating session's own `R CMD check` is the authoritative one.

frmtmb (core): **Status: 1 NOTE**, with
`[ FAIL 0 | WARN 17 | SKIP 173 | PASS 11840 ]` in `tests/testthat.Rout`.
The NOTE is this machine's standing one, "Skipping checking math
rendering: package 'V8' unavailable", and `checking PDF version of
manual ... OK`. The base commit's own release run of core is also
`Status: 1 NOTE` (`dev/release/check.log` line 101), so core is
unchanged in what it reports.

## 8. Versions, and what this lane does not choose

No version number is changed anywhere, and none is proposed. What the
consolidating session needs to know:

- **frmtmb**: a MINOR bump. Two names are added to the exported sampling
  API, which is new surface, and `rescor_row_loglik()` changes what it
  returns at a large `nu`, which is a fix to an exported function. No
  existing call gets a different answer anywhere the old one was
  correct.
- **frmtmb.sample**: a MINOR bump. `log_lik()`, `loo()`, `waic()`,
  `psis()` and `loo_compare()` accept a class of model they refused, and
  one refusal message changed wording.
- **frmtmb.sample's floor on frmtmb RISES.** It now reads
  `arma_cond_resp()` and `arma_cond_dpars()`, which frmtmb 0.64.0 does
  not export, and it depends on the `rescor_row_loglik()` fix for a
  Student-t `rescor` model's `log_lik()` to be finite. `Imports:
  frmtmb (>= 0.64.0)` in `DESCRIPTION` must become whatever number core
  takes. It was NOT edited here, and neither was the NEWS heading, which
  states the dependency in words instead.

## 9. The scripts, and what each one measured

Every number in this file comes from one of these, run from the
worktree, with its log in `dev/arcovsample-log/`.

- `arcovsample-before.R`: the reference build's refusals (section 1).
  The ONLY script that runs against `rellib-r3` by default.
- `arcovsample-install.R`: roxygenise plus `R CMD INSTALL` into the lane
  library, core then frmtmb.sample.
- `arcovsample-validate.R`: the row sum against the objective, `loo()`,
  and `log_lik()`'s mu against `posterior_epred()`'s (section 3a).
- `arcovsample-rescor.R`: the `set_rescor(TRUE)` plus ARMA case (3a2).
- `arcovsample-rescor-t.R`: localizing the Student-t disagreement to a
  model with no autocorrelation, on both builds (section 6).
- `arcovsample-tnu.R`: the three forms of the multivariate-t density
  against nu, away from any fit (section 6).
- `arcovsample-brms.R`: the row-by-row comparison with brms 2.23.0 and
  the two `loo()` comparisons (3b, 3c).
- `arcovsample-ppcheck.R`: the `pp_check()` `loo_*` types on both builds.
- `arcovsample-rd.R`: renders the changed Rd topics and greps the
  rendered text.
- `arcovsample-rd2pdf.sh`: `R CMD Rd2pdf` on base and lane sources, the
  discriminator for the PDF-manual WARNING.
- `arcovsample-run.R`: one test file, one process, `lane` or `ref`.
  `arcovsample-run-gated.R` is the same with `FRMTMB_BRMS_FIT_TESTS`,
  and `ARCOVSAMPLE_REF=true` points it at the base commit.
- `arcovsample-skips.R` prints a file's skip REASONS, and
  `arcovsample-detail.R` its failure text; a count alone is not one.
- `arcovsample-suite.sh`: the whole frmtmb.sample suite, 8 processes at a
  time, with an aborted-file check and generated totals.
- `arcovsample-check.sh`: `R CMD check --as-cran`, `core` or `sample`.

Added in punch round 1:

- `arcovsample-firstrows.R`: which within-group positions carry the
  shift, and which lag first reaches which position. The measurement
  behind the corrected sentence.
- `arcovsample-rd2.R`: renders the two corrected topics and asserts
  the refuted wording is gone from the RENDERING, not the source.
- `arcovsample-punch1.R`: the mechanical nits, in one pass.

## 10. Punch round 1, 2026-09-29

Review: `dev/reviews/2026-09-29-arcovsample.md`, verdict MERGEABLE with
one blocking documentation error. Reviewer scripts
`dev/arcovsample-rev-*`, logs `dev/arcovsample-rev-log/`.

### The blocking error, and what it really is

`?sample-log_lik` said "with each group's first `max(p, q)` rows taking
the unshifted mean". **False, and measured false.** The reviewer's
`dev/arcovsample-rev-06-firstrows.R` (seed 6161, 4 groups x 7 rows)
found the unshifted positions to be `{1}` at every order it tried
(`ar(1)`, `ar(2)`, `ma(2)`, `arma(2,2)`, `arma(3,1)`), with position 2
shifted in 4 of 4 groups.

Reproduced on this lane's own design and then pushed further, because a
correction needs the RULE and not only the refutation.
`dev/arcovsample-firstrows.R`, log `dev/arcovsample-log/firstrows.txt`,
seed 6162, N = 26 in 5 groups of UNEQUAL length 9, 7, 3, 1, 6, so one
group is shorter than `max(p, q) + 1` at `p = 3` and one has a single
row. Each coefficient is perturbed ALONE and the positions that move are
recorded, which names the lag without reading the recursion:

```
term                      max(p,q)  unshifted  first position each
                                    positions  lag reaches
ar(t, g, p = 1)           1         1          ar[1] -> 2
ar(t, g, p = 2)           2         1          ar[1] -> 2, ar[2] -> 3
ar(t, g, p = 3)           3         1          ar[1..3] -> 2, 3, 4
ma(t, g, q = 2)           2         1          ma[1] -> 2, ma[2] -> 3
arma(t, g, p = 2, q = 2)  2         1          ar,ma[1] -> 2; [2] -> 3
arma(t, g, p = 3, q = 1)  3         1          ar[1..3] -> 2,3,4;
                                               ma[1] -> 2
```

Position 1 is shifted in 0 of 5 groups; every group present at position
2 is shifted, at every order. So the rule is: **lag `i` first reaches
the row at within-group position `i + 1`.** Only the first row of a
group is unshifted; a row at position `k` up to `max(p, q)` carries the
lags up to `k - 1`, that is SOME but not all; from `max(p, q) + 1` on a
row carries all of them. This is brms's `Err` initialized to zero and
then filled from `J_lag[n]`, read directly rather than inferred.

Fixed in `extensions/frmtmb.sample/R/loo.R`'s section
"Autocorrelation, and what a row conditions on", and verified by
RENDERING: `dev/arcovsample-rd2.R`, log `dev/arcovsample-log/rd2.txt`,
greps the `Rd2txt` output and asserts the refuted wording is gone from
the RENDERING rather than from the source.

### The same claim in CORE, which the review did not flag

The refutation applies to pre-existing text this lane did not write, so
it was checked and corrected too rather than left standing:

- `R/autocor.R` `?frmtmb-autocor`, the `cov = FALSE` `\itemize` bullet:
  "A group's first rows get no lagged term" became the rule above.
- `R/autocor.R`'s internal comment above `autocor_cond_positions()`:
  "a group's first rows get no lagged term AT ALL" was the strongest
  form of the false claim and is now the rule.
- `R/compat.R`'s `ar()`/`ma()`/`arma()` group note, which a user reads
  through the compatibility registry: the same correction.

The refusal in `autocor_matrix()` says the form "conditions on each
group's first rows", which is TRUE, and was left alone.

These three are a wording defect shipped in 0.64.0. They are corrected
here because the measurement that refutes them is this round's; the
consolidating session can revert them if it would rather the arcov lane
own them.

### One more thing the rendering caught

The first spelling of the corrected core bullet used
`\eqn{\max(p, q) + 1}`, which `Rd2txt` renders as the literal
`\max(p, q) + 1`. Changed to the two-argument form
`\eqn{\max(p, q) + 1}{max(p, q) + 1}`, which renders as
`max(p, q) + 1`. `tools::checkRd()` gives 0 messages on both changed
topics. This is exactly why the check is on the rendering: the source
read correctly.

### The nits, all five

1. `tests/testthat/test-mv-gaps.R`: the new `test_that()` title was 86
   columns, now 63 (`dev/arcovsample-punch1.R`).
2. Same file: a blank line now separates the new block from the
   `test_that()` after it, as every other pair in the file has.
3. `extensions/frmtmb.sample/tests/testthat/test-loo.R`: the new brms
   block no longer calls `skip_sampler()`. It builds its draws from a
   fixed matrix with `stanfit = NULL`, so tmbstan is not on its path and
   the skip could only have hidden the block on a machine whose tmbstan
   is broken for reasons the block does not touch. A comment says so,
   because the neighbouring block DOES need the pair and a reader will
   wonder why they differ.
4. Core NEWS quoted 5.6e-6 at `nu = 1e10` and 227 log units at
   `nu = 1e20` with no construction named, and the reviewer measured
   1.4e-5 at 1e9 and 2268 at 1e20 on its own design
   (`dev/arcovsample-rev-log/13-rrl.txt`, seed 777, K = 2 with a
   distributional `sigma`), so those figures are design dependent. The
   bullet now states the RULE, that the leading digits of the lgamma
   difference cancel as `nu` grows and overflow to `Inf - Inf` above
   `nu = 3.6e305`, and names both logs for the magnitudes instead of
   quoting one design's numbers as if they were the property.
5. `?frmtmb-autocor` asserted what `frmtmb.sample` does, which is false
   in a mixed install. It now says a POINTWISE log-density exists for
   this form and not for the covariance one, names the two exports an
   extension reads to build it, and sends the reader to
   `frmtmb.sample`'s own `?log_lik` for whether the installed one does.
   No version number, and no cross-package claim.

Nit 6, that `loo()` prints no message about the quantity on a
`cov = FALSE` fit, is NOT changed. brms is silent there too,
`?sample-loo` carries the paragraph, and the `unit` attribute mechanism
exists for columns that are not observations, which these are. Parity is
the argument, and the review makes it itself.

### The floor, with the reviewer's construction

Not fixed here: `frmtmb.sample/DESCRIPTION`'s
`Imports: frmtmb (>= 0.64.0)` is the consolidating session's to raise
and this lane does not choose numbers. What the artifact does today,
from `dev/arcovsample-rev-08-floor.R` (the reviewer's; it installs into
its own scratch library `wt-arcovsample-rev-floorlib`, which this lane
neither created nor wrote to), log
`dev/arcovsample-rev-log/08-floor.txt`:

- `R CMD INSTALL` of this lane's frmtmb.sample against frmtmb 0.64.0
  **succeeds**, exit 0, "testing if installed package can be loaded"
  included.
- `library(frmtmb.sample)` **succeeds**. `NAMESPACE` has
  `import(frmtmb)`, a whole-namespace import, and the explicit
  `importFrom()` names 30 functions, none of them the two new ones, so
  nothing fails at load.
- The failure is LOUD and immediate at the first call:
  `could not find function "arma_cond_resp"` from `log_lik()`, and the
  SAME error from `posterior_predict()`. So the floor bites a method
  that worked on 0.64.0 for this model class and was not part of the
  gap.
- `test-loo.R` on that pairing: **pass=32 fail=0 err=16 skip=2**. 16 of
  the file's blocks error, including every pre-existing `log_lik()` and
  `loo()` block, not only the five new ones.

Nothing silently returns a wrong number; the pairing is unusable rather
than subtly broken. It still must not ship.

### Filed, not fixed

`dev/test-backlog.md` has a new section, "Filed by wt-arcovsample after
punch round 1, 2026-09-29", with four entries: the two
`laplace = TRUE` defects the reviewer found, the `pp_check()` `loo_*`
failure, and the remaining `NULL`-`stanfit` reads.

The two laplace ones were re-run by this lane on BOTH builds
(`dev/arcovsample-log/punch1-laplace-{lane,ref}.txt`, seed 1212) and the
output is byte-identical between them, so neither is caused by this
change: `posterior_predict()` and `posterior_epred()` give 4500
non-finite cells of 4500 for `y ~ x + ar(t, g) + (1 | g)` and 3000 of
4500 for `y ~ x + (1 | g)`, with only a base "NAs produced" warning,
while `log_lik()` on the same object refuses properly; and
`frm_sample(laplace = TRUE)` with no random effect dies with
`Error in -obj$env$random : invalid argument to unary operator`.

### What the review confirmed, which this lane did not re-derive

Recorded because it is evidence this lane need not repeat. The reviewer
wrote its own designs and seeds throughout and every shared count
matched: 7 of 7 installed functions `identical()` to the worktree
source; 15 constructions agreeing with the objective to 7.1e-14
absolute, 8.5e-16 relative, on seeds 5701 and 501, including a group of
ONE row and a mixed `cov = FALSE` / `cov = TRUE` model; 6 further brms
cell-by-cell comparisons on seed 9301 to 5.33e-15, on constructions this
lane did not run (`arma(2,2)`, `ma(2)` with `weights()`, `arma(2,2)`
with `cens()`); and `posterior_predict()` / `posterior_epred()`
`identical()` between the two builds on 6 model shapes with a fixed seed
per call, which is stronger than a tolerance.

The `rescor_row_loglik()` fix is LARGER than this lane reported: on the
reviewer's seed 777 design the base build is wrong at **8 of 15** `nu`
values, by up to **34504 log units**, against 0 of 15 on the lane build.
This lane's `tnu.txt` design (K = 2, n = 5) understated it. The small-nu
end did not regress: 0 to 4 ulp against `mvtnorm::dmvt` at `nu` from
1.05 to 50.

### Tests rerun in this round

One file per process, `NOT_CRAN=true`, lane library first
(`dev/arcovsample-log/p1-*.txt`). `R/compat.R` is in the diff now, so
the compatibility suites run too.

```
RESULT lane frmtmb test-compat.R pass=639 fail=0 err=0 skip=0
RESULT lane frmtmb test-autocor.R pass=192 fail=0 err=0 skip=0
RESULT lane frmtmb test-autocor-cond.R pass=56 fail=0 err=0 skip=0
RESULT lane frmtmb test-mv-gaps.R pass=57 fail=0 err=0 skip=0
RESULT lane frmtmb test-conditions.R pass=150 fail=0 err=0 skip=0
RESULT lane frmtmb.sample test-loo.R pass=103 fail=0 err=0 skip=2
RESULT lane frmtmb.sample test-compat-preflight.R pass=55 fail=0 err=0 skip=5
RESULT lane frmtmb.sample test-message-uniqueness.R pass=6 fail=0 err=0 skip=0
RESULT gated frmtmb.sample test-loo.R pass=111 fail=0 err=0 skip=0
```

`R CMD check` was NOT rerun: it ran once per package in round 0
(section 7), and this round changes documentation, two test-file
cosmetics and one compat string. The consolidating session runs the
authoritative one.
