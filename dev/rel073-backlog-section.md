## Filed at the 0.73.0 release (2026-10-21)

Four lanes off 0.72.0 (`dev/round-20261021.md`). What each lane, or its
reviewer, filed and did not fix, with the record to read. The items
this round closed are marked "CLOSED at 0.73.0" in place above. Bugs in
other packages are in `dev/upstream-bugs.md` only (brms-30 to brms-35,
RTMB-8, RTMB-9, RTMBdist-2 to RTMBdist-4, rstan-3, Stan-1). Nothing was
posted.

### Nested quadrature (the user wants it in a future round)

Record: `dev/fuzzfix-findings.md`, "To file: nested quadrature, a
design note for a later lane"; `dev/reviews/2026-10-20-fuzzfix.md`,
findings 2 to 4 and "Item 1".

- **An exact route for a row that reads two random-effect blocks.**
  0.73.0 refuses it (`quad_check_blocks()`), as lme4 refuses
  `nAGQ > 1` for more than one term. The exact iterated route is
  `dev/fuzzfix-nested-quad.patch` against 4854e425: `quad_passes()`,
  `quad_check_depth()`, `quad_crossed_stop()`, their call sites in
  `fit_assembled()` and `quad_fit()`, the round-0 `quadrature` docs, and
  the round-0 test as `test-quadrature-nested-exact.R`. On top of 0.73.0,
  remove the `quad_check_blocks()` call first.
  - Accuracy: quadrature minus an iterated Gauss-Hermite integral in R
    over frmtmb's own joint objective, 4.794e-09 (cumulative fuzz spec)
    and 4.505e-08 (Bernoulli); the review's five models 4.1e-09 to
    1.0e-05 in the log-likelihood, SEs equal to a numerical Hessian of
    the exact integral (ratios 0.99998 to 1.00000).
  - Memory is linear in rows at depth 2: Bernoulli 12.3 to 19.4 MB per
    row (1356 MB at 72 rows, 11016 MB at 576); a censored Beta 551 MB
    per row (6612 MB at 12 rows). One block: 405 MB at 576 rows. Depth
    3: 378 s and 49 GB on 32 Poisson rows.
  - The patch orders the passes by level count; a by-split inner factor
    then goes second and costs 5 GB on 12 rows, with a false "crossed"
    refusal at 24. Order by nesting instead: a block goes before another
    when the rows of each of its levels lie inside one level of the
    other (inner first: 1.9 s and 439 MB, 4.7 s and 534 MB, 12.1 s and
    694 MB; `dev/fuzzfix-rev-order.R`).
  - It needs a memory guard that predicts before the iterated tape: the
    Laplace tape's operation count per row, times a per-family factor
    per depth calibrated on the points above, refused above a budget the
    user can set (for example `frmtmb_control(quad_max_gb = )`). Not
    built or measured.
  - The root cause of the old `summary()` failure ("no 'dimnames'
    attribute for array"): one Gauss-Kronrod pass over every effect
    integrated the outer block and left the inner one to the Laplace
    approximation; `RTMB::sdreport()` then warned "Cannot mix inner and
    outer integration methods (yet)" and returned no random-effect
    summary. A route that keeps any effect `random` meets it again.
  - Unmeasured ideas: one tape per outer cluster bounds memory by the
    largest cluster; a smaller inner rule may cut the per-row growth.

### Filed by lane fuzzfix (`dev/fuzzfix-findings.md`, "Defects found, not fixed", "To file"; its review)

#### Open - medium

- **The weights check fails open where it could decide.** Each leaves
  0.72.0's runaway, and none refuses a valid model
  (`dev/reviews/2026-10-20-fuzzfix.md`, "Check after fix round 2"; N3):
  - a `b` prior on only some cell means of an intercept-free
    dispersion (`normal(0, 1)` on `b_sigma_faca` only: logLik 363.151,
    log sigma -36.89, conv 1);
  - `student_t()` on every cell mean (logLik 330.942, log sigma
    -37.54): every column is "free", and the verdict should say so;
  - `student_t()` on Gamma's `shape` or Beta's `phi` (log shape 33.76,
    logLik 352.966; log phi 39.11, logLik 29663.6). The growth is a
    power, about `k^(excess / 2)`, and a Student-t tail is
    `k^-(nu + 1)`, so the verdict follows from `nu`.
- **`mvbf()` with `(1 | obs)` per response runs away without a
  refusal** (two gaussian responses, weights summing above n: logLik
  790.247 on 0.72.0 and 0.73.0). The check runs for one response only.
- **The Laplace allowance for censored rows is derived, not proved.**
  Since the release review's finding 1 the check runs on the uncensored
  rows and, under the Laplace approximation, takes 1/2 off the excess
  per censored row (a censored row's own effect sits where its
  probability turns from 0 to 1, where the curvature grows like
  1 / width). It fails open where that is too much: on 80 gaussian
  rows with 16 right censored, an uncensored excess of 8 is not
  refused, and the unguarded fit ran away (log sigma -13.635, conv 1);
  excesses of 2 and 4 stopped inside (log sigma -0.862, -1.923), 16
  and 32 ran away and are refused (`dev/rel073-wcens-probe.R`,
  `dev/rel073-log/wcens-probe-before.txt`, `-after.txt`). The same for
  asym_laplace with 72 of 80 rows censored: 8 uncensored rows less 36
  is below 0, so it is not refused, and it stopped at logLik 36.51, log
  sigma -7.24 with false convergence (0, 8 and 40 censored are
  refused). A measurement of the censored rows' Laplace term would
  settle the constant.
- **Quadrature stops a runaway at a spurious maximum.** In the same
  probe, `quadrature = TRUE` stopped at log sigma -3.65 to -5.21 (conv
  0) at every excess from 2 to 32, where the exact integral has no
  maximum: TMB's Gauss-Kronrod rule fixes its nodes at tape time
  (TMB-3) and cannot follow a spike narrower than their spacing. The
  check treats quadrature as exact for that reason; a fit that the
  check does not reach (several responses, `se()`) gets that
  spurious maximum.

#### Open - low

- **A quadrature fit with a variance at zero gets a finite Wald SE for
  the sd and an interval below zero**, where the Laplace fit gives NaN
  and the boundary message (`dev/fuzzfix-quad-boundary.R`, seed 3: sd
  4.773e-05, Est.Error 0.1355, Q2.5 -0.2654). `sdr_nan_nonpos()` and the
  boundary rule do not run on a quadrature objective.
- **Weights summing to exactly n with an observation-level effect**
  drift to a boundary dispersion (log shape 15.7 on the fuzz data;
  gaussian `scale = TRUE` log sigma -8.27). Bounded, so not refused; the
  boundary message should reach a dispersion at its limit.
- **`gamma(2, 1)` on a log-scale `b` of `shape ~ 0 + fac`** stops at
  "NA/NaN gradient evaluation" on 0.72.0 and 0.73.0 alike (a
  positive-support prior on a coefficient of either sign; review,
  "Check after fix round 2").

### Filed by lane cdfs (`dev/cdfs-findings.md`, sections 10 and 12; its review)

#### Open - medium

- **Fuzz spec 20478772: `zero_inflated_poisson()` with `trunc(lb = 1)`
  and `zi ~ 1 + z`** returns 12 non-finite vcov entries and warns about
  nothing. Under `trunc(lb >= 1)` the window holds no zeros, so `zi`
  cancels from the likelihood exactly and is not identified. brms fits
  it with `zi` set by its prior alone. Refuse a zero-inflated family
  under `trunc(lb >= 1)`, or warn by name (the review, "Not closed").
- **`asym_laplace()` with a free quantile stops at nlminb's false
  convergence** on 160 rows drawn from it with no censoring, on 0.72.0
  and 0.73.0 alike (`dev/cdfs-dbg14.R`); the likelihood has a kink at
  every observation. Related: lane adhess's filed "no `asym_laplace()`
  fit has a valid Hessian standard error" (0.72.0 section).
- **The cost of a CDF pair per row** (`dev/cdfs-fix1-cost.R`, one sweep
  each): `RTMB::ppois()` 3.00 ms, `poisson_cdf_log()` 207 to 245 ms,
  `log_igamma_two()` 422 ms, `gamma_cdf_log()` 419 ms,
  `nbinom_cdf_log()` 677 ms, so the incomplete gamma and beta pairs are
  140 to 226 times `ppois()`. The poisson pair takes the closed form at
  0 and finite sums to a count of 30, and a zero-truncated poisson under
  importance sampling is back to 0.72.0's time; a large censored or
  truncated fit on the other families pays the full factor.
  - A censored count costs more, most of all under quadrature (the
    release review, finding 4: `dev/rel073-rev-quadmem.R`, seed
    20261021, `(1 | g)` with 6 rows per group, right censored, each run
    alone under a memory cap, peak working set):

    | model, 1440 rows | 0.72.0 | 0.73.0 |
    |---|---|---|
    | poisson `cens()` at 5, Laplace | 0.6 s, 281 MB | 2.6 s, 455 MB |
    | poisson `cens()` at 5, quadrature | 2.3 s, 495 MB | 25.0 s, 2996 MB |
    | gaussian `cens()` at 1, Laplace | 0.6 s, 260 MB | 0.6 s, 267 MB |

    Poisson quadrature at 360 and 720 rows: 0.72.0 399 and 393 MB,
    0.73.0 901 and 1528 MB, about 2 MB per row. A censored Beta under
    single-block quadrature, new at 0.73.0, uses 1272 MB at 180 rows,
    4266 MB at 360 and 5971 MB at 540.

#### Open - low

- **A censored Gamma at shape 3000 stops with "false convergence (8)"
  at the maximum** (logLik equal to the independent one, AD and FD
  gradients agreeing at -0.0116): nlminb's tolerance at a CV of 1.8
  percent. It reported conv 0 on 0.72.0; the warning now reaches a user
  (the review, N6).
- **frmtmb.eam's `wiener()` declares no `cdf_log`**, so it keeps the
  probability-scale path, its own interval workaround and its refusal
  of left or interval censoring with `trunc()`.
- Student-t tails at `nu = 1e12` hold 3.6e-9 relative, and a
  beta-binomial CDF at shapes (1e5, 1e5) 2e-10 to 3.8e-10 (`lgamma`
  cancellation): inputs that do not occur in practice (the review, re-
  check of fix round 1).
- `residuals(type = "osa")` stays refused on a censored discrete fit
  (its window would shift by one; no reference measures it).

### Filed by lane nb2dm (`dev/nb2dm-findings.md`, "To file", "Not done", "Defects found and not fixed"; its review)

#### Open - medium

- **`beta_binomial()` on binomial data stops with "false convergence"
  in 4 of 8 seeds** (the review's `runoff.txt`, seeds 3 to 6): the plain
  `lgamma` difference that `lgamma_shift_diff()` replaced for
  `dirichlet_multinomial()`. At `phi = 1e8` the K = 2
  dirichlet_multinomial agrees with `beta_binomial()` only to 4.22e-07,
  `beta_binomial()`'s own error.
- **The joint-density translator of the gated brms tests
  (`stan_pars_from_fit()`, `tests/testthat/helper-brms.R`) misreads
  uncorrelated per-category group terms**: `multinomial()` with
  `(1 | g)` fails check C by 8.04 log units with a gradient of 2.45, and
  `dirichlet_multinomial()` by 126.7 (seed 65), while `(1 | ID | g)`
  passes. The model is right: `frm_sample()` of the same model agrees
  with brms's posterior. Pre-existing for `multinomial()`.

#### Open - low

- **frmtmb.latent's `hmm()` refuses a `dirichlet_multinomial()` state
  distribution** with its "needs a 'mu' parameter" message: it names
  `multinomial` only. No brms counterpart.
- At the Poisson limit `negbinomial()` and `negbinomial2` can stop at
  different points (seed 4: log shape 16.83, log sigma -17.53; seed 6:
  false convergence against an SE warning). The objectives agree at
  mirrored points to every digit, gradients included
  (`dev/nb2dm-rev-mirror.R`): nlminb's path on a flat ridge.
- brms's `tests.posterior_epred.R` and `tests.stancode.R` call both
  families, and neither file is in the ported-suite ledgers.

### Filed by lane loodraws (`dev/loodraws-findings.md`, "Defects found, not fixed"; its review)

#### Open - medium

- **`posterior_epred()`, `posterior_linpred()` and `posterior_predict()`
  of multivariate draws with `resp = NULL` return the first response's
  S x N matrix**; brms returns S x N x nresp (frmtmb 4000 x 60, brms
  4000 x 60 x 2). The `loo_*()` functions build brms's layout
  themselves.
- **`posterior_linpred()` of draws of an ordinal model with a `cs()`
  term returns S x N**; brms returns one layer per threshold (bin-2 row
  `brmsfit-methods:523`: brms's third dimension is 3). `loo_linpred()`
  inherits it.

#### Open - low

- **`posterior_linpred.frmtmb_draws()` ignores an unknown argument**
  (`posterior_linpred(ds, junk = 1)` returns the matrix): it never calls
  `frm_check_dots()`, where `posterior_epred()` refuses.
- **`fit$data` is NULL on draws**: a `frmtmb_draws` holds `stanfit`,
  `draws` and `fit`, so brms code that reads `fit$data` gets NULL in
  silence (bin-2 row `brmsfit-methods:469` holds hollowly on it).
- `log_lik()`, `posterior_*()` and `loo()` of draws take frmtmb's
  response spelling (`y_a`) only; the `loo_*()` functions and `psis()`
  take both, and `loo_R2()` and `bayes_R2()` brms's (`ya`) only.
