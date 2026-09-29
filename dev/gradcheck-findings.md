# Lane wt-gradcheck: the convergence warning on the gradient

The claim under test: "Large maximum absolute gradient at the optimum"
fires on correct fits. It does. On 720 correct fits over nine designs at
four sample sizes the 0.64.0 build warned on 289 of them, and the rate
rises with the sample size, from 0.111 at 200 rows to 0.750 at 20,000.
The threshold was absolute, and what it was compared against is not.

Everything below was measured on this box with R 4.6.1. "base" is the
read-only reference build of 0.64.0 at `C:/Users/adf44/source/r/rellib-r3`
and "lane" is this worktree installed into
`C:/Users/adf44/source/r/wt-gradcheck-lib`. Every script takes the arm as
its first argument and every log is under `dev/gradcheck-log/`.

## 1. What changed

`R/fit.R`

- `fit_outer_box()` (new). The box the optimizer ran under, over the
  outer parameter vector. `set_prior(lb =, ub =)` is the only way to put
  a bound on a frmtmb fit, and that box is merged in after `fit$lower`
  and `fit$upper` are filled, so those slots never hold anything and the
  resolved box is now recorded in `fit$cache$bounds` by
  `fit_assembled()`. Without that the check could not see the exact
  bounds the arcov lane reported false alarms on. The reviewer showed the
  `resolve_bounds()` fallback is unreachable, because `frm()` has no
  `lower`/`upper` argument and every internal caller forwards
  `fit$lower`; it stays as the open box an object assembled elsewhere
  should get, and the comment now says so.
- `grad_bound_active()` (new). Which gradient components a bound holds in
  place, by the KKT sign test, with "on the bound" judged relative to the
  bound's own magnitude (`sqrt(.Machine$double.eps)` times it).
- `grad_headroom()` (new). `0.5 * g' H^-1 g` over the parameters no bound
  holds: the log likelihood one exact Newton step from the stopping point
  would still buy. `fit$obj$he()` where RTMB offers it, otherwise
  `stats::optimHess()` on the gradient. `NA` when the Hessian cannot be
  built or is not positive definite there.
- `grad_verdict()` (new). The two stages in order, cached on the fit. It
  reports `proj_par`, the parameter the projected gradient belongs to,
  which is the one the verdict is ABOUT: on a bounded fit the argmax of
  the raw gradient is the parameter the projection removed.
- `grad_warning_msg()` (new). The warning text. The opening phrase
  "Large maximum absolute gradient" is kept verbatim because
  `helper-fuzz.R:1763` greps it; the NUMBER after it is the maximum
  absolute gradient, which is what the phrase says, and the projected
  gradient, its parameter and the headroom follow in their own clauses.
  Punch item 2.
- `check_convergence()`. Evaluates the gradient once as a VECTOR and
  hands it to `grad_verdict()`; the `else if` arm now warns on the
  verdict instead of on the absolute gradient. The importance branch
  above it is untouched.
- `frmtmb_control()`'s `@param grad_tol` rewritten; `@srrstats {RE3.0}`
  and `{RE3.2}` updated.

`R/confint.R`

- `diagnose()` calls `grad_verdict()` and returns four new fields,
  `grad_proj`, `grad_proj_par`, `grad_bound_held` and `grad_headroom`.
  `max_grad` and `worst_grad` are untouched, because seven extension test
  suites read `max_grad` as a ratio to the log likelihood.
- NOT on an importance-corrected fit: the verdict is skipped there, the
  four fields are `NA`, and the printed block says why (punch item 5).
  Where there is no verdict the OLD gate decides the clean line, which is
  the raw gradient against `grad_tol` (punch round 2, item 11).
- The printed block names the bounds that hold a gradient, the parameter
  the projected gradient belongs to, the headroom when one was measured,
  and, on a standardized fit, which lines are on which scale.
- The "No convergence problems detected" condition used a hardcoded
  `out$max_grad < 1e-3`, which ignored `frmtmb_control(grad_tol =)`
  entirely. That was a defect in its own right. It now uses the fit's own
  verdict where there is one, and that same gate with `grad_tol` honoured
  where there is not.

`R/autoscale.R`

- The comment at `autoscale_prefit()` justified muffling the pre-fit's
  warnings partly by a FALSE one: a bound on the centered intercept made
  the pre-fit warn on a clean fit. That case is gone at its source now,
  so the comment says the muffling stands for the duplication alone.
- `autoscale_small_sd`'s comment now names the trip-wire explicitly,
  because that comment records the one gap the trip-wire leaves (section
  6).

`tests/testthat/test-grad-verdict.R` (new), `NEWS.md`,
`vignettes/diagnostics.Rmd`, `man/diagnose.Rd`,
`man/frmtmb_control.Rd`, and `dev/test-backlog.md` (the bound-aware
covariance defect, punch item 6).

`optimize_obj()`, `quad_fit()` and `importance.R` were NOT changed. See
section 7.

`docs/` was not touched. The rendered pkgdown copies still carry the old
sentence (`docs/articles/diagnostics.md:36`,
`docs/reference/frmtmb-scales.md:146`,
`docs/frmtmb.coupling/articles/coherence.md`) and the release's docs
rebuild is what replaces them; `docs/` is in `.Rbuildignore` and no check
reads it.

## 2. The criterion

Two readings of `grad_tol`, in this order, because they cost different
amounts.

1. The gradient trip-wire, free: the largest absolute gradient component
   in the natural units `autoscale` works in. A fit under `grad_tol` says
   nothing. This is exactly the old test.
2. The verdict, one Hessian: the components a bound holds are dropped,
   and the fit warns only if `0.5 * g' H^-1 g` over the rest exceeds
   `grad_tol` as well, now read in log-likelihood units. An unusable
   Hessian counts as a confirmation, so the check never goes quiet
   because a measurement failed.

Why the Newton decrement and not the alternatives:

- It is what the user wants to know. It is the log likelihood still on
  the table, and 1e-3 of those cannot move a standard error, an AIC
  difference or a likelihood-ratio p-value.
- It is invariant under reparameterization, including the diagonal
  rescaling `autoscale` applies, so it does not have to be told which
  units the gradient is in.
- It does not grow with the sample size, which is the mechanism behind
  the false alarms (section 4).
- It is calibrated, and that was measured rather than asserted (section
  5).

Rejected, with the measurement:

- **A relative gradient, `max|g| / max(|f|, 1)`.** Handles the sample
  size but not the bounds: on the two bounded constructions it reads
  0.190 and 0.370 while the fits are exactly right
  (`dev/gradcheck-01-construct.R`, seeds 101 and 102). It also misses a
  real stall on a large data set, because a gradient of 1 against an
  objective of 1e5 is 1e-5 whatever is left on the table.
- **`grad_tol^2` as the headroom tolerance** (the form
  `frmtmb.latent::hmm_starts()` uses, and dimensionally the natural
  companion of a gradient tolerance). 1e-6 at the default, and measured
  too tight: a correct mixed fit with a collapsed variance component has
  a headroom of 2.105e-6 at 20,000 rows and 3.348e-5 at 400,000
  (`dev/gradcheck-03-design.R` design `nearzero`, seed 413459;
  `dev/gradcheck-06-truepos2.R` U7). The near-singular direction inflates
  `H^-1` there.
- **`grad_tol^2 * max(|f|, 1)`**, the same extension's relative form. Too
  loose the other way: it would miss a poisson GLMM that is 8.7e-4 short
  with an objective of 992 (tolerance 9.9e-4) while catching an ordinal
  fit that is 6.3e-2 short. A tolerance that moves with the objective
  makes the verdict depend on the sample size again, in the opposite
  direction.
- **`max|g_i| / sqrt(H_ii)`**, the diagonally scaled gradient. Costs the
  same Hessian and separates the classes about as well, but it is not
  invariant under reparameterization and it does not name a quantity the
  user can act on. It is worth recording that it is the ONE candidate
  that sees the defect in section 6: 40.0 where every correct fit reads
  below 1e-2.

## 3. The false alarms, constructed

`dev/gradcheck-01-construct.R`, logs `dev/gradcheck-log/01-base.txt` and
`01-lane.txt`. Eight of the eleven rows warned on base; one warns on
lane, and that one is a fit whose user set `grad_tol = 1e-12`.

| construction | seed | max\|grad\| | projected | headroom | base | lane |
|---|---|---|---|---|---|---|
| `ub=0.1` on `b`, gaussian n=300 | 101 | 121.2 | 4.9e-07 | 1.6e-15 | warns | silent |
| `ub=0.1` on `ar`, `cov = TRUE` | 102 | 374.1 | 6.1e-04 | 2.0e-10 | warns | silent |
| cumulative, n = 6000 | 103 | 7.76e-03 | same | 9.1e-09 | warns | silent |
| cumulative, n = 20000 | 1030 | 5.54e-03 | same | 3.0e-09 | warns | silent |
| gaussian, n = 50000 | 1040 | 6.45e-03 | same | 1.9e-09 | warns | silent |
| gaussian, n = 200000 | 1040 | 5.84e-02 | same | 2.2e-08 | warns | silent |

The evidence that each fit is correct, one construction per row:

- **Seed 101.** The estimate of `x` is exactly 0.1, the bound. The
  objective at `x` = 0.1, 0.099, 0.09, 0.05 and 0 is 637.543786510,
  637.665065800, 638.759450250, 643.685901640, 649.987507844, monotone
  increasing away from the bound, so the boundary point IS the
  constrained maximum. A refit from it with five restarts moves the
  objective by 3.4e-11 and no parameter by more than 3.3e-07. The
  gradient of 121.2 is the KKT multiplier and nothing else.
- **Seed 102.** The same shape on `ar(cov = TRUE)` with `ub = 0.1`,
  which is the arcov lane's second row.
- **Seed 103.** `logLik` is -7448.63737723 against `MASS::polr`'s
  -7448.63737722, a relative gap of 1.223e-12; the two slopes agree to
  1.871e-06 relative. An independent implementation reaches the same
  optimum.
- **Seed 103 again**, re-optimized from its own optimum with eight
  restarts and `grad_tol = 1e-12`: the objective moves by -8.8e-09 and no
  parameter by more than 1.5e-06. The optimizer will not leave the point.
- **Seed 1040.** Coefficients agree with `lm()` to 5.2e-07 and 4.1e-07
  relative at n = 50,000 and 200,000. The closed form is exact.

The arcov lane's figures reproduce in kind, not in value: it reported
16.7 and 32.9, this lane's designs give 121.2 and 374.1. The mechanism is
the same and the magnitude is a property of the design, so nothing in the
arcov record is wrong.

## 4. The false-alarm rate on a design set

`dev/gradcheck-03-design.R`, one process per group of three designs, 20
replicates per cell, four sample sizes, one seed per (design, n,
replicate) and no seed shared between cells. Summarised by
`dev/gradcheck-05-summarise.R` into the block below, which is pasted from
its output. Log: `dev/gradcheck-log/05-both.txt`.

```
rows per build: 720  distinct seeds: 720  fit errors: 0 0

Every fit in the design set, both builds, same seed:
  objectives bitwise identical : TRUE
  max|grad| bitwise identical  : TRUE
  optimizer codes identical    : TRUE
  largest objective difference : 0

design           n  reps  base  lane  proj  head    max gmax   max gproj    max head
binomial       200    20     0     0     0     0   8.645e-04   8.645e-04   8.910e-09
binomial      1000    20     0     0     0     0   9.587e-04   9.587e-04   5.261e-09
binomial      5000    20     0     0     0     0   6.653e-04   6.653e-04   5.181e-10
binomial     20000    20     4     0     4     0   1.264e-03   1.264e-03   3.946e-10
bounded        200    20    20     0     0     0   6.633e+01   2.461e-04   2.575e-10
bounded       1000    20    20     0     6     0   2.860e+02   2.478e-03   1.538e-09
bounded       5000    20    20     0    16     0   1.395e+03   9.939e-03   9.567e-09
bounded      20000    20    20     0    20     0   5.348e+03   1.547e-02   3.018e-09
cumulative     200    20     0     0     0     0   8.361e-04   8.361e-04   7.013e-09
cumulative    1000    20     2     0     2     0   2.183e-03   2.183e-03   7.929e-09
cumulative    5000    20    14     0    14     0   3.849e-03   3.849e-03   6.841e-09
cumulative   20000    20    17     0    17     0   7.941e-03   7.941e-03   3.759e-09
gaussian       200    20     0     0     0     0   9.302e-04   9.302e-04   4.351e-09
gaussian      1000    20     7     0     7     0   3.141e-03   3.141e-03   5.443e-09
gaussian      5000    20    18     0    18     0   3.058e-02   3.058e-02   4.848e-08
gaussian     20000    20    20     0    20     0   5.472e-02   5.472e-02   4.771e-08
nearzero       200    20     0     0     0     0   8.127e-05   8.127e-05   2.091e-08
nearzero      1000    20     0     0     0     0   6.818e-04   6.818e-04   1.014e-07
nearzero      5000    20     0     0     0     0   8.291e-04   8.291e-04   4.090e-07
nearzero     20000    20     5     0     5     0   1.136e-02   1.136e-02   2.105e-06
nonlinear      200    20     0     0     0     0   4.544e-04   4.544e-04   1.483e-10
nonlinear     1000    20     0     0     0     0   9.867e-04   9.867e-04   6.642e-11
nonlinear     5000    20     8     0     8     0   3.347e-03   3.347e-03   9.899e-11
nonlinear    20000    20    18     0    18     0   1.997e-02   1.997e-02   8.132e-10
poisson        200    20     0     0     0     0   8.996e-04   8.996e-04   1.597e-09
poisson       1000    20     5     0     5     0   1.480e-03   1.480e-03   7.341e-10
poisson       5000    20    13     0    13     0   8.234e-03   8.234e-03   5.293e-09
poisson      20000    20    18     0    18     0   4.525e-03   4.525e-03   3.681e-10
rescor         200    20     0     0     0     0   9.813e-04   9.813e-04   7.553e-09
rescor        1000    20     5     0     5     0   2.022e-03   2.022e-03   1.286e-09
rescor        5000    20    13     0    13     0   4.166e-03   4.166e-03   1.297e-09
rescor       20000    20    17     0    17     0   1.054e-02   1.054e-02   2.091e-09
sratio         200    20     0     0     0     0   7.642e-04   7.642e-04   1.237e-08
sratio        1000    20     1     0     1     0   1.205e-03   1.205e-03   5.258e-09
sratio        5000    20     8     0     8     0   3.374e-03   3.374e-03   3.382e-09
sratio       20000    20    16     0    16     0   8.922e-03   8.922e-03   9.448e-09

TOTAL                720   289     0   251     0

false-alarm rate on the reference build: 289 / 720 = 0.401
false-alarm rate on the lane build     : 0 / 720 = 0
the two stages scored here, for the mechanism: 251 after the bound projection, 0 after the headroom
largest headroom over all 720 correct fits: 2.105e-06 at nearzero n = 20000 seed 413459
  tolerance 0.001 , margin 475 x
headrooms that could not be computed: 0

     n   fits    med gmax   med |obj|    med head   warn rate
   200    180   1.895e-04   2.854e+02   1.997e-10       0.111
  1000    180   5.283e-04   1.422e+03   3.172e-10       0.222
  5000    180   1.125e-03   7.082e+03   1.901e-10       0.522
 20000    180   2.249e-03   2.840e+04   2.778e-10       0.750
```

Read the last block for the mechanism. The median absolute gradient grows
by 11.9x from 200 rows to 20,000, in step with the objective's 99.5x, and
the median headroom does not move at all (1.9e-10 to 2.8e-10). Both
readings are `nlminb` stopping on a relative change in a growing
objective. The warn rate follows the gradient.

The `bounded` column separates the two stages: the bound projection alone
takes 80 warnings to 0 at 200 rows, and the headroom is what removes the
rest, including the 20 of 20 at 20,000 rows where the projected gradient
is still 1.5e-02.

**The designs.** gaussian, poisson, bernoulli, cumulative, sratio, a
two-response gaussian with `set_rescor(TRUE)`, a `set_prior(ub = 0.1)`
box that binds, a two-parameter exponential-decay `nl = TRUE` formula,
and a mixed model whose variance component is truly zero. Sample sizes
200, 1000, 5000, 20000. 0 fit errors on either build.

## 5. Is the headroom real? Calibrated twice

`dev/gradcheck-02-truepos.R` and `dev/gradcheck-06-truepos2.R` compare
the predicted headroom against the log-likelihood shortfall of the same
fit against a properly converged one. The agreement is 3 to 4 digits on
GLM, ordinal and Laplace GLMM fits:

| construction | predicted | actual shortfall |
|---|---|---|
| gaussian `iter.max = 10` | 7.652e-07 | 7.65161e-07 |
| gaussian `iter.max = 5` | 1.586e-02 | 1.578e-02 |
| poisson GLMM `iter.max = 8` | 8.708e-04 | 8.6918e-04 |
| cumulative n=6000 `iter.max = 8` | 6.296e-02 | 6.27912e-02 |
| poisson `rel.tol = 1e-2` | 3.825e-01 | 3.82518e-01 |
| poisson GLMM `rel.tol = 1e-4` | 3.972e-03 | 3.94279e-03 |
| poisson GLMM `rel.tol = 1e-2` | 1.0428 | 1.04218 |

This is a MEASUREMENT and not an identity: the two sides are a quadratic
model of the objective and the objective itself, at two different points.

`dev/gradcheck-08-calib.R` (log `08-lane.txt`) repeats it by
re-optimizing from each stopping point with `rel.tol = 1e-14` and twelve
restarts. Nine rows agree to within 0.3 percent, two correct fits come in
at ratios 0.739 and 0.808 on headrooms of 2e-09, and ONE row disagreed by
seven orders of magnitude: a collinear design at `eps = 1e-6` predicted
0.528 and the harder optimization realized 2.3e-08.

`dev/gradcheck-09-newtonstep.R` settled that by taking the step and
evaluating the objective, which needs no optimizer:

```
collinear eps = 1e-05      H from he()       rcond 2.064e-11  |step|  4.259e+03
                           predicted     0.529474  best step drop     0.528352 at t = 1  ratio    0.9979
collinear eps = 1e-06      H from he()       rcond 2.064e-13  |step|  4.258e+04
                           predicted     0.529293  best step drop     0.528027 at t = 1  ratio    0.9976
collinear eps = 1e-07      H from he()       rcond 1.555e-15  |step|  5.927e+05
                           predicted     0.736638  best step drop     0.479873 at t = 0.5  ratio    0.6514
```

The headroom was right and the optimizer was wrong. `nlminb` with
`rel.tol = 1e-14`, `x.tol = 1e-14` and twelve restarts could not find
0.528 log-likelihood units that a single Newton step recovers, because
the step is 4.3e+04 long in a direction with reciprocal condition number
2.1e-13. That is a true positive the absolute gradient reported only as
5.1e-03, and it is the strongest argument for the criterion: the
quantity the check now reports is the one the user is missing.

## 6. The true positives, and the miss rate

Two sets, because they are not the same problem.

**Fits the optimizer's own code already condemns**
(`dev/gradcheck-02-truepos.R`, seeds 201 to 204; logs `02-base.txt`,
`02-lane.txt`). 17 rows, all with `convergence = 1`. Base raised the
gradient warning on 17, lane on 15. The two misses are the gaussian fit
at `iter.max = 10`, which is 7.65e-07 log-likelihood units short, and the
poisson GLMM at `iter.max = 8`, which is 8.69e-04 short. Both are under
the tolerance the user set, and both still warn, through "Optimizer did
not report convergence". **Miss rate for a fit warning at all: 0 of 17.**

**Fits `nlminb` calls converged** (`dev/gradcheck-06-truepos2.R`, seeds
401 to 407; logs `06-base.txt`, `06-lane.txt`). These are what a
gradient check is FOR, because nothing else catches them. 15 rows: 8 are
genuinely short, 4 are correct fits at large sample sizes, 1 is a
correct fit with a badly scaled column, and 2 (a flat nonlinear peak and
complete separation) have gradients of 1.5e-05 and 1.2e-257 and warn on
neither build.

| row | headroom | base | lane |
|---|---|---|---|
| poisson `rel.tol = 1e-2` | 3.825e-01 | warns | warns |
| poisson `rel.tol = 1e-3` | 1.323e-03 | warns | warns |
| poisson `rel.tol = 1e-4` | 1.323e-03 | warns | warns |
| poisson `rel.tol = 1e-5` | 2.256e-04 | warns | silent |
| poisson GLMM `rel.tol = 1e-2` | 1.0428 | warns | warns |
| poisson GLMM `rel.tol = 1e-3` | 2.940e-02 | warns | warns |
| poisson GLMM `rel.tol = 1e-4` | 3.972e-03 | warns | warns |
| collinear ridge `eps = 1e-5` | 1.003e-02 | warns | warns |
| badly scaled column (spread 1e6), correct | 6.95e-10 | warns | silent |
| gaussian n = 1e5, correct | 2.07e-11 | warns | silent |
| gaussian n = 4e5, correct | 1.13e-10 | warns | silent |
| `nearzero` n = 4e5, correct | 3.348e-05 | warns | silent |

**Miss rate on this set: 1 of 8**, the poisson at `rel.tol = 1e-5`, whose
shortfall is 2.257e-04 log-likelihood units. **False alarms removed: 4 of
4.**

Where the two classes are NOT separable, stated plainly. The largest
headroom on a correct fit anywhere in this lane's measurements is
3.348e-05 (`nearzero`, n = 400,000, seed 407). The smallest headroom on a
`convergence = 0` fit that is genuinely short is 2.256e-04. They are 6.7x
apart, so no threshold separates them; 1e-03 sits above both, which is
why one row is missed and none is a false alarm. Counting the
`convergence = 1` set too, the smallest short-fit headroom is 7.652e-7,
below the largest correct one, and that fit warns through the optimizer's
own code. Both borderline fits are inside 2.3e-04
log-likelihood units of their optimum, where no reported quantity moves.
The choice follows the standing rule that a check firing on a correct
model is worse than no check, and the margins are 30x above the worst
correct fit and 1.3x below the weakest true positive that is caught.

## 7. Found, not fixed

- **The trip-wire's gap, costing 1609.6 log-likelihood units.**
  `dev/gradcheck-07-gaps.R` (logs `07-base.txt`, `07-lane.txt`), seed
  501: `y ~ xs` with `sd(xs) = 1e-7` and `autoscale = FALSE` stops with
  `convergence = 0`, a maximum absolute gradient of 8.076e-05, and a log
  likelihood 1609.6 units below the well-scaled fit of the same data.
  NEITHER build warns, because 8.076e-05 is under the trip-wire. This is
  the same gap `autoscale_small_sd`'s comment records from the other
  side, and the `autoscale` default is what closes it in practice: a
  user has to set `autoscale = FALSE` to reach it. Closing it in the
  check would mean a Hessian on EVERY fit rather than only on the warning
  path, which is 4 percent of a whole fit at 4 parameters and 73 percent
  at 63 (`dev/gradcheck-04-cost.R`). Not paid here. The diagonally
  scaled gradient reads 40.0 on that fit against below 1e-2 on every
  correct one, so it is the cheapest known instrument for it if the
  package ever wants one.
- **The standard errors are not bound-aware.** On the seed 101 bounded
  fit the UNCONSTRAINED Hessian at the stopping point is indefinite
  (eigenvalues 693.40, 73.10, -29.63), so `diagnose()` reports
  `pdHess = FALSE` and NaN standard errors for `x` and
  `sigma_(Intercept)`. Restricted to the two parameters the bound does
  not hold, the same Hessian is positive definite, which is exactly why
  the headroom could be measured. A bound-aware covariance is its own
  change, now FILED in `dev/test-backlog.md` with this construction, the
  cap sweep and the test it wants.

  Whether it BITES depends on how far the bound sits from the
  unconstrained optimum, not on bound-awareness. The reviewer's sweep
  (`dev/gradcheck-rev-10-flip.R`) holds the bound active at every cap
  from 0.1 to 1.99 while `pdHess` is FALSE at 0.1 and 0.5 and TRUE from
  1 up, so the first spelling of the FLIP block, a bare
  `expect_false(pdHess)`, pinned this design's cap and not the defect:
  it fails in 6 of 8 caps where the defect is still there. Punch round 1
  rewrote it to assert the implication, which is what flips: the full
  Hessian here is indefinite (measured, with `skip_if_not()` rather than
  assumed), the Hessian over the parameters no bound holds is positive
  definite, and the report follows the full one.
- **A bound-active fit still burns a restart.** `optimize_obj()` restarts
  while the RAW maximum absolute gradient is above `grad_tol`, so the
  seed 101 fit restarts once for nothing. Making that loop
  bound-aware was deliberately NOT done: a restart that currently
  replaces `opt` moves the parameters by up to 3.3e-07 on that fit, so
  removing it would change fitted values in the tests that compare
  objectives to brms at 1e-12. Leaving the loop alone is what buys the
  bitwise identity in section 4, and the cost is one wasted optimizer run
  on a fit that was going to be silent.
- **A refit loop pays for a warning nobody sees: 1.91x at worst.**
  `frm_bootstrap()`, `influence()` and `frm_allfit()` wrap their refits in
  `suppressWarnings()`, which muffles the warning but does not stop it
  being BUILT, so every replicate that trips the trip-wire now pays for a
  Hessian.

  The first measurement of this was WITHDRAWN in punch round 1, and the
  number it reported was too low. It timed the two BUILDS in two
  processes, which the instrument rules forbid, and it offered an n = 400
  arm as a control that runs 0.13 s, 13 ticks of a 10.0 ms clock, and
  reads 0.857 in the effect arm AND in the control arm. That is the
  clock, not the code. `dev/gradcheck-12-refitcost.R` and its logs are
  kept as the record of the withdrawn measurement.

  Remeasured inside ONE process on ONE build
  (`dev/gradcheck-13-refitcost2.R`, log `13-lane.txt`): `restarts = 0` in
  every arm so `grad_tol` governs only the warning path, the "off" arm
  sets `grad_tol` to ten times the fit's own gradient so stage 2 never
  runs, arms interleaved per round, blocks past 1.2 s, minimum of 5
  rounds, and a SECOND "off" arm as the control.

  | design | stage 2 on | off | ratio | control |
  |---|---|---|---|---|
  | gaussian n=20000, `nsim = 20` | 2.660 s | 1.390 s | **1.914** | 1.022 |
  | cumulative n=20000, `nsim = 10` | 5.090 s | 3.480 s | 1.463 | 0.989 |
  | gaussian n=6000, tol raised | 0.447 s | 0.453 s | 0.985 | 0.993 |

  The third row is the control the first attempt did not have: both arms
  run identical code because `grad_tol` is above the fit's gradient in
  both, every block is over 0.44 s, and it reports 0.985 against a
  same-code control of 0.993. The reviewer's independent run of the same
  design got 1.833 with a control of 1.064
  (`dev/reviews/2026-09-29-gradcheck.md` section 8), so two instruments
  put the penalty between 1.8x and 1.9x on that design.

  Not fixed. The shape of the fix is an
  internal control field set by the five places that already do
  `ctl$verbose <- FALSE` for a refit loop, which would make those refits
  judge convergence by the absolute gradient again. That is a behavioural
  change to five call sites for a performance gain on a subset of
  designs, and `frm_allfit()` reads a gradient of its own, so it is not
  the one-line change it looks like.
- **A perturbed fit loses its box.** `fit_set_outer()` in
  `R/brms-shapes.R` replaces `fit$cache` with a fresh environment when it
  writes new estimates, which is right for a stored sdreport and drops the
  recorded box with it. `fit_outer_box()` then falls back to the caller's
  bounds alone, so a prior-spelled bound would not be excluded on such an
  object. Nothing routes a perturbed fit through `check_convergence()`,
  and the fallback errs toward warning, so this is recorded rather than
  fixed; recording the box outside the cache would mean a new slot on the
  fit object next to `bform`, which partial matching makes a worse
  trade.
- **`diagnose()$max_grad` ignores `par_units`.** `check_convergence()`
  judges the gradient in the natural units `autoscale` works in and
  `diagnose()` reports the raw one, so on an autoscaled fit the two
  numbers differ: 1001 against 1.685 on the reviewer's B2 fit
  (`dev/gradcheck-rev-03-probe.R`, seed 9202). `max_grad` is not changed,
  because seven extension test suites read it as a ratio to a log
  likelihood. Punch round 1 fixed the reporting of it: the claim that
  the difference was "documented on the manual page" was FALSE when it
  was written, so `?diagnose` now carries a TWO SCALES paragraph and the
  printed block says which lines are on which scale when the fit was
  standardized. The test that asserted
  `expect_identical(d$grad_proj, d$max_grad)` asserted the fixture; it
  now asserts `expect_null(f$par_units)` first, and a second block
  measures the autoscaled case where the two differ by exactly
  `par_units` at the argmax.

## 7b. Punch round 1, 2026-09-29

Review: `dev/reviews/2026-09-29-gradcheck.md`, scripts
`dev/gradcheck-rev-*`. Verdict was NOT MERGEABLE on four text defects,
with the criterion holding on every probe the reviewer built (0 lane-only
warnings over 180 fits of 9 new designs, 17 model shapes and 17 miss
rows; bitwise fits; the marginal Hessian confirmed; the 0.528 collinear
step confirmed to 15 digits against a hand-written gaussian likelihood
outside RTMB; the fuzz invariants unmuted with 0 firings). What changed.

1. **`?frmtmb_control` asserted a guarantee this lane's own measurement
   contradicts.** It ended "every fit stopped short of its optimum by
   more than `grad_tol` still warned". The counterexample is the lane's
   own section 7 defect: the reviewer reproduced it at their seed 501
   (`dev/gradcheck-rev-09-miss.R` row M1), `sd(xs) = 1e-7` with
   `autoscale = FALSE`, `convergence = 0`, max|grad| 1.6e-5, **297.175
   log-likelihood units short, silent on both builds**. The clause is
   gone and a WHAT IT DOES NOT COVER paragraph replaces it, naming the
   trip-wire's gate, that measurement, and the `autoscale` default as
   what closes it. The vignette carries the same paragraph. Verified by
   rendering (`dev/gradcheck-rd.R`, log `rd.txt`): 0 matches for the old
   clause, 0 stray percent signs.

2. **The warning named one quantity and printed another, and pointed at
   the wrong parameter.** `grad_warning_msg()` formatted `v$proj` inside
   a clause kept verbatim for grep-compatibility. On the reviewer's seed
   9301 (`dev/gradcheck-rev-16-msg.R`) it read "Large maximum absolute
   gradient at the optimum (52.5)" where the maximum was 5971.25, then
   said "diagnose() names the offending parameter" and `diagnose()` named
   `x`, the one parameter the verdict had excluded.

   The prefix "Large maximum absolute gradient" STAYS, because
   `tests/testthat/helper-fuzz.R:1763` greps it in
   `FUZZ_NONCONVERGENCE` to decide whether an invariant was measured on a
   fit that warned. Checked before touching it: that line, `R/fit.R` and
   `test-grad-verdict.R` are the only matches in `R/`, `tests/` and every
   extension's tests, sources and vignettes.

   The NUMBER after the prefix is now the maximum absolute gradient,
   which is what the prefix says, with its parameter; the projected
   gradient, ITS parameter and the headroom follow in their own clauses;
   and the pointer is "diagnose() reports all three numbers". Same fit
   after (`dev/gradcheck-log/rev16-lane-after.txt`):

   ```
   Large maximum absolute gradient at the optimum (5971 at x), of which
   1 component is held by a bound (x), leaving 52.5 at z as the largest
   no bound holds: one Newton step over the parameters no bound holds
   would still gain 0.131 in log-likelihood, more than
   frmtmb_control(grad_tol = 0.001). ...

   Max |gradient|: 5971 at x
     1 parameter is held by a bound (x); the largest gradient no bound
     holds is 52.49 at z
   ```

   `grad_verdict()` gains `proj_par` and `diagnose()` gains
   `grad_proj_par` for it.

   The grep compatibility is ASSERTED, not assumed: a passing gated fuzz
   run proves nothing here, because it passes either way when nothing
   fires. `dev/gradcheck-15-fuzzgrep.R` reads `FUZZ_NONCONVERGENCE` out
   of `helper-fuzz.R` and matches it against all five branches of the
   rewritten message, including the one where no parameter name is
   available: 5 of 5 (log `15-lane.txt`). The gated tier itself is green
   as well, `FRMTMB_FUZZ=true FRMTMB_FUZZ_N=300`, pass=2 fail=0.

3. **`vignettes/diagnostics.Rmd` cited a 400-row cell that does not
   exist.** The grid is 200 / 1000 / 5000 / 20000
   (`dev/gradcheck-03-design.R`); the 0-of-20 figure is the 200-row cell.
   One character.

4. **`?diagnose` omitted the units caveat section 7 claimed it carried.**
   It now has a TWO SCALES paragraph: `max_grad` and `worst_grad` are
   raw, `grad_proj`, `grad_proj_par` and the warning are in `par_units`,
   and `grad_headroom` is invariant to the scaling because a Newton
   decrement is. The printed block says so on a standardized fit, verified
   on the reviewer's own B2 probe
   (`dev/gradcheck-log/rev03-lane-after.txt`).

5. **`diagnose()` published a headroom on an importance-corrected fit
   that `check_convergence()` refuses to judge.** `grad_verdict()` is now
   skipped when `fit$importance` is present, the three fields are `NA`,
   the printed block says the gradient is a Monte Carlo estimate no
   criterion applies to, and the clean line is withheld rather than
   inherited from the absence of a verdict, which is the one way this
   change could have become LESS diagnostic than the `max_grad < 1e-3`
   it replaced. A test block asserts all of it.

6. **The FLIP block pinned a cap, not bound-awareness.** See the
   standard-errors bullet in section 7 and the new
   `dev/test-backlog.md` entry.

7. **The refit-loop cost was remeasured and the old control withdrawn.**
   See the refit-loop bullet in section 7. 1.914 with a control of 1.022,
   against 1.65 with no usable control.

8. **The one absolute tolerance in the test file is gone.**
   `abs(headroom/realized - 1) < 0.05` is replaced by a yardstick the run
   measures. Along the Newton direction an exactly quadratic objective
   drops `D * (2t - t^2)`, so the half-step drop is exactly `0.75 * D`;
   the departure from 0.75 is this objective's non-quadraticity at that
   point, and the headroom's own error is that same third-order term.
   Measured over six constructions (`dev/gradcheck-14-quadyard.R`, log
   `14-lane.txt`) the ratio of the two is 0.959, 1.200, 1.200, 1.288 and
   1.199 on unbounded rows, so the assertion is `error <= 10 * yardstick`
   with a floor at `sqrt(.Machine$double.eps)`. The bounded row reads
   6.03 and is NOT used: there the full-Hessian Newton step walks off the
   bound, so it is not the step the criterion takes. The two hardcoded
   `0.1` values now read `r$cap`.

9. **`diagnose()`'s own Hessian cost is recorded**, in section 9 and on
   `?diagnose`: free on a fit that already warned, because the verdict is
   cached, and 0.13 to 0.45 s on the large latent designs when the cache
   was replaced.

10. **`fit_outer_box()`'s comment described a bound no caller can pass.**
    `frm()` has no `lower`/`upper` argument, so the `resolve_bounds()`
    fallback can only return the open box the line above it returns. The
    comment now says that, says `set_prior(lb =, ub =)` is the only route,
    and says why an unreachable open-box fallback is still the right
    answer for an object assembled elsewhere.

## 7c. Punch round 2, 2026-09-29: the importance guard failed CLOSED

One blocking item, and it was MINE, introduced by punch item 5 rather
than found in the original change. Everything else from the review is
verified fixed in its "Re-check" section.

11. **The clean line was withheld from every importance-corrected fit.**
    `diagnose()` sets `gv <- NULL` on such a fit, and the clean-line
    condition read
    `(degenerate || (!is.null(gv) && !isTRUE(gv$warn)))`, whose second
    conjunct is FALSE for every one of them. So no importance fit could
    print "No convergence problems detected", however well it had
    converged. 0.64.0 gated the same line on `out$max_grad < 1e-3`, which
    a well-converged importance fit passes, so this was a REGRESSION
    against the build it replaces, in the opposite direction from the one
    I had argued for: I reasoned that withholding is the conservative
    side, and the conservative side is not the same as the correct one.
    Reporting a problem on a fit that has none is the defect this whole
    lane exists to remove.

    The reviewer's construction, `dev/gradcheck-rev-20-impclean.R`, is a
    20-group poisson GLMM with `importance = 2000L` and `se = TRUE`. It
    also carries the ABSENT case, two draws whose Monte Carlo gradient
    lands ABOVE `grad_tol`, so a fix that printed the line
    unconditionally would be caught. On the fixed build
    (`dev/gradcheck-log/rev20-lane-after.txt`):

    | seed | max\|grad\| | under `grad_tol` | clean line |
    |---|---|---|---|
    | 9503 | 0.00094754 | TRUE | TRUE |
    | 9501 | 0.0034419 | FALSE | FALSE |
    | 9502 | 0.013354 | FALSE | FALSE |

    All three have `convergence = 0`, `pdHess = TRUE`, no bad standard
    errors, and every other clean conjunct satisfied, so the gradient is
    the only thing that can decide. The same three designs fitted WITHOUT
    `importance` land at 7.9e-05, 8.1e-05 and 5.3e-05 and print the line,
    which is the reference for what an otherwise identical fit does.

    The rule the code now states is "no verdict means the old gate
    decides": with `gv` NULL the condition falls back to
    `out$max_grad < control$grad_tol`, the 0.64.0 gate with the setting
    honoured instead of hardcoded. The `is.null(gv) ||` spelling was
    rejected because seeds 9501 and 9502 sit above the tolerance and must
    still withhold the line.

    Two test blocks, one per side. Seen to fail on the pre-fix lane
    library (`dev/gradcheck-log/10-lane-prefix.txt`): **pass=71 fail=1**,
    `an importance fit under grad_tol still prints the clean line`,
    `Expected any(grepl("No convergence problems detected", r$out)) to be
    TRUE. actual: FALSE`. Both blocks PASS on the reference build, which
    is what makes this a regression rather than a new feature, and the
    over-tolerance block passed on the pre-fix library too, so it is a
    genuine complement and not a copy of the first.

## 8. Decided against

- **Changing the restart criterion** or `quad_fit()`'s `stationary` flag.
  Both would move fitted values. See section 7.
- **Touching the importance branch.** Out of scope by the brief, and its
  reasoning is separate: an importance-corrected objective has a Monte
  Carlo gradient with an O(N^-1/2) error, so no gradient criterion of any
  shape applies to it.
- **Adding a second control argument** for the headroom tolerance. One
  knob that moves both readings together is what the brief asks for and
  what keeps the documentation honest; a second would have to be
  explained in terms of the first anyway.
- **Using `sdr$cov.fixed` as the Hessian when `se = TRUE`** (free, since
  sdreport already ran). Rejected: `cov.fixed` is the inverse of the FULL
  Hessian, so its free-set block is not the inverse of the free-set
  Hessian, and under `autoscale` it has already been mapped back to
  another parameterization. Two unit traps for a saving on a path that
  only a warning reaches.

## 9. Cost

`dev/gradcheck-04-cost.R`, log `04-base.txt`. Blocks grown past 1.5 s,
minimum of 3 rounds, per call.

| fit | np | one gradient | `he()` | `optimHess` | whole fit | Hessian / fit |
|---|---|---|---|---|---|---|
| gaussian GLM n=500 | 4 | 3e-5 s | 1.4e-4 s | 2.8e-4 s | 0.007 s | 0.042 |
| cumulative n=20000 | 5 | 5.6e-3 s | 4.7e-2 s | 8.3e-2 s | 0.413 s | 0.202 |
| poisson GLMM n=2000 q=50 | 3 | 4.2e-4 s | none | 8.9e-3 s | 0.065 s | 0.139 |
| gaussian LMM | 63 | 9.7e-4 s | none | 0.191 s | 0.263 s | 0.726 |

`he()` is about half the cost of `optimHess` where RTMB offers it, which
is every model without random effects; a Laplace objective answers
"Hessian not yet implemented" and falls back. Only a fit whose projected
gradient already tripped `grad_tol` pays anything, so a healthy fit pays
nothing, and a fit that trips was already paying a restart.

**Three places pay, not one.** The table above is the FIT path, which is
the only one the first pass measured.

- The fit path, once per fit that trips, the table above.
- `diagnose()`, which did not build a Hessian before. Free on a fit that
  already warned, because the verdict is cached on `fit$cache`, and 0.13
  to 0.45 s on the large latent-variable designs when the cache was
  replaced, which is what `fit_set_outer()` does to a perturbed fit
  (measured by the reviewer, `dev/reviews/2026-09-29-gradcheck.md`
  sections 1(b) and 9.3). Recorded on `?diagnose`.
- A refit loop, where `suppressWarnings()` muffles the warning but does
  not stop it being built: 1.914x on a gaussian bootstrap at n = 20000.
  See the refit-loop bullet in section 7.

## 10. Tests

`tests/testthat/test-grad-verdict.R`, 18 blocks, 72 expectations: 12
blocks and 46 at the end of the first pass, 16 and 63 after punch round
1, 18 and 72 after punch round 2.

Seen to fail on the reference build (`dev/gradcheck-10-seenfail.R`, log
`10-base.txt`): **pass=38 fail=18 error=8** on the current file, and
**pass=23 fail=10 error=7** on the file as it stood at the end of the
first pass, which is the figure the reviewer reproduced. Punch round 2's
two blocks PASS on base, because the thing they pin was a regression
against base rather than a defect in it; they were seen to fail on the
pre-fix LANE library instead (section 7c). The behavioural
failures on base, not the missing-symbol ones, are:

- "a parameter held by a bound does not raise the gradient warning":
  `Expected r$grad to be FALSE. actual: TRUE`
- "a gaussian GLM at n = 50000 trips the trip-wire and is clean": same
- "an ordinal fit at n = 6000 agrees with polr and stays clean": same
- "a start far away warns even when the curvature is unusable":
  the warning does not say what could not be measured
- "a flat ridge warns on curvature the gradient does not show": the
  warning carries no log-likelihood number
- "an importance-corrected fit gets no verdict at all": base publishes a
  `grad_proj` and a `grad_headroom` there and prints no explanation
- "diagnose() names the bound that holds the gradient": base prints
  neither "held by a bound" nor the parameter the verdict is about

On the lane build: **pass=72 fail=0 error=0 skip=0**
(`dev/gradcheck-log/punch2/test-grad-verdict.R.txt`), against pass=63 at
the end of punch round 1 (`10-lane.txt`).

The blocks that must keep firing are there too: a loosened optimizer on a
GLM and a GLMM (both with `convergence = 0`, so the gradient check is the
only thing standing between the user and a short fit), a fit at its
iteration cap, a far start whose Hessian is unusable, and the collinear
ridge. `grad_bound_active()` is also unit-tested with the guarded thing
ABSENT: an inward-pointing gradient at a bound, an open box, and a
parameter well inside its bound all return FALSE.

No absolute numeric tolerance appears in the file. Every assertion is a
ratio to the fit's own `grad_tol`, to its own standard errors, to the
bound it was given, to this objective's own departure from a quadratic,
to `.Machine$double.eps`, or to a quantity the run measures twice by two
routes; every bare number left is a MARGIN on one of those ratios. The
one exception the reviewer found, `abs(headroom / realized - 1) < 0.05`,
is gone (punch item 8).

### The core suite

All 180 files in `tests/testthat/`, one per R process, eight at a time
(`dev/gradcheck-suite.ps1`, counted by `dev/gradcheck-suite-count.sh`,
logs under `dev/gradcheck-log/suite-lane/`):

Run twice, once at the end of the first pass and again on the build punch
round 1 produced:

| run | files | pass | fail | error | skip |
|---|---|---|---|---|---|
| first pass | 180 | 11924 | 0 | 0 | 152 |
| after punch round 1 | 180 | 11941 | 0 | 0 | 152 |

The 17 extra passes are the four new blocks in `test-grad-verdict.R`. The
second run's per-file record is
`dev/gradcheck-log/suite-lane-punch1-flat.txt`, the driver's is
`suite-lane-punch1.txt`, and 0 of 180 files came back without a result.

The 152 skips are the gated tiers, which need
`FRMTMB_BRMS_FIT_TESTS=true` (`test-brms-likelihood.R` 37,
`test-brms-methods.R` 46, `test-brms-priors.R` 12,
`test-drmtmb-agreement.R` 13, thirteen `test-bcm-*.R` files 39,
`test-brms-agreement.R` 2, `test-rl-example.R` 2, `test-fuzz.R` 1). Every
file with a skip is named in `dev/gradcheck-log/suite-lane-count.txt`.
`test-importance.R` lost its log to the driver and was rerun on its own
(`dev/gradcheck-log/suite-lane-importance-rerun.txt`, pass=226); its 226
are in the total.

One edit landed after that run: `grad_headroom()` puts the tape back at
the optimum after `optimHess()`, because `optimHess()` leaves the last
evaluation at a perturbed point and on a Laplace objective that is a
displaced inner solve (`diagnose_flat()` has the same pattern and does
not restore, and `fit_end_checks()` runs after the convergence check).
The eight files that can reach a convergence warning were rerun against
the reinstalled build (`dev/gradcheck-log/post-edit/`):

```
test-autoscale.R                           pass=46 fail=0 error=0 skip=0
test-conditions.R                          pass=150 fail=0 error=0 skip=0
test-confint-anova.R                       pass=35 fail=0 error=0 skip=0
test-diagnostics-ux.R                      pass=127 fail=0 error=0 skip=0
test-grad-verdict.R                        pass=46 fail=0 error=0 skip=0
test-predfix.R                             pass=92 fail=0 error=0 skip=0
test-thres.R                               pass=67 fail=0 error=0 skip=0
test-verbose.R                             pass=40 fail=0 error=0 skip=0
```

Punch round 1 reran the same set plus the four files the reviewer added
to it, before the whole-suite rerun above
(`dev/gradcheck-log/punch1/`): PASS 953, FAIL 0, ERROR 0, SKIP 0 over
`test-autoscale.R` 46, `test-bracket-access.R` 33, `test-conditions.R`
150, `test-confint-anova.R` 35, `test-diagnostics-ux.R` 127,
`test-grad-verdict.R` 63, `test-importance.R` 226, `test-predfix.R` 92,
`test-review-v28.R` 74, `test-thres.R` 67 and `test-verbose.R` 40.

Punch round 2 reran the four files its one-line change can reach
(`dev/gradcheck-log/punch2/`): `test-grad-verdict.R` pass=72,
`test-importance.R` pass=226, `test-diagnostics-ux.R` pass=127,
`test-confint-anova.R` pass=35, all with fail=0 error=0 skip=0, 460
passing in total. The last two are the files that assert the clean line
on an ordinary fit, which is the line the regression touched.

The GATED fuzz tier matters here more than anywhere else, because
`helper-fuzz.R` greps the warning's opening phrase to decide which
invariants a warning mutes, and punch item 2 rewrote that warning.
`FRMTMB_FUZZ=true FRMTMB_FUZZ_N=300`: pass=2 fail=0 error=0 skip=0
(`dev/gradcheck-log/punch1/test-fuzz-gated.txt`), matching the
reviewer's own gated run.

## 11. R CMD check

One run, `--as-cran`, on a tarball built WITH vignettes so the rebuilt
`diagnostics` article is checked (`dev/gradcheck-log/check-build3.txt`,
`dev/gradcheck-log/check3.txt`, and the check log kept as
`dev/gradcheck-log/check3-00check.log`; the .Rcheck directory and the
tarball were removed afterwards).
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, and
`_R_CHECK_FORCE_SUGGESTS_=false` because `drmTMB` is not in the lane
library and installing it there is out of the question.

An earlier attempt was discarded rather than counted: it was built with
`--no-build-vignettes`, which makes the check report a vignette-index NOTE
that is an artifact of the build, and it stopped at
"Package suggested but not available: 'drmTMB'" before checking anything.

Every R source, manual page, `NEWS.md` and the vignette predate the
tarball (23:15:32 against 22:38 to 23:10). `test-grad-verdict.R` was
edited afterwards (23:20:46), loosening two assertions: an exact `0.1` on
the bound became a relative comparison against
`sqrt(.Machine$double.eps)`, and a `polr` log-likelihood gap was compared
to `grad_tol` instead of to the headroom, which is a ratio with a 1e5
margin rather than a 1.1x one. Both versions read pass=46 fail=0 error=0.

**Status: 1 NOTE**, and it is the expected one:

```
* checking examples ... [43s] OK
* checking examples with --run-donttest ... [46s] OK
* checking tests ... [12m] OK
  Running 'testthat.R' [12m]
* checking re-building of vignette outputs ... [179s] OK
* checking PDF version of manual ... OK
* checking HTML version of manual ... [20s] NOTE
Skipping checking math rendering: package 'V8' unavailable
Status: 1 NOTE
```

No WARNING, no ERROR. The V8 math-rendering NOTE is the one
`dev/lane-rules.md` names as expected on this box.

NOT rerun after punch round 1, because the rule is one run per lane and
the consolidating session runs the authoritative one. What changed since
it: the warning's wording, one new `diagnose()` field, the importance
guard, the printed lines, two manual pages, the vignette prose and the
test file. The hazard a check catches in that set is an Rd rendering
defect, and both pages were rendered with `tools::Rd2txt`
(`dev/gradcheck-rd.R`, log `rd.txt`): 261 and 129 lines, 0 stray percent
signs, every new paragraph present. The test and vignette changes are
covered by the whole-suite rerun above.

## 12. Scripts and logs

| script | what |
|---|---|
| `dev/gradcheck-helpers.R` | scoring, shared by every arm |
| `dev/gradcheck-01-construct.R` | the false alarms, with their evidence |
| `dev/gradcheck-02-truepos.R` | true positives with `convergence = 1` |
| `dev/gradcheck-03-design.R` | 720 replicates, three process groups |
| `dev/gradcheck-04-cost.R` | gradient, `he()` and `optimHess` timings |
| `dev/gradcheck-05-summarise.R` | the counts in section 4, generated |
| `dev/gradcheck-06-truepos2.R` | true positives with `convergence = 0` |
| `dev/gradcheck-07-gaps.R` | the two gaps in section 7 |
| `dev/gradcheck-08-calib.R` | headroom against a harder optimization |
| `dev/gradcheck-09-newtonstep.R` | headroom against the step itself |
| `dev/gradcheck-10-seenfail.R` | the new tests failing on the base build |
| `dev/gradcheck-11-messages.R` | the warning text and diagnose block |
| `dev/gradcheck-12-refitcost.R` | WITHDRAWN: refit cost, two processes |
| `dev/gradcheck-13-refitcost2.R` | refit cost, two readings in one process |
| `dev/gradcheck-14-quadyard.R` | the measured yardstick the test uses |
| `dev/gradcheck-15-fuzzgrep.R` | the fuzz harness still greps the warning |
| `dev/gradcheck-backlog.R` | the dev/test-backlog.md entry |
| `dev/gradcheck-news.R` | the NEWS edit |
| `dev/gradcheck-runtests.R` | one test file per process, with counts |

Logs are in `dev/gradcheck-log/`, named after the script and the arm.

## 13. Version

The default warning changes on about 40 percent of large-sample fits, and
`diagnose()` gains three return fields and two printed lines. No estimate
and no exported signature changes. That is a MINOR bump, not a patch: the
number is not chosen here.
