# Item 2.5: the Royston-Parmar frailty against rstpm2

Lane `frailty`, worktree `frmtmb-wt-frailty`, branch `wt-frailty`,
based on the round 3 release (frmtmb 0.55.2, frmtmb.spline 0.5.0).

The row, verbatim from `dev/extension-gaps-plan.md`:

> spline RP: `(1 | centre)` frailty against rstpm2's log-normal frailty
> on the same data; a random effect on `gamma1`. Reference: rstpm2.
> Current state: the current frailty test is a smoke test for
> finiteness.

**The claim holds.** frmtmb's `(1 | centre)` on `mu` and rstpm2's
`stpm2(cluster =, RandDist = "LogN")` are the same model, and the
estimates agree to a fraction of a thousandth of a standard error over
200 replicates at the plan's realistic scale. The two report DIFFERENT
log likelihoods, and that difference is entirely the way each one
integrates the frailty out. A random effect on `gamma1` also works,
recovers, and carries one property a user has to be told about: it is
anchored at `t = 1`.

## What the two packages are fitting, side by side

This was settled before anything was compared, because the two write
the same curve in different coordinates.

| | frmtmb.spline | rstpm2 |
|---|---|---|
| baseline | `eta = gamma0 + gamma1 x + sum_j gamma_{j+1} v_j(x)`, `x = log t`, `v_j` the Royston-Parmar natural cubic basis | `eta = a0 + sum_j c_j nsx_j(x)`, `nsx` a natural spline with the same knots |
| spline space | natural cubic splines on the knots, dimension `df + 1` | the same space, dimension `df + 1` |
| covariates | any dpar takes a formula; `mu` is `gamma0` so it is proportional hazards, `gamma1` upward is time varying | `formula` is proportional hazards, `tvc` is time varying |
| frailty | `(1 \| centre)` on `mu`, a random-effect block with `sd` on the log scale | `cluster =`, `RandDist = "LogN"`, reported as `logtheta` |
| what theta means | `VarCorr()` holds the covariance; the printed number is the standard deviation | theta is the VARIANCE, so `sqrt(exp(logtheta))` is frmtmb's `sd` |
| how the frailty is integrated out | Laplace approximation (RTMB `random =`) | adaptive Gauss-Hermite, `nodes = 9` by default |
| what `logLik()` reports | the Laplace value | the quadrature value |

The two spline spaces really are the same space, and that is measured
rather than asserted: rstpm2's fitted `eta` on a 400-point grid is
regressed on the Royston-Parmar basis at the same knots, and the
residual is at most **2.96e-14 relative to the linear predictor's own
scale** over the 200 replicates. Move one interior knot by 0.6 and the
same residual is 1.59e-02, so the check is not vacuous
(`frailty-falsify.R`).

Knots are pinned on both sides in every comparison. `royston_parmar()`'s
default rule and the pinned values agree exactly: `max abs diff 0` at
seed 20260910 (`frailty-recover.R`, its preamble).

## The recovery table

`frailty-recover.R`, seeds 20260910 to 20261109, 200 replicates. Design:
2000 subjects, 40 centres of 50, 40 percent censored, `df = 3`. The data
come from a Weibull proportional-hazards model with a shared log-normal
frailty, which IS a Royston-Parmar model whose spline is linear in log
time, so every truth is known at any `df`: a globally linear function
has one representation in the RP basis and its interior coefficients are
zero.

| quantity | truth | mean | sd | coverage |
|---|---|---|---|---|
| `beta(trt)` | 0.6 | 0.5984 | 0.0587 | 186/200, 0.930 |
| `sd(frailty)` | 0.5 | 0.4905 | 0.0689 | 184/200, 0.920 |
| `gamma0` | -2.0923 | -2.093 | 0.386 | 192/200, 0.960 |
| `gamma1` | 1.3 | 1.312 | 0.188 | 190/200, 0.950 |
| `gamma2` | 0 | 0.0069 | 0.0664 | 190/200, 0.950 |
| `gamma3` | 0 | -0.0146 | 0.129 | 191/200, 0.955 |

Binomial Monte Carlo error at 0.95 with 200 replicates is 0.0154, so
this run can resolve a coverage shortfall of about 0.031 at two Monte
Carlo errors. `beta` at 0.930 is 1.3 of those below nominal and
`sd(frailty)` at 0.920 is 1.9 below: neither is a finding at this
replicate count, and the pair is what a reader should expect to see
again rather than a defect to chase.

**The punch round settled that at 800 replicates and the reading holds**
(`dev/reviews/20260910-frailty.md`, section 4, whose first 200 seeds are
this row's own, every figure reproducing to the digit). Over 800, `beta`
covers 0.9450 at a Monte Carlo error of 0.0081 and `sd` covers 0.9425,
z of -0.65 and -0.97 against nominal; the 600 seeds this row never ran
give exactly 0.9500 each. There is no under-coverage here. As a
by-product the same run gives mean `se(beta)` 0.05956 against an
empirical `sd(beta)` of 0.05991, a ratio of 0.9942, so the standard
errors are calibrated too.

What IS real is the LEVEL of `sd(frailty)`, and 200 replicates
understated it. 0.4905 against 0.5 is 1.95 Monte Carlo standard errors
here; over 800 it is 0.4889 at a standard error of 0.0023, which is
**4.8**. It is not the integration rule: rstpm2's quadrature estimate is
0.1386 percent higher on 200 of 200 replicates, a fourteenth of the
shortfall. It is the ordinary downward bias of a maximum likelihood
variance component, which for `m = 40` groups and `p = 2` fixed effects
predicts `0.5 * sqrt(38 / 40) = 0.4873` against the observed 0.4889 and
leaves 0.3 percent over. The misses are asymmetric in the same
direction: of the 46 `sd` intervals that missed, 36 lay entirely below
0.5 and 10 above.

Convergence: `nlminb` code 0 on 199 of 200, a positive definite Hessian
on 200 of 200, and a largest maximum absolute gradient of 0.049 (seed
20260995). Four replicates were flagged on code or gradient
(20260938, 20260940, 20260995, 20260997), and their agreement with
rstpm2 is no worse than the rest: 2.5e-04 to 3.6e-04 of a standard error
on `beta`, and 1.2e-03 to 1.4e-03 relative on `sd`.

Floors: 0 non-monotone rows over 200 fits, and the deepest censored
`-log S` is 5.79 against `rp_floored()`'s threshold of 19.2. Nothing on
this design is a floor artifact.

## The agreement with rstpm2

Same 200 replicates, same data, matched knots.

| quantity | mean | worst |
|---|---|---|
| `beta(trt)`, absolute | 1.93e-05 | 5.27e-05 |
| `beta(trt)`, signed: rstpm2 higher | on 200 of 200 | |
| `beta(trt)`, in units of the run's own standard error | 3.23e-04 | 8.82e-04 |
| `sd(frailty)`, relative | 1.38e-03 | 1.96e-03 |
| `sd(frailty)`, signed: rstpm2 higher | on 200 of 200 | |
| `se(beta)`, relative | 2.56e-06 | 1.31e-05 |
| basis map residual, relative | | 2.96e-14 |

rstpm2 fitted every one of the 200 without failing.

## The log likelihoods differ, and it is the rule and not the model

This is the part a reader is most likely to mistake for a defect.
frmtmb reported -2645.410 and rstpm2 -2645.346 on the first dataset
tried, a gap of 0.064, on a model both agreed about to five decimal
places in every coefficient.

The gap was attributed by writing a THIRD log likelihood that belongs to
neither package: the marginal likelihood integrated per cluster by
`stats::integrate()` at `rel.tol = 1e-12`, in plain R, from the basis
written out by hand (`frailty-common.R`, `frailty_exact_ll()`). Scored
at each package's own optimum, over the 200 replicates:

| rule | value minus the exact integral, mean | range |
|---|---|---|
| frmtmb, Laplace | -0.0577 | -0.0866 to -0.0155 |
| rstpm2, 9-node adaptive Gauss-Hermite | -4.90e-06 | -1.3e-05 to 0 |

So the whole reported gap is the integration rule. The exact integral
also says what that costs the ESTIMATE, which is the question that
matters: the exact log likelihood at frmtmb's estimate is **5.53e-05
units below** the exact log likelihood at rstpm2's on average, worst
7.7e-05, which is **0.105 percent of the Laplace offset itself** (worst
0.457 percent). That magnitude is the informative number, and at that
size it is invisible next to a standard error of 0.0587 on `beta`.

frmtmb's point was the better one under the exact criterion on 0 of 200
replicates, and that unanimity is STRUCTURAL rather than 200
independent confirmations: rstpm2's rule is accurate to 4.9e-06, so
rstpm2's estimate sits essentially at the exact maximizer, and any
displaced point must then score lower. What 200 replicates buy is the
distribution of the magnitude, not the direction.

rstpm2's own rule was checked against the exact integral at 3, 9, 21 and
51 nodes on one dataset (`frailty-ident.R`, seed 20260910): -0.0800,
-5.7e-06, +7.5e-09, +7.5e-09. Its 9-node default is converged for this
design and its 3-node setting is worse than the Laplace approximation.

### Where the Laplace approximation stops being free

`frailty-sweep.R`, seeds 20260910 to 20260914, 5 replicates per cell,
2000 rows throughout, `df = 2`. Cluster size is the axis: the Laplace
approximation is exact for a Gaussian conditional and a cluster's
conditional log likelihood approaches one as its event count grows.

| centres | per centre | events per centre | Laplace minus exact | GH minus exact | exact ll gap | `beta` in se | `sd` relative |
|---|---|---|---|---|---|---|---|
| 10 | 200 | 120 | -0.0066 | -7.6e-08 | -1.3e-07 | 1.8e-05 | 1.2e-04 |
| 40 | 50 | 30 | -0.0666 | -6.2e-06 | -4.8e-05 | 3.3e-04 | 1.3e-03 |
| 200 | 10 | 6 | -0.0875, range -0.0105 to -0.1923 | -1.5e-04 | -6.7e-03 | 9.7e-03 | 9.7e-03 |
| 500 | 4 | 2.4 | +1.38 | -3.6e-04 | -3.3e-03 | 7.5e-03 | 5.2e-03 |

The 200-centre cell is quoted with its range because it is the one that
needs it: five seeds span a factor of 18, -0.0105 to -0.1923, since that
cell straddles the sign change and per-cluster errors of both signs
partly cancel. The other three cells are tight, -0.00685 to -0.00641,
-0.0724 to -0.0613 and +1.337 to +1.411.

At four subjects per centre the Laplace marginal likelihood is 1.38
units ABOVE the exact one, and the sign has flipped: the offset is not a
bound. The estimates still agree, to 0.0075 of a standard error on
`beta` and 0.5 percent on `sd`. Two independent rules landing 3.6e-04
apart at 500 clusters is also the evidence that the exact reference is
right, since a broken reference would have to break both.

Frailty size matters much less than cluster size. At 40 centres of 50,
the Laplace offset is -0.0105 at `sd = 0.25`, -0.0666 at 0.5 and -0.118
at 1.0, while `beta` agrees to 3.8e-04 of a standard error at every one.

Practical reading, which is what went into `?royston_parmar`: **compare
`logLik()` across frmtmb fits, not across packages**, and treat the
frailty variance as the parameter the offset touches first.

## A random effect on gamma1: what it means

`gamma1` multiplies `x = log t`. A centre deviation `u` therefore gives
that centre a cumulative hazard `t^(gamma1 + u)` and a hazard ratio
against another centre of `t^(u - u')`, which MOVES WITH TIME. With
`df = 1` it is exactly a per-centre Weibull shape. It is not a frailty:
a frailty scales the cumulative hazard by a constant, and this tilts it.
rstpm2 has no counterpart, so this half has no third-party reference and
is measured against simulation and against the same exact integral,
extended to a random slope (`frailty_exact_ll_slope()`).

### It recovers

`frailty-gamma1-recover.R`, seeds 20260910 to 20260969, 60 replicates,
2000 subjects in 40 centres, `df = 1`, true `sd(gamma1 | centre)` 0.2
against `gamma1 = 1.3`, 40 percent censored.

| quantity | truth | mean | sd | coverage |
|---|---|---|---|---|
| `sd(gamma1 \| centre)` | 0.2 | 0.1994 | 0.0368 | 54/60, 0.900 |
| `gamma1` | 1.3 | 1.300 | 0.0449 | 57/60, 0.950 |
| `gamma0` | -2.0923 | -2.089 | 0.0625 | 56/60, 0.933 |
| `beta(trt)` | 0.6 | 0.599 | 0.0614 | 55/60, 0.917 |

Binomial Monte Carlo error at 0.95 with 60 replicates is 0.0281, so this
run resolves a shortfall of about 0.056. The 0.900 on `sd` is 1.8 of
those below nominal and is not a finding at 60 replicates; saying more
about it needs more replicates than this row is sized for.

### And it earns its place

A component is not validated by being estimable. The per-centre Weibull
shape is what it exists to estimate, so the error on that is the test,
against two nulls:

| arm | mean per-centre `abs(shape_hat - shape_true)` | beats it on |
|---|---|---|
| the block ON | 0.0892 | |
| a pooled fit with no block at all | 0.1629 | 60 of 60 |
| the SAME fit with its own deviations dropped | 0.1626 | 60 of 60 |

The second null is the one that matters, and it caught a weak assertion.
At the test-file size, the ON fit's population `gamma1` alone already
beats the pooled fit, 0.2336 against 0.2359, so `err_on < err_off` would
pass for a block that contributed nothing (`frailty-falsify.R`,
"an inert gamma1 block still beats the pooled fit": broken TRUE). The
shipped test compares against the fit's own population value instead.

Correlation between the fitted and the true per-centre shape: mean
0.838, minimum 0.676. Log likelihood with the block minus without:
mean 23.0, minimum 3.68 over 60 replicates, for two extra parameters.

Convergence: code 0 and a positive definite Hessian on 60 of 60, largest
maximum absolute gradient 1.3e-03.

Here the Laplace offset has the OTHER sign: +0.087 on average, range
+0.034 to +0.131. A random slope is a different integrand from a random
intercept and the sign of the offset is a property of the integrand, not
of the package.

### The monotonicity floor, which did not fire, and why

`gamma1 + u` has to stay positive or that centre has no hazard, and the
identity link does not hold it there. Over the 60 replicates:
0 non-monotone rows, smallest fitted per-centre slope 0.610, and the
smallest margin to zero was 2.61 estimated standard deviations (mean
4.94).

A design built to reach the floor did not reach it either
(`frailty-floor.R`, 15 seeds): 400 subjects in 40 centres of 10, a
population shape of 0.6 and a centre `sd` of 0.35, with the truth
reject sampled to keep every true shape above 0.05. The smallest fitted
per-centre slope over the 15 was 0.0332 and the floor fired 0 times.

That is not luck, and the reason is worth writing down: at `df = 1` the
log density carries `log(d eta / d x)`, which is `log(gamma1 + u)`
exactly, so a centre with any event at all pays a penalty that grows
without bound as its slope approaches zero. `sp_floor_pos()` makes the
penalty finite rather than infinite, and measured on the shipped
`eps2 = 1e-14` it is 16.81 log units per event at a slope of exactly
zero, 29.02 at -0.01 and 32.01 at -0.2 (`frailty-barrier.R`), against a
normal prior on `u` that would spend a fraction of one.

### The centre with no events, which this lane got WRONG

The first draft of this document said the residual risk was a centre
with no events, "whose deviation is unidentified and shrinks toward
zero anyway". **That is false, and the punch round falsified it by
construction** (`dev/reviews/20260910-frailty.md`, section 5c). It is
recorded here rather than deleted, because it is the one claim in this
row that was asserted instead of constructed, and `dev/lane-rules.md`
says exactly that: "'No test reaches it' is NOT evidence that code is
unreachable."

The deviation of an all-censored centre is not unidentified and it does
not shrink toward zero. That centre contributes no density term, so
there is no `log(gamma1 + u)` barrier for it at all; what it does
contribute is `-H_i` per row, and the score in `u` is `-sum(x_i H_i)`,
which is NEGATIVE for rows past `t = 1`. The fit therefore pushes that
centre's slope DOWN and only the normal prior stops it. The
administrative time in this construction is 2.8199, so `log t` is
+1.0367 on every one of those rows, and 65.5 percent of all the
analysed times exceed 1, so the sign is not a hypothetical.

Re-derived here as `frailty-floor3.R`, seeds 20260910 to 20260915, on
`frailty-floor.R`'s own design (400 subjects, 40 centres of 10,
population shape 0.6, centre `sd` 0.35) with five centres followed to a
common administrative time and no deaths in any of them, which is an
ordinary thing for a small centre in a multi-centre trial:

| seed | population `gamma1` | min centre slope | centres at or below 0 | `n_nonmonotone` | refuses |
|---|---|---|---|---|---|
| 20260910 | 0.5623 | **-0.2852** | 5 | **0** | **FALSE** |
| 20260911 | 0.6848 | **-0.2583** | 5 | **0** | **FALSE** |
| 20260912 | 0.5498 | **-0.3056** | 5 | **0** | **FALSE** |
| 20260913 | 0.6446 | 0.1418 | 0 | 0 | FALSE |
| 20260914 | 0.5143 | **-0.2785** | 5 | **0** | **FALSE** |
| 20260915 | 0.5811 | 0.2106 | 0 | 0 | FALSE |

4 of 6 seeds. The control, the same six seeds with no such centre, is
`frailty-floor.R`'s own result and reproduces inside the same script:
minimum slope 0.0375 to 0.4286, nothing at or below zero.

At seed 20260910 the fit is silent the whole way through: `nlminb` code
0, maximum absolute gradient 2.74e-05, positive definite Hessian, no
warning and no message; `rp_floored(action = "report")$n_nonmonotone`
is **0**; `rp_floored()` does not refuse; `frm_curve()` passes it
through. Centre 1's fitted slope is -0.2127, so its fitted survival
RISES with time, 4.16e-50 at `t = 1e-12` to 0.774 at `t = 2.82`. **A
survival function that increases is not one, and the package's own
checker calls that fit clean.**

The mechanism is one PRE-EXISTING line in `R/rp-check.R`, which this
lane did not write and did not change:

    mono_rows <- which(f$cens == 0 & f$detadx <= 0)

Only uncensored rows are tested, and such a centre has none. Filed as
item 3 under "What is filed but not fixed"; NOT fixed here, for the
reason given there. What this lane owned was the sentence in
`?royston_parmar` telling a user that the likelihood holds the slope
positive and that `rp_floored()` reports it when it does not. That
sentence is corrected: the help page now says which rows the check
tests, shows this construction, and gives the two-line check to run
instead.

Which leaves the guard itself unexercised, and a zero that no run
produced is not a measured zero: a check that only ever read the
population `gamma1` would report the same 0 on every fit. So the case
is CONSTRUCTED (`frailty-floor2.R`, and it is now an assertion in
`test-frailty.R`). Push one centre's deviation past the population
slope, from `gamma1 = 1.3529` to a centre slope of -0.2, and
`rp_floored()` returns 15 non-monotone rows; those rows are exactly the
15 event rows of that one centre out of the 800 observed, and
`rp_floored()` refuses. The population value stayed positive
throughout, so nothing in that answer could have come from reading it.

### It is anchored at t = 1, and t = 1 is a unit

Rescale time by `c`. Then `x = log t` becomes `x' = x + log c` and

    gamma0 + (gamma1 + u) x = (gamma0 - (gamma1 + u) log c)
                              + (gamma1 + u) x'

so in the new units the same model carries a random INTERCEPT of
`-u log c` beside the slope. A block on `gamma1` alone therefore asserts
that every centre has the same cumulative hazard at `t = 1` exactly.

Measured, `frailty-gamma1-units2.R`, seed 20260910, 2000 rows, 40
centres, one dataset read in years, months and days. Log likelihoods are
compared as `logLik + n_event log c`, which is the Jacobian a density
owes a rescale: an event row contributes `log f` and `f` scales by
`1/c`, while a censored row contributes `log S`, which does not move.
Leaving that correction out was the first thing this measurement got
wrong, and it made both arms look non-invariant by 2982 units.

| model | spread over the three units | `sd(gamma1 \| centre)` at 1, 12, 365.25 |
|---|---|---|
| `gamma1 ~ (1 \| centre)` | 8.191 | 0.2201, 0.0756, 0.0370 |
| `(1 \| c \| centre)` on `mu` AND `gamma1` | 4.264e-08 | 0.225898, 0.225900, 0.225898 |

The paired block is invariant to 4.3e-08 log likelihood units over a
365x rescale and holds the slope standard deviation to 8.0e-06 relative,
while its intercept standard deviation moves from 0.151 to 1.381, which
is what the algebra above says it must. The slope-only model's estimate
of the very thing it is fitted for falls by a factor of 5.9 because the
user changed the units of the response.

Two CONTROLS were added in punch round 1, both in the same script. The
first licenses the Jacobian: a model with NO centre term is the same
model in any time unit by construction, so its corrected log likelihood
must be flat, and `~ trt` gives a spread of **2.1782e-09**. Without
that control a reader cannot separate a real anchoring from a Jacobian
I got wrong, and the first draft of this measurement did get it wrong,
by 2982 units. The second says what the anchoring is a property OF:
fixed per-centre slopes with a shared intercept, `gamma1 ~ centre`,
move **11.453** units over the same rescale. So the anchor has nothing
to do with the block being RANDOM. It follows from letting the slope
vary while holding the intercept common, however that is spelled, and
`?royston_parmar` now says so.

`?royston_parmar` says to pair the block unless the anchor is meant.

## What changed in the package

| file | change |
|---|---|
| `tests/testthat/test-frailty.R` | NEW. The measurement that replaces the smoke test: the change of basis against rstpm2, the estimate agreement relative to the run's own standard errors, the pin that `theta` is a variance, the decomposition of the log likelihood gap, the exact-integral optimality check (which needs no reference package at all), the `gamma1` recovery, the constructed monotonicity-floor case with the REACH of that check asserted beside it, and the anchoring ratio. 28 assertions in 4 tests. |
| `tests/testthat/test-surface.R` | the frailty smoke test is REMOVED from "a random effect and a smooth both reach a spline coefficient", which is now "a smooth reaches a spline coefficient". Two compat assertions added for the new rows. |
| `R/royston-parmar.R` | roxygen only: a new section, "A frailty, and a random effect on gamma1", carrying the numbers above. Punch round 1 corrected its monotonicity paragraph, which had told a user the likelihood holds a per-centre slope positive and that `rp_floored()` reports it when it does not. |
| `R/zzz.R` | two compat rows: royston_parmar against `us`, and against the ID-syntax key. |
| `DESCRIPTION` | `rstpm2` added to Suggests. |
| `vignettes/royston-parmar.Rmd` | one paragraph pointing at the new help-page section. No new fit, for the reason under "Not done". |
| `NEWS.md` | the entry, under an unreleased heading. |
| `dev/frailty/` | NEW. Every measurement script, plus a README, plus the RAW ROWS of the three replicate runs as `.tsv`, so the tables above can be spot-checked without a rerun. |
| `dev/frailty-findings.md` | NEW. This file. |
| `man/royston_parmar.Rd` | roxygenised. |

No `R/` code changed except `zzz.R`'s registration table, and exactly
one user-visible behaviour moved with it: `frm_compat("royston_parmar",
"|ID|")$status` goes from `untested` to `works`, and the `us` row now
carries the frailty note instead of the blanket covstruct one. The
punch round caught the first draft saying "no `R/` code changed" and
"no `R/` behaviour changed" flatly, which is false as written and is
contradicted by this document's own compat rows. No DENSITY, no
estimate and no fitted value changed: the review diffed every one of
the 56 objects in the installed namespace between the base build and
this one and found the two added `r(...)` calls and nothing else, and
five model families came back `identical()` across the two builds
(`dev/reviews/20260910-frailty.md`, section 1).

Nothing in this row needed a code fix to the family. The row was a
measurement, and what it found is in the THREE filed items below rather
than in a patch. The third of those is a real defect, in
`R/rp-check.R`, which this lane neither wrote nor changed and is
deliberately not fixing here.

### The suite, one test file per R process

`dev/frailty/frailty-suite.sh`, `NOT_CRAN=true`, against the package
installed in the lane's own library.

| file | pass | fail | error | skip |
|---|---|---|---|---|
| `test-bracket-access.R` | 1 | 0 | 0 | 0 |
| `test-curve.R` | 62 | 0 | 0 | 0 |
| `test-deriv.R` | 38 | 0 | 0 | 0 |
| `test-difference.R` | 59 | 0 | 0 | 0 |
| `test-frailty.R` | 28 | 0 | 0 | 0 |
| `test-gratia.R` | 13 | 0 | 0 | 0 |
| `test-message-uniqueness.R` | 4 | 0 | 0 | 0 |
| `test-royston-parmar.R` | 67 | 0 | 0 | 0 |
| `test-rp-floored.R` | 42 | 0 | 0 | 0 |
| `test-scale.R` | 0 | 0 | 0 | 1 |
| `test-span.R` | 73 | 0 | 0 | 0 |
| `test-surface.R` | 47 | 0 | 0 | 0 |

A count is not a count if the runner produced it wrongly, and the first
run of this suite produced four ERRORs that were the RUNNER's:
`testthat::test_file()` evaluates a file in the calling environment,
where a package's internal functions are not visible, while
`test_check()` uses `test_env(package)`, whose parent is the namespace.
Three files in this suite reach an internal (`sp_sim_crit()` and two
more), and they reported ERROR under the naive runner and pass under
the corrected one. Nothing in the package was involved.

### The instrument the sweep rests on, now asserted

`frailty-sweep.R` and `frailty-ident.R` score the Laplace objective at
a named parameter point through
`setNames(c(gam[1], beta, gam[-1], log(sd)), names(dry$obj$par))`.
`setNames()` relabels and does not reorder, so that line ASSERTED an
ordering which nothing checked, and the objective's own names are
`beta, beta, betad, betad, theta`: duplicated, so they give no
protection. Every offset in this document is scored through it. The
punch round found the assumption right and unasserted.

Both scripts now carry `check_par_order()`, which compares the Laplace
objective at the fitted estimates with the fit's own `logLik()` and
stops if they differ by more than 1e-8 relative. As shipped it reads
5.16e-15. `frailty-parorder.R` runs it against the case it exists to
catch, at seed 20260910:

| ordering | `-obj$fn` | relative to `logLik()` |
|---|---|---|
| as the scripts build it | -2645.40957653 | **5.16e-15** |
| the covariate ahead of the intercept | -3891.77738236 | 4.71e-01 |
| the two gammas swapped | -85054.21066222 | 3.12e+01 |
| the `mu` and `gamma` blocks exchanged | -64839.43992926 | 2.35e+01 |

Seven orders of magnitude of daylight on either side of the threshold,
so the guard fails closed.

### R CMD check --as-cran

`dev/frailty/frailty-check.sh`: `R CMD build`, then
`R CMD check --as-cran --no-multiarch` on the tarball, with pandoc and
TinyTeX on PATH and `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, as
`dev/lane-rules.md` prescribes.

    Status: 1 NOTE
    * checking HTML version of manual ... NOTE
      Skipping checking math rendering: package 'V8' unavailable

That is the expected NOTE and the only one. Inside the check the suite
reports `FAIL 0 | WARN 0 | SKIP 2 | PASS 430`, the two skips being
`test-scale.R`'s gate and the one `test-message-uniqueness.R` assertion
that needs the package sources rather than an installed package.
Vignettes rebuild with no output.

It was run TWICE, not once, and the SECOND is the one that counts. The
first was the lane's own final pass, at `PASS 426`; punch round 1 then
changed `man/royston_parmar.Rd` and added four assertions to
`tests/testthat/test-frailty.R`, which is package content the check
reads, so leaving the first result standing would have been reporting a
stale count. Both give `Status: 1 NOTE`, and the count moved 426 to 430
exactly as the four new assertions predict.

## The guards, run against the case they exist to catch

`frailty-falsify.R` and `frailty-parorder.R`. A check that cannot fail
is not a check.

| guard | broken arm holds | as shipped |
|---|---|---|
| `theta` read as a standard deviation | FALSE | TRUE |
| the two spline spaces differ (one knot moved 0.6) | FALSE | TRUE |
| a point one standard error off the estimate called an optimum | FALSE | TRUE |
| an inert `gamma1` block still beats the pooled fit | **TRUE** | TRUE |
| a time-unit-invariant model called anchored | FALSE | TRUE |
| `rp_floored()` reading only the population `gamma1` | FALSE | TRUE |
| the objective's parameter ordering, three wrong ones | FALSE | TRUE |

Rows 1 to 5 are `frailty-falsify.R`; row 6 is the constructed floor
case in the same script; row 7 is `frailty-parorder.R`, added in punch
round 1.

Row 4 failed open on its first spelling and the test was changed, as
described above. The others fail closed. Row 6 needs its construction
because `rp_floored()` returned 0 on all 75 fits this lane made: with
nothing to observe, the broken arm has to be the check reading the
population value, and the constructed arm is what separates the two
readings.

One thing in `test-frailty.R` is NOT a guard and should not be read as
one: the assertion that an all-censored group comes back with
`n_nonmonotone = 0`. That pins the GAP from the blocker, so that the
reach of the check is written down rather than inferred, and it is
meant to FAIL when the filed widening of `mono_rows` lands. The test
says so beside it, because a test that silently outlives the behaviour
it pins is worse than none.

It does fail when that lands, and that was checked rather than hoped:
on the test's own `bad2` the narrow rule finds 0 rows and the widened
rule finds 40, every row of centre 3, so `expect_equal(..., 0L)` sees
40 and `expect_silent()` sees a refusal, and neither line can keep
passing quietly after the fix (`dev/reviews/20260910-frailty.md`,
re-check). Worth recording for a reason beyond this row: the re-check
calls it **the first guard in four rounds of this project to fail
closed on its first spelling.** Every guard built in the 0.55.1 round
failed open on its first try, and `dev/lane-rules.md` carries that as a
standing warning. It is possible to get one right first time, and what
made the difference here was writing the inverse case at the same time
as the assertion rather than after it.

## What I did NOT do, and why

- **No vignette fit for the frailty.** `bc`, which the vignette uses,
  has no centre variable, and the twelve labels the old smoke test dealt
  round robin carry no centre effect. A simulated dataset inside the
  vignette would be a second, weaker copy of a measurement that now
  lives in the help page and in this file. The vignette gets a
  cross-reference instead.
- **No refusal or warning for `gamma1 ~ (1 | centre)` without a block
  on `mu`.** The model is legitimate, not wrong, and a user who has
  centered log time deliberately means it. Punch round 1 strengthened
  this: the anchoring is not a property of the block being RANDOM at
  all. Fixed per-centre slopes with a shared intercept,
  `gamma1 ~ centre`, move 11.453 log likelihood units over the same
  365x rescale where the paired random block moves 4.3e-08. A refusal
  keyed on the random-effect spelling would therefore miss the
  fixed-effect one and refuse a model whose fixed twin it accepts. The
  behaviour belongs to a varying slope with a shared intercept in ANY
  model, so a check belongs in core if it belongs anywhere. It is
  documented in all three places a user looks, and the measurement that
  would justify a check is in this file if core ever wants it.
- **No scale-tier row.** The plan's Phase 0 table says the spline row
  measures the curves design and that the survival design belongs to
  this item. The fits here are 0.7 to 3.1 s at 2000 rows (frmtmb) and
  1.6 to 5.1 s (rstpm2), so there is no cost finding to record and a
  `test-scale.R` row would only repeat the recovery table more slowly.
- **No `check_laplace()` run.** It lives in `frmtmb.sample` and needs
  Stan. The exact per-cluster integral answers the same question here
  for a few seconds of R, and it answers it at a named parameter point
  rather than through a second fit.
- **No delayed entry.** `trunc()` on this family is still `conditional`
  for the reason its compat row gives, and a frailty does not change
  that. It is item 4.1's.
- **No widening of `rp_floored()`**, which is the code fix the blocker
  above points at. `mono_rows <- which(f$detadx <= 0)` is one line and
  it is deliberately NOT in this row. It changes a shipped refusal, and
  a check that fires on a correct model is worse than no check, so it
  needs a measured false-alarm rate on designs the field produces
  before it can ship. That is its own lane, and the coordinator is
  scheduling it. Item 3 under "What is filed but not fixed" states it
  so that it can be lifted into the plan verbatim, with the
  construction and with what the fix would have to measure.

## What is filed but not fixed

1. **The Laplace maximizer is systematically displaced**, by 5.5e-05
   exact log likelihood units at 50 subjects per centre, on 200 of 200
   replicates in the same direction, that unanimity being structural
   for the reason given above. It is 0.1 percent of the Laplace
   offset and 1e-3 of a standard error, so it is not worth a remedy at
   this design, but it is the quantity that would grow first if someone
   fits many small clusters. At 4 subjects per centre it is 3.3e-03
   units. A higher-order correction is a core question (RTMB's
   `random =`), not a family one.
2. **`sd(frailty)` sits low**, and the punch round settled how far.
   Over 200 replicates it is 0.4905 against 0.5, 1.95 Monte Carlo
   standard errors, which is all 200 replicates can resolve. Over 800
   (`dev/reviews/20260910-frailty.md`, section 4, the first 200 seeds
   being this lane's own) it is 0.4889 at a Monte Carlo standard error
   of 0.0023, so **4.8 standard errors**: the shortfall is real and it
   is not the integration rule, since rstpm2's quadrature estimate is
   only 0.1386 percent higher on 200 of 200. It is the ordinary
   downward bias of a maximum likelihood variance component, which for
   `m = 40` groups and `p = 2` fixed effects predicts
   `0.5 * sqrt(38 / 40) = 0.4873` against the observed 0.4889, leaving
   0.3 percent unexplained. The asymmetry of the misses says the same:
   of the 46 `sd` intervals that missed, 36 lay entirely below 0.5 and
   10 above. Not a defect; recorded with its explanation so that the
   next reader does not re-derive it.

3. **`rp_floored()` does not see a group with no events, and a
   Royston-Parmar fit can put one on a decreasing hazard.** This is a
   defect in `R/rp-check.R`, pre-existing, not written or changed by
   this lane, and stated here so that it can be lifted into
   `dev/extension-gaps-plan.md` verbatim.

   *What is wrong.* `mono_rows <- which(f$cens == 0 & f$detadx <= 0)`
   tests only uncensored rows. A random effect on `gamma1` gives each
   group its own slope in log time, and for a group whose rows are all
   censored past `t = 1` the score in the deviation is
   `-sum(x_i H_i) < 0`, so the fit drives that slope negative with no
   `log(gamma1 + u)` barrier to stop it, because the barrier is in the
   density and the group contributes none. The check then has no row to
   test and returns a count of zero.

   *The construction.* `dev/frailty/frailty-floor3.R`, seeds 20260910
   to 20260915: 400 subjects, 40 centres of 10, population shape 0.6,
   centre `sd` 0.35, five centres followed to a common administrative
   time with no deaths. On 4 of 6 seeds those five land at slopes of
   -0.21 to -0.31. At seed 20260910 `nlminb` returns code 0 at a
   maximum absolute gradient of 2.74e-05 with a positive definite
   Hessian, `frm()` warns nothing, `rp_floored(action = "report")`
   returns `n_nonmonotone = 0`, `rp_floored()` does not refuse,
   `frm_curve()` passes the fit through, and centre 1's fitted survival
   rises rather than falls, 4.16e-50 at `t = 1e-12` against 0.774 at
   `t = 2.82`. The control, the same design with no such centre, keeps
   every fitted slope between 0.0375 and 0.4286 on all six seeds. The
   administrative time is 2.8199 at this seed, so `log t` is +1.0367
   for every one of those rows and the sign argument above is not a
   hypothetical: 65.5 percent of the analysed times exceed 1.

   *What a fix would have to measure.* The one-line widening,
   `mono_rows <- which(f$detadx <= 0)`, changes a shipped REFUSAL, and
   a check that fires on a correct model is worse than no check. So it
   needs a false-alarm rate on designs the field produces, not a
   patch: the rate on ordinary right-censored survival data where the
   fitted spline is monotone everywhere, at the `df` values flexsurv
   users pick. (A draft of this item also named a fit monotone across
   the EVENT times but bending past the last one, with a censored row
   beyond the boundary knot. That case cannot occur under the default
   knot rule, and the reason the draft gave is why: beyond the boundary
   the basis is LINEAR, so its derivative is constant, and the default
   rule puts `bknots = range(log t | event)`, which places the last
   event row AT `kmax` with the same `detadx` as every row past it. The
   narrow check already covers that region. Pinning `bknots[2]` 1.5
   above the last event makes the region genuinely cubic, and a
   25,921-point sweep over `gamma2` by `gamma3` with 320 censored rows
   in the window found no vector with every event row positive and any
   censored row non-positive. The remaining route, a censored row in a
   dip between two event times inside the interior knots, was not
   reached either.) It also has to say what the refusal means for a
   fit whose reported likelihood is NOT floored: `rp_floored()`'s
   charter is "this fit's reported likelihood is not the model's", and
   `sp_floor_pos()` never touches a censored row, so under the current
   charter this fit's likelihood IS the model's and the widened check
   would be answering a different question. That question, "is the
   fitted hazard a hazard everywhere the data are", may deserve its own
   name rather than an extra row in this one.

## The punch rounds, and what they changed

`dev/reviews/20260910-frailty.md`, 842 lines over two rounds. Round 1
returned "mergeable with one named fix", round 2 "mergeable". Every
headline figure in this document reproduced to the digit under the
reviewer's own library and its own scripts, including the twelve test
counts, the three shipped `.tsv` files and the parameter-ordering
guard. What changed here:

| finding | what it was | what it is now |
|---|---|---|
| BLOCKER, the all-censored centre | asserted safe, "shrinks toward zero anyway" | falsified by construction, `frailty-floor3.R`; the help page corrected, the check filed as item 3 |
| "no `R/` behaviour changed" | flatly false for `frm_compat()` | says what moved, in the findings and in NEWS |
| the `sd` shortfall | 1.95 Monte Carlo errors, which is all 200 replicates could see | 4.8 over 800, with the ML bias factor that predicts it to 0.3 percent |
| coverage 0.930 and 0.920 | read as what 200 replicates show | CONFIRMED: 0.9450 and 0.9425 over 800, and exactly 0.9500 on the 600 seeds this lane never ran |
| the 200-centre sweep cell | a bare mean | quoted with its range, -0.0105 to -0.1923 |
| the "200 of 200" | read as 200 confirmations | said to be structural, with the magnitude leading |
| the parameter ordering | assumed | asserted in both scripts, and falsified against three wrong orderings |
| the time-unit measurement | two arms | two CONTROLS added, one of which shows the anchor is not a random-effect property |
| `confint_varcorr()` | frmtmb 0.13.0 | 0.12.0, and the conclusion strengthened |
| `bad$estimates$b` | `$` on a container with three `b` names | `[["b"]]` |
| the raw replicate rows | regenerated | shipped, three `.tsv` files |
| the two-line remedy (round 2) | silent about `df` | said to be the `df = 1` form, with the exposed set named |
| the `ranef()` column (round 2) | "that column", singular | named, since the paired block gives two |
| filed item 3's "bends beyond the last one" (round 2) | a case that cannot occur | struck, with the reason and the sweep that failed to reach it |

Punch round 2 was three prose clauses and the verdict was MERGEABLE.
The two-line remedy in `?royston_parmar` now says it is the `df = 1`
form and names the set it leaves exposed, a group with NO events at
`df >= 2`, because above `df = 1` the derivative carries the interior
coefficients and the two lines err OPTIMISTIC: on a `df = 3` fit they
read 1.207 to 2.162 where the real derivative runs 0.768 to 2.017,
0.636 to 0.796 of it per centre. The remedy also names the `ranef()`
column, `time.gamma1:(Intercept)`, since the same page recommends the
paired block and that gives two columns of which the `mu` one is
irrelevant to monotonicity. And filed item 3 lost a clause that was
wrong: see the parenthesis in it.

That round also confirmed the corrected paragraph does not overclaim in
the other direction, on 24 fits this lane did not run: giving those
five centres 1, 2 or 3 events each turns every seed positive, minimum
slope +0.0367 with no negative centre on 6 of 6, against -0.3056 and
five negative centres on 4 of 6 at zero events. ONE event is enough,
which is stronger than "a centre with any event at all pays a barrier"
was able to claim from the barrier arithmetic alone.

One thing the reviews raised and this lane deliberately left for the
consolidating session: the base compatibility table covers `us`,
`diag`, `homdiag`, `cs` and the rest with one blanket note, and this
row replaces that note for `us` alone. A reader comparing `us` with
`diag` now sees one measured row and one blanket row for two features
with the same evidence behind them. Nothing is wrong; whether the
blanket note should point at the measured one is a table-wide decision
and not this row's.

## Reproducing every number

The scripts are in `dev/frailty/`, all prefixed `frailty-`, with a
README beside them, and so are the raw rows: `frailty-recover.tsv`
(200 replicates), `frailty-gamma1.tsv` (60) and `frailty-sweep.tsv`
(30 cells), one `key=value` line per replicate. The punch round asked
for those, because a table that cannot be spot-checked without a
26-minute rerun is a table a reader has to take on trust.

All of it is in the repository rather than in the session scratchpad on
purpose: `dev/lane-rules.md` records a process
document that was lost when a scratchpad was cleaned mid-round, and a
number whose construction is gone gets re-derived from scratch. Each
sources `frailty-lib.R` for the library order and `frailty-common.R`
for the simulator, the hand-written basis, the exact integral and the
change of basis, all by bare name, so run them from inside that
directory.

| script | what it produces | seeds |
|---|---|---|
| `frailty-ident.R` | the first side-by-side, the three rules at one point, rstpm2 at 3/9/21/51 nodes | 20260910 |
| `frailty-sweep.R`, `frailty-sweep-report.R` | the cluster-size and frailty-size table, `frailty-sweep.tsv` | 20260910 to 20260914 |
| `frailty-recover.R`, `frailty-recover-report.R` | the 200-replicate recovery and agreement tables, `frailty-recover.tsv` | 20260910 to 20261109 |
| `frailty-gamma1-recover.R`, `frailty-gamma1-report.R` | the 60-replicate `gamma1` table, `frailty-gamma1.tsv` | 20260910 to 20260969 |
| `frailty-gamma1-units2.R` | the time-unit measurement | 20260910 |
| `frailty-falsify.R` | the five guards against their broken arms | 20260910 |
| `frailty-testcost.R` | what a test-sized design costs, and the covstruct key a `(1 \| centre)` block resolves to | 20260910 |
| `frailty-floor.R` | 15 fits on a design built to reach the monotonicity floor, which did not reach it | 20260910 to 20260924 |
| `frailty-floor2.R` | the constructed floor case, which is what pins `rp_floored()` to the centre deviations | 20260910 |
| `frailty-barrier.R` | the size of the monotonicity penalty per event, from `sp_floor_pos()`'s own arithmetic | none |
| `frailty-floor3.R` | the all-censored centre, which is the case this lane asserted and the punch round falsified, with its control | 20260910 to 20260915 |
| `frailty-parorder.R` | the parameter-ordering guard against three wrong orderings | 20260910 |

Library: a private `C:/Users/adf44/source/r/frailtylib` first, then the
round's shared `rellib-r3`, then `pinlib`, then the user library.
Nothing was installed into any of the last three. **rstpm2 1.7.1 was
already present in the user library**, so no install was needed at all;
the lane's instruction to install it was not acted on because acting on
it would have been a second copy of a package that was already there.

One trap for whoever runs rstpm2 next: `stpm2()` rewrites its own call
to `gsm(...)` and evaluates it in the CALLER's frame, so it fails with
"could not find function gsm" unless rstpm2 is attached. The test file
calls `rstpm2::gsm()` directly, which is what `stpm2()` calls, and a
test file should not attach a suggested package. `logLik()` on the
result also needs bbmle attached to dispatch; the `@min` slot holds the
same number and does not.

## Version

This row adds tests, documentation and two compat rows, and changes no
`R/` behaviour. It needs a frmtmb.spline version bump and a NEWS
heading to go with it; the number is the consolidating session's to
choose, and the NEWS entry is written under an unreleased heading for
it to fill in.

No new floor on frmtmb is needed, and the punch round gave a stronger
route to that than naming symbols: frmtmb 0.55.1 and 0.55.2 are BOTH
documentation-only releases in core's `NEWS.md`, so nothing added after
the declared floor `frmtmb (>= 0.55.1)` exists to be reached. Naming
the symbols agrees: the newest one anything here touches is
`confint_varcorr()`, which core's `NEWS.md` puts at frmtmb **0.12.0**
(the first draft of this document said 0.13.0 and was wrong), and the
rest is `frm()`, `bf()`, `cens()`, `fixef()`, `ranef()`, `VarCorr()`,
`frm_compat()`, `diagnose()` and the `(1 | c | g)` spelling.

`rstpm2` joins Suggests, and every test that needs it is behind
`skip_if_not_installed("rstpm2")`. The one comparison test is the only
thing in the suite that needs it; the exact-integral test that measures
the Laplace offset needs no reference package at all.
