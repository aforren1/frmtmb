# Lane wt-shapes: brms's return shapes and brms's predict()

Items 2.6f and 2.6d together. Worktree
`C:\Users\adf44\source\r\frmtmb-wt-shapes`, base `f8b45ef` (frmtmb
0.60.0, frmtmb.sample 0.8.0). Private library
`C:/Users/adf44/source/r/shapes-lib`, built by `dev/shapes-install.R`;
the round's shared read-only base build is
`C:/Users/adf44/source/r/rellib-r3` and nothing is installed into it.
There is no `pinlib` on any path: StanHeaders 2.39.1 with tmbstan 1.2.1
is correct now, and `dev/brmsport-run.R` asserts tmbstan >= 1.2.1 and
the `-std=gnu++17` flag in the user Makevars instead of the old pin.

**Punch round 1 (2026-09-18).** `dev/reviews/20260918-shapes.md`
returned NOT MERGEABLE with two blockers and five majors. Everything
it named is fixed, and three of the corrections are CORRECTIONS OF
THIS DOCUMENT rather than of the code: section 5's "every number
agrees with the base build" was false and the instrument that produced
it was the reason; section 8's claim that multivariate `fitted()`
"stays refused" was wrong either way; and section 8's one-sentence
account of the sdreport cost undercounted it. Each is marked
**CORRECTION** where it stands, with what was wrong and how it was
measured, rather than quietly rewritten. Sections 3, 5a, 6 and 9 carry
the punch round's own measurements.

## 1. What the two items are, in one sentence each

**2.6d.** `predict()` carried glmmTMB's `type` vocabulary and returned
the LINEAR PREDICTOR; brms's `predict()` summarizes the PREDICTIVE
distribution. brms is the tiebreaker, so `predict()` is brms's and the
old function moved to `frm_linpred()` under its own name.

**2.6f.** `fitted()`, `residuals()`, `fixef()`, `ngrps()`, `vcov()` and
`summary()` returned lme4 and glmmTMB shapes; they now return brms's.
frmtmb.sample gains `fitted()`, `predict()` and `residuals()` on draws,
answers `nsamples()`, `posterior_samples()` and `parnames()` instead of
refusing them, and honors `point_estimate`.

## 1a. What each method returns now

| method | before | after |
|---|---|---|
| `predict()` | the linear predictor, a vector (a list under `se.fit`) | brms's `n x 4` predictive summary; `n x K` `P(Y = k)` proportions for a category response; brms's `n x 4 x nresp` array for a multivariate fit; the draws under `summary = FALSE` |
| `predict(allow_new_levels = TRUE)` | the unseen level's effect held at ZERO, so the interval was a known level's | that effect DRAWN from the block's estimated covariance per replicate, which is brms's `sample_new_levels = "gaussian"` |
| `frm_linpred()` | did not exist | `predict()` of 0.60.0, body and arguments unchanged |
| `fitted()` | a vector; `n x K` probabilities for a category response | brms's `n x 4` matrix; `n x 4 x K` array for a category response |
| `residuals()` | a vector | brms's `n x 4` matrix; `n x 4 x K` for a matrix-valued response |
| `fixef()` | a list per dpar, named by design column | brms's summary matrix, brms's row names, brms's ORDER |
| `fixef_by_dpar()` | did not exist | `fixef()` of 0.60.0 |
| `vcov()` | every estimated coefficient, internal names, 10 x 10 on the port's fixture | brms's population-level block, brms's names, 9 x 9 |
| `fixef()`, `vcov()`, `summary()$fixed` on an ORDINAL fit | the coefficients only, `x` and 1 x 1 | the thresholds and `cs()` coefficients too: `Intercept[1]`, `Intercept[2]`, `x` and 3 x 3, which is brms's |
| `coef()` | design-column names, `(Intercept)` | brms's names, `Intercept`, the same ones `fixef()` and `vcov()` use |
| `coef()` on a mixed ordinal fit | thresholds repeated per group plus a stray `(Intercept)` column (round 1) | each threshold moved by the group's intercept with brms's sign (punch round 2) |
| `ranef()` | `(Intercept)` | brms's names, `Intercept` (punch round 2) |
| `predict(summary = FALSE)` on a rescor fit | each response drawn alone (round 1) | the joint multivariate normal (punch round 2) |
| `insight::get_residuals()` | a length-n vector, through insight's default | the same vector, through a method of frmtmb's own, with insight's class and the matrix on a `"full"` attribute |
| `hypothesis()` class | `c("frmtmb_hypothesis", "brmshypothesis")` | `"frmtmb_hypothesis"` |
| `vcov_estimated()` | did not exist | `vcov()` of 0.60.0, for the interop seams |
| `ngrps()` | lme4's named integer vector, counting smooth-like blocks | brms's named list, or `NULL`, grouping factors only |
| `summary()` | frmtmb's slots | brms's `$fixed`, `$random`, `$spec_pars`, `$cor_pars` beside them; brms's print headings |
| `print(fit)` | frmtmb's own layout | the summary, as brms's does |
| `variables()` on an ordinal fit | 3 of 9 names | every one, thresholds as `b_Intercept[k]` and `cs()` as `bcs_<term>[k]` |
| `fitted()`, `predict()`, `residuals()` on DRAWS | `NULL`, through `stats::fitted.default()` | brms's summaries of `posterior_epred()`, `posterior_predict()` and `predictive_error()` |
| `nsamples()`, `posterior_samples()`, `parnames()` | refused | brms's answer and brms's deprecation warning |
| `point_estimate`, `ndraws_point_estimate` | accepted and ignored | honored, as a parameter-space collapse |

## 2. The frequentist construction of each brms column

brms fills `Estimate`, `Est.Error` and two `Q` columns from posterior
draws. A maximum-likelihood fit has no posterior, so each method's
columns are the frequentist analogue, and which analogue is a decision
per method rather than one rule:

| method | Estimate | Est.Error | Q columns |
|---|---|---|---|
| `fixef()`, `summary()$fixed` | the coefficient | its standard error, `sqrt(diag(vcov()))` | Wald |
| `fixef()`, an ordinal fit's `Intercept[k]` | the threshold on the model's own scale | the delta method over the joint covariance, because the family may estimate `(tau_1, log increments)` rather than the thresholds | Wald |
| `fitted()` | the modelled response | the delta-method standard error `frm_linpred(se.fit = TRUE)` reports | Wald |
| `fitted()`, ordinal or categorical | each category probability | a finite-difference delta method over the whole outer parameter vector | Wald, held inside [0, 1] because the quantity is a probability |
| `residuals()` | `y` minus the fitted value | the standard error of the fitted value, since `y` is fixed | Wald |
| `residuals(type = "pearson")` | that, over the residual SD | the same, over the same SD | Wald |
| `residuals(type = "deviance"/"osa")` | the residual | `NA`: neither has a standard error here | `NA` |
| `residuals()`, matrix-valued response | `y` minus the fitted cell | the finite-difference delta method, as for a category probability | Wald |
| `summary()$random`, `$spec_pars`, `$cor_pars` | the natural-scale value | the delta method through the transform | the transformed-scale Wald interval mapped back, so it cannot leave the parameter's range |
| `predict()` | the mean of the simulated draws | their standard deviation | the draws' own empirical quantiles |

`predict()` is the one with real draws, because it simulates, so its
interval carries the family's skew and its discreteness rather than
being symmetric: a poisson predictive interval is on the counts.

## 3. predict(): the construction, and what it covers

Each replicate draws the outer parameter vector from
`N(theta_hat, vcov(fit, full = TRUE))`, evaluates every distributional
parameter there and draws a response from the family's own simulator,
the one `simulate()` and `posterior_predict()` use. That is a
parametric bootstrap of the predictive distribution.
`param_uncertainty = FALSE` drops the parameter draw, which is the
plug-in predictive distribution.

The random effects of a level the fit SAW stay at their conditional
modes, the convention every other method in this package follows, so
`Var(b | y)` is NOT in the interval. That is the term the mixed arm
below is sensitive to, and it is stated in `?predict.frmtmb_fit`.

A level the fit did NOT see has no mode to condition on, so
`allow_new_levels = TRUE` draws that level's effect from the block's
own estimated covariance, once per replicate, at THAT replicate's
parameters. This is brms's `sample_new_levels = "gaussian"`. The
design pieces come from `lp_extra_var()` and `extra_var_blocks()`, the
two functions `frm_linpred(se.fit = TRUE)` already used for the same
term, so the variance the draw uses and the variance the standard
error reports cannot drift apart. The `newlevel` arm below is the
measurement.

**The measurement.** Out of sample, which is the predictive question:
fit on 60 rows, form the interval at 10 FRESH design points, draw the
responses there from the same truth, and count. 220 replicates of 10
points is 2200 predictions per arm; 202 replicates is where 0.80 power
is first SUSTAINED for a one-sample count against a rate known exactly
(`dev/lane-rules.md`). The predictions inside one replicate share a
fit, so the pooled binomial interval is too narrow and the
replicate-level standard error is reported beside it.

Four designs and up to four arms each. Two of the arms are CONTROLS
and neither is optional:

- `wald` is `fitted()`'s interval, which carries no observation noise
  at all. It MUST under-cover badly, so that a coverage near 0.95 in
  the other arms cannot be an artifact of the harness.
- `exact` (gaussian design) is the textbook prediction interval from
  `lm()`, `t_{n-p} sigma sqrt(1 + h)`. It MUST cover at the nominal
  rate, so that the harness is shown to be able to REPORT 0.95 and not
  only to report a shortfall. A negative control alone cannot show
  that.

and one is a TEST of an attribution rather than a proposal:

- `condvar` (mixed design) is the `joint` arm widened by the analytic
  conditional variance of the group mode,
  `sigma^2 tau^2 / (sigma^2 + n_g tau^2)`. The documentation says the
  mixed shortfall IS the missing `Var(b | y)`. If adding exactly that
  term moves the coverage by about the shortfall, the attribution is
  measured; if it does not, the attribution is wrong.

**The Monte Carlo figures below are punch round 2's.** Round 2 changed
how `predict()` spends random numbers (all parameter draws and one seed
per replicate taken up front, so that a row's parameter draws do not
depend on how many rows there are), so every figure moved by Monte Carlo noise, within one
replicate-level standard error on every arm; round 1's outputs are kept
in `dev/shapes-log/cov-r1/`. The controls did not move: `wald` and
`exact` do not simulate.

<!-- BEGIN GENERATED: dev/shapes-coverage.R -->

```
design: gauss  replicates: 220  failed fits: 0
n: 60  new points per replicate: 10  ndraws: 1000  seed0: 20260917  allow_new_levels: FALSE
script: dev/shapes-coverage.R

joint     2077 of  2200 = 0.9441  (0.9337, 0.9533) binomial; replicate mean 0.9441 se 0.0049  median width 3.9451  non-finite widths 0 of 220
plugin    2067 of  2200 = 0.9395  (0.9288, 0.9491) binomial; replicate mean 0.9395 se 0.0053  median width 3.8261  non-finite widths 0 of 220
wald       577 of  2200 = 0.2623  (0.2440, 0.2812) binomial; replicate mean 0.2623 se 0.0092  median width 0.6730  non-finite widths 0 of 220
exact     2101 of  2200 = 0.9550  (0.9455, 0.9633) binomial; replicate mean 0.9550 se 0.0043  median width 4.0737  non-finite widths 0 of 220
```

```
design: pois  replicates: 220  failed fits: 0
n: 60  new points per replicate: 10  ndraws: 1000  seed0: 20260917  allow_new_levels: FALSE
script: dev/shapes-coverage.R

joint     2162 of  2200 = 0.9827  (0.9764, 0.9877) binomial; replicate mean 0.9827 se 0.0029  median width 4.6000  non-finite widths 0 of 220
plugin    2161 of  2200 = 0.9823  (0.9758, 0.9874) binomial; replicate mean 0.9823 se 0.0028  median width 4.5000  non-finite widths 0 of 220
wald       540 of  2200 = 0.2455  (0.2276, 0.2640) binomial; replicate mean 0.2455 se 0.0084  median width 0.8811  non-finite widths 0 of 220
```

```
design: mixed  replicates: 220  failed fits: 0
n: 60  new points per replicate: 10  ndraws: 1000  seed0: 20260917  allow_new_levels: FALSE
script: dev/shapes-coverage.R

joint     2078 of  2200 = 0.9445  (0.9341, 0.9537) binomial; replicate mean 0.9445 se 0.0059  median width 4.0616  non-finite widths 0 of 220
plugin    2047 of  2200 = 0.9305  (0.9190, 0.9407) binomial; replicate mean 0.9305 se 0.0063  median width 3.9068  non-finite widths 0 of 220
wald      1170 of  2200 = 0.5318  (0.5107, 0.5528) binomial; replicate mean 0.5318 se 0.0118  median width 1.5747  non-finite widths 0 of 220
condvar   2114 of  2200 = 0.9609  (0.9519, 0.9686) binomial; replicate mean 0.9609 se 0.0049  median width 4.2690  non-finite widths 0 of 220
```

```
design: newlevel  replicates: 220  failed fits: 0
n: 60  new points per replicate: 10  ndraws: 1000  seed0: 20260917  allow_new_levels: TRUE
script: dev/shapes-coverage.R

joint     2095 of  2200 = 0.9523  (0.9425, 0.9608) binomial; replicate mean 0.9523 se 0.0048  median width 4.9822  non-finite widths 0 of 220
plugin    2046 of  2200 = 0.9300  (0.9185, 0.9403) binomial; replicate mean 0.9300 se 0.0060  median width 4.5988  non-finite widths 0 of 220
wald      1531 of  2200 = 0.6959  (0.6762, 0.7151) binomial; replicate mean 0.6959 se 0.0144  median width 2.6241  non-finite widths 0 of 220
```

<!-- END GENERATED -->

**What the four designs say.**

- The negative control fails as it must: 0.2623, 0.2455, 0.5318 and
  0.6959 against a nominal 0.95. `fitted()`'s interval is around the
  expected response and is not a predictive interval; the numbers say
  how far apart the two questions are.
- The positive control passes as it must: `exact`, the textbook
  `lm()` prediction interval on the gaussian design, covers 0.9550
  with a binomial interval of (0.9455, 0.9633) that contains 0.95. The
  harness can report 0.95.
- `joint` covers at the nominal rate on the gaussian design, 0.9441
  with a binomial interval of (0.9337, 0.9533) that contains 0.95. It
  sits just below the `exact` arm's 0.9550 and its interval is just
  narrower, which is the normal-versus-t difference at n = 60 and the
  Monte Carlo error of 1000 draws, not a defect in the construction.
- `joint` covers at least as well as `plugin` on every design, by
  0.0046, 0.0004, 0.0140 and 0.0223. That is the direction the theory gives: the plug-in
  interval ignores the uncertainty in the estimates, so it is too
  narrow, and the gap grows with the number of parameters relative to
  n. It is why that draw is on by default. (The argument was
  `param_uncertainty` when this was measured and is `propagate_error`
  since 2026-09-23; `dev/reunc-findings.md` section 13 has the rename.
  The other mentions in this document are measurements taken under the
  old name and are left as they were recorded.)
- The poisson design OVER-covers, 0.9827, and that is not a defect of
  the construction. A predictive interval for a count is read off a
  LATTICE: the 2.5 and 97.5 percent empirical quantiles land on whole
  numbers, so the realized coverage is at least nominal and usually
  above it. brms's `predict()` has the same property for the same
  reason. The `plugin` arm is 0.9823, still above nominal, so the
  conservatism is the discreteness and not the parameter draw.

**The unseen-level arm, `newlevel`.** Every fresh point sits at a
grouping level the fit never saw, and the truth draws a fresh group
effect per level, so the correct predictive spread is
`sqrt(sigma^2 + tau^2)` rather than `sqrt(sigma^2)`. Before punch
round 1 the effect was held at ZERO: the reviewer measured coverage
**0.8618** there, and `Est.Error` at an unseen level bit-identical to
a known level's, 0.84358 against 0.84842
(`dev/reviews/20260918-shapes.md`, BLOCKER 2, seed base 883000).
**CORRECTION (punch round 2):** the review gave 1.02794 as the right
value there, and the recheck found that target itself wrong: it
counted `tau^2` twice. The plug-in value is `sqrt(se^2 + sigma^2)`
with the block variance inside `se^2` once, 0.936 and 0.944 on the
review's two fits; the finding stands, since the unseen level still
had no between-group variance at all. With
the draw in, this design's `joint` arm covers **0.9523**, with a
binomial interval of (0.9425, 0.9608) that contains 0.95 and a
replicate-level standard error of 0.0048, so 0.5 standard errors
above nominal. The `wald` control on the same design is 0.6959, which
is higher than the mixed design's 0.5318 for the right reason:
`fitted()` already added the block's marginal variance at an unseen
level, and the gap that remains is the observation noise.

**The mixed shortfall, and its attribution, now TESTED.** The mixed
design is where the construction is weakest, and the shortfall is
real rather than a rounding of it:

| arm | coverage | replicate se | z against 0.95 | binomial CI |
|---|---|---|---|---|
| `joint` | 0.9445 | 0.0059 | **-0.93** | (0.9341, 0.9537) |
| `plugin` | 0.9305 | 0.0063 | -3.10 | (0.9190, 0.9407) |
| `condvar` | 0.9609 | 0.0049 | **+2.22** | (0.9519, 0.9686) |

The reviewer ran an independently constructed version of the same
design, SEED0 = 771000 rather than this lane's, and got `joint`
0.9336 with z = -2.80 and a binomial interval that EXCLUDES 0.95, and
`joint_plus_condvar` 0.9477 with z = -0.44
(`dev/shapes-rev-cov.R`). Two designs, two seed streams, one
conclusion: the `joint` interval is NARROW at a known level of a
mixed fit, and the term that is missing is the conditional variance
of the group mode. Adding exactly `sigma^2 tau^2 / (sigma^2 + n_g
tau^2)` moves the coverage by +0.0164 here and +0.0141 there; it
closes the gap on the reviewer's design and overshoots on this one.
The attribution is therefore MEASURED, not asserted, and its
magnitude is right within a design's worth of noise.

**This section is not closed by "the binomial interval contains
0.95".** It does, on this design, and on the reviewer's it does not.
A pooled binomial interval assumes the 2200 predictions are
independent and they are not: they share 220 fits, so it is not the
test. That is why the replicate-level standard error is printed beside
it and is the one quoted above.
What is settled is the DIRECTION and the SIZE of the term, and that
`?predict.frmtmb_fit` names it.

**CORRECTION (punch round 2).** The sentence that stood here gave, as
the reason to leave `Var(b | y)` out, that carrying it "would make it
a marginal predictive interval rather than one conditional on the
modes, which is a different quantity from the one every other method
in this package reports, and from brms's". The brms half is WRONG.
brms's predictive draws at a level it SAW carry the posterior of that
level's effect, `r_g`, draw by draw: its interval at a known level
includes the group-effect uncertainty, and `Var(b | y)` is the
frequentist analogue of exactly that term. So leaving the term out is
a divergence FROM brms, not a way of matching it, and the mixed
shortfall above (0.9445 here, 0.9336 on the reviewer's design) is the
size of the divergence. What remains true is the package half: every
other method here (`fitted()`, `frm_linpred()`, `residuals()`) is
conditional on the modes, and a `predict()` that carried
`Var(b | y)` would be the one method that is not. Behavior is
unchanged in this round, as the coordinator directed; the decision,
with both halves stated, is filed under `## Open - high priority` of
`dev/test-backlog.md`.

## 3a. What predict() refuses, and why each one

| argument | why |
|---|---|
| `type`, `se.fit`, `dpar`, `scale` | the old vocabulary; each message names `frm_linpred()` or `fitted()` |
| `sample_new_levels = "uncertainty"`, `"old_levels"` | both resample the POSTERIOR draws of the levels the fit saw, and there are none. `"gaussian"` is ANSWERED and is what happens |
| `draw_ids` | there are no stored draws to index; `ndraws` sets how many are taken |
| `cores` | the simulation is not parallelized |
| `sort = TRUE` | rows come back in the order of the data, always |
| `negative_rt` | brms's sign convention for its `wiener` family |
| `re.form`, `allow.new.levels` | lme4's spellings, retired at 0.57.0 |

`ndraws`, `summary`, `robust`, `probs`, `transform`, `ntrys` and
`sample_new_levels = "gaussian"` are ANSWERED, because `predict()` has
real draws: it simulates them.
`ntrys` reaches the family simulator's `max_iter` through a new
optional field on `sim_context()`, which is how a `trunc()`ed draw's
rejection limit is set.

## 4. Where the glmmTMB vocabulary went, and every caller that moved

`frm_linpred()` is `predict.frmtmb_fit()` of 0.60.0 with its body
unchanged and its arguments unchanged: `type` (`link`, `response`,
`conditional`, `zprob`, `zlink`, `disp`), `dpar`, `resp`, `re_formula`,
`se.fit`, `allow_new_levels`. `predict()` refuses `type`, `se.fit`,
`dpar` and `scale` BY NAME and each message says where that argument
went.

**Internal callers, all moved.** Every call site that asked
`predict()` for a linear predictor now asks `frm_linpred()`:

| file | what reads it |
|---|---|
| `R/conditional-effects.R` | `ce_boot_one()`, the categorical grid, the mean-display band, the link band (5 sites) |
| `R/interop.R` | `get_predict.frmtmb_fit()`, which is marginaleffects' seam |
| `R/spectral.R` | `frm_series_draw()` |
| `extensions/frmtmb.coupling/R/coupling.R` | `frm_coherence()` and `frm_phase()` (2 sites) |
| `extensions/frmtmb.eam/R/ddm-shared.R` | `ndt_time()` |
| `extensions/frmtmb.spline/R/curve.R`, `curve-cov.R`, `rp-check.R` | `frm_curve()`'s transform check, `sp_predict_eta()`, the covariance check, the Royston-Parmar monotonicity check (4 sites) |
| `extensions/frmtmb.sample/R/methods-draws.R` | `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`, one per draw (3 sites) |
| `extensions/frmtmb.sample/R/conditional-effects-draws.R` | the category names of a polytomous grid |

`fitted()` and `residuals()` gained `fitted_point()` and
`residual_values()`, the plain-vector helpers the public methods now
wrap, and `plot.frmtmb_fit()` and `dharma_residuals()` read those: a
scatter plot needs one number per row, not a four-column summary.

**The covariance seam.** `vcov()` is brms's population-level block, and
five internal callers need the covariance of EVERY estimated
coefficient instead, because each is paired with a coefficient vector
that has the extra rows. They read the new `vcov_estimated()`:
`frm_multiple()`'s pooling (twice), `emmeans`'s `emm_basis()`, and
`frmtmb.latent`'s HMM start check. A covariance one row short of the
coefficient vector is the WRONG matrix, not a differently named one;
section 5 has what it did to `avg_slopes()` before this.

`vcov()` itself is UNCHANGED on every fit that has nothing outside
the coefficients, and that is asserted rather than assumed: on a
gaussian with a factor, a bernoulli, a mixed and a distributional fit,
`vcov(fit)` is `identical()` to the exact subset of
`vcov_estimated(fit)` its rows name, so the reviewer's bitwise result
still holds after the delta path was added beside it.

Punch round 1 split that seam once more. `marginaleffects::get_coef()`,
`set_coef()` and `get_vcov()`, `insight::get_parameters()` and
`get_varcov()`, and `influence()`'s `dfbetas()` and
`cooks.distance()` read `interop_coef_vector()` and `interop_vcov()`
instead, which are `vcov_estimated()` plus an ordinal fit's
thresholds and `cs()` coefficients. On every other fit they ARE
`vcov_estimated()`, returned unchanged; on an ordinal fit they are
the joint matrix `hypothesis()` assembles, cut to those rows. Section
5a has the factor of 69 this closed.

## 5. Interop, measured by RUNNING it

`dev/shapes-interop.R` RUNS emmeans, insight and marginaleffects on
four fits (gaussian with a factor, bernoulli, gaussian mixed, and
ordinal `cumulative()`) and prints, for EVERY accessor, its class, its
length, its dim and its dimnames before any value:
`emmeans()` at a grid, `emmeans(pairwise ~ f)`,
`insight::get_predicted()`, `get_parameters()`, `find_parameters()`,
`get_varcov()`, `model_info()`, `get_residuals()`, `get_response()`,
`stats::fitted()`, `stats::residuals()`, `coef()`, `vcov()`,
`lmtest::coeftest()`, `marginaleffects::avg_slopes()` and
`marginaleffects::predictions()` over a `datagrid()`, plus an outside
judge: the same ordinal average marginal effects from `MASS::polr`.
The same script runs against the base build `rellib-r3` and against the
lane, and the two outputs are diffed:
`dev/shapes-log/interop-base.txt` against
`dev/shapes-log/interop-lane.txt`.

**CORRECTION (punch round 1).** The sentence that stood here, "Every
number agrees with the base build ... Nothing in the diff is a changed
value", was FALSE, and the instrument that produced it was the reason.
It printed

    round(head(as.numeric(insight::get_residuals(f)), 3), 4)

and `as.numeric()` flattens a matrix COLUMN FIRST. When
`get_residuals()` began returning the whole 150 x 4 summary matrix,
because `residuals()` is brms's matrix now and insight's default hands
back whatever `residuals()` gives, the first three numbers were still
the first three residuals, bit for bit, and the diff was empty. A value
probe cannot see a shape change. `mean()`, `sd()`, `sum(r^2)` and
`qqnorm()` on the result all changed value in silence.

The script is now a HARNESS and not a log. Every accessor carries a
contract that is checked, and the script exits 1 and prints
`HARNESS FAILED` when one breaks. Run against the build the review read
it reports **13 failures**
(`dev/shapes-log/interop-lane-seen-failing.txt`): `get_residuals()` not
a length-n vector on four fits, `get_residuals()` without insight's
class on four fits, `coeftest()` losing a `vcov()` row and `coef()` and
`vcov()` naming different parameters on two fits, and the ordinal
`avg_slopes()` standard error disagreeing with `MASS::polr`. Against
the fixed build it reports **0**. The base build reports 7, which is
where the lane is now strictly better than it: base's
`get_residuals()` carries no `insight_residuals` class, base's
`vcov()` has a row `coef()` does not name, and base has the same
ordinal defect.

The three warnings the old text pointed at are still gone, and for the
reason it gave: twice insight's "Could not compute standard errors or
confidence intervals because the model and variance-covariance
matrices are non-conformable", and once "Could not apply Delta method
to transform standard errors". Both were the old `get_predicted()`
fallback, which had no method for this class and no standard errors;
`get_predicted.frmtmb_fit()` computes them.

**It did not start that way, and the first run is the reason two seams
moved.** Against the lane before those two fixes:

- `insight::get_predicted()` returned NOTHING on all four fits, with
  the warning "Could not compute predictions for model of class
  `frmtmb_fit`". insight's default calls `predict(x, type = )` and then
  `fitted(x, type = )`, neither of which takes `type` any more, and its
  default WARNS and returns `NULL` rather than erroring. Every
  prediction was silently lost.
- `marginaleffects::avg_slopes()` and `predictions()` returned
  DIFFERENT standard errors, and some `NA`: on the gaussian fit the
  first row moved from 0.1597 to 0.0804 and row 5 became `NA`.
  `get_vcov()` returned `vcov()`, which is brms's population-level
  block, while `get_coef()` returns every estimated coefficient, so the
  numeric jacobian and the covariance had different lengths. That is a
  wrong matrix, not a renamed one, and it is the strongest argument for
  keeping `vcov_estimated()` as a separate seam.

## 5a. The ordinal hole, and the outside judge that closes it

An ordinal fit keeps its THRESHOLDS and its `cs()` coefficients in
`extra_names` components of the parameter template, not in `beta` or
`betad`, so the coefficient machinery never saw them. Three surfaces
were short because of it, and one of them gave a wrong number:

- `fixef()` reported `x` where brms reports `Intercept[1]`,
  `Intercept[2]`, `x`; `vcov()` was 1 x 1 where brms is 3 x 3;
  `print(fit)` showed a one-row coefficient block. Measured against
  brms on the same data, `dev/shapes-p1-brmsord.R`, which also settles
  where a `cs()` term goes: brms's rows on `bf(ord ~ x + cs(z)) +
  sratio()` are `Intercept[1]`, `Intercept[2]`, `x`, `z[1]`, `z[2]`,
  and its `variables()` names them `b_Intercept[k]` and `bcs_z[k]`.
- `marginaleffects::avg_slopes()` reported **0.00029** for the middle
  category where `MASS::polr` reports **0.02016**, a factor of 69, on
  a point estimate that agreed to five decimals. marginaleffects
  builds its Jacobian by perturbing `get_coef()` one entry at a time,
  and the thresholds were not in that vector, so they contributed
  nothing to the standard error. This one PREDATES the shapes work:
  the base build gives the same 0.00029.

The estimates and their standard errors are judged from OUTSIDE, by a
package whose numbers do not come through frmtmb
(`dev/shapes-p1-ordjudge.R`, `dev/shapes-log/ordjudge.txt`,
`bf(ord ~ x) + cumulative()` against `MASS::polr`, n = 120, seed
20260917):

| row | frmtmb estimate | polr | frmtmb SE | polr SE | SE ratio |
|---|---|---|---|---|---|
| `Intercept[1]` | -2.3894838 | -2.389310 | 0.3296123 | 0.3295998 | 1.000038 |
| `Intercept[2]` | -0.2167022 | -0.216763 | 0.2282469 | 0.2282429 | 1.000017 |
| `x` | 1.7018563 | 1.701719 | 0.2780606 | 0.2780486 | 1.000043 |

Maximum absolute estimate difference 1.7e-04, maximum relative
standard-error difference 4.3e-05. The threshold standard error is a
DELTA METHOD, not a subset: `cumulative()` and `sratio()` estimate
`(tau_1, log increments)` so that the thresholds cannot cross, and the
derivative is taken numerically from the map the FAMILY declares in
`post$ord_thresholds`, so an ordinal family shipped by another package
is covered.

The interop seam agrees with `polr` too, on all three categories,
which is the check the first round's harness could not make because it
printed frmtmb's number with nothing beside it
(`dev/shapes-interop.R`, now a contract):

| category | frmtmb `avg_slopes()` SE | `MASS::polr` | ratio |
|---|---|---|---|
| 1 | 0.02414848 | 0.02414857 | 1.0000 |
| 2 | 0.01934503 | 0.01934454 | 1.0000 |
| 3 | 0.02378709 | 0.02378688 | 1.0000 |

A cluster-robust `vcov(cluster = )` cannot cover the threshold rows,
because there is no joint cluster-robust matrix behind them. It
reports the coefficient rows it does have and WARNS, naming
`vcov(object)` as the route to the thresholds, rather than returning a
matrix whose dimension is the only clue.

## 6. The ported tier

The acceptance is `dev/brmsport-ledger.tsv` and the generated files it
writes. The pipeline is lane `wt-brmsport`'s, unchanged in what it
does, and re-run here end to end:

```
sh dev/brmsport-record.sh          # 15 files, one R process each
Rscript dev/shapes-verdicts.R      # the verdict moves, idempotent
Rscript dev/brmsport-ledger.R      # refuses a stale verdict; writes the ledger
Rscript dev/brmsport-gen.R         # regenerates the 14 files
sh dev/brmsport-tier.sh            # runs the tier with the verdicts asserted
```

**A fixed defect is a DELETED verdict row, never a typed "pass".**
`dev/brmsport-ledger.R` takes a pass from the RUN and stops if a manual
verdict claims otherwise, so the 40 rows below were found by running,
not by editing: the first ledger pass refused 41 problems, 40 of them
"HOLDS but has manual verdict", and `dev/shapes-verdicts.R` drops
exactly those 40. The 41st is the one that went the other way, `:326`,
and the pass count therefore rises by 39 rather than 40.

The runner scripts were tied to the `wt-brmsport` worktree and to the
StanHeaders 2.32.10 pin, which is gone. `dev/brmsport-run.R`
now takes its library from `FRMTMB_PORT_LIB`, defaults it to the lane
library, resolves the Stan cache from `getwd()`, and asserts
`tmbstan >= 1.2.1` and the `-std=gnu++17` flag in the user Makevars
instead of the pin; `dev/brmsport-record.sh` and `dev/brmsport-tier.sh`
refuse to run outside a worktree root rather than `cd`-ing to a fixed
one.

### Totals, before and after

| outcome | before (`wt-brmsport`) | after |
|---|---|---|
| pass | 192 | **231** |
| defect | 112 | **83** |
| divergence | 34 | **35** |
| pending 2.6d | 12 | **0** |
| cannot transfer | 144 | **145** |
| total | 494 | 494 |

Bin 1 passes move from 192 of 494 (38.9 percent) to 231 of 494 (46.8
percent). The `pending 2.6d` class is gone: 9 of its 12 rows are passes,
and 3 are defects for a reason that is not the shape. After punch round
1 those three are `:772` (`sample_new_levels = "old_levels"` is still
refused, and correctly), `:775` (the call runs now and the FIXTURE's
`newdata` has 2 rows, the `fit$data` partial match `:764` records) and
`:747` (still dies on `fit1$data` for the same reason).

frmtmb.sample's half moves from 88 of 165 (53.3 percent) to 94 of 165
(57.0 percent).

**The 40 rows whose verdict was dropped**, by what fixed them. All 40
are passes now; the net gain is 39 because `:326` went the other way,
from pass to `cannot transfer`:

| rows | what |
|---|---|
| `:292 :293 :325 :328 :334 :337 :339 :340 :342 :348 :355` | `fitted()` is brms's four-column summary, the `n x 4 x K` array on an ordinal fit, and takes `nlpar` and `allow_new_levels` |
| `:364 :369` | `fixef()` is brms's matrix, in brms's row ORDER, and takes `pars` |
| `:584 :585` | `ngrps()` is brms's named list |
| `:784 :887 :888 :894 :896 :897` | `print(fit)` is brms's layout; `summary()` carries `$fixed` and `$random` and honors `priors = TRUE` |
| `:825 :835` | `residuals()` is brms's four-column summary and takes `probs` |
| `:999 :1000` | `vcov()` is brms's 9 x 9 block and takes `correlation` |
| `:726 :727 :729 :750 :753 :758 :761 :762 :767` | `predict()` is brms's predictive summary (item 2.6d) |
| `:593 :594 :633 :634` | `nsamples()` and `posterior_samples()` answer as brms's do |
| `:713 :714` | `point_estimate` and `ndraws_point_estimate` are honored |

**The rows that changed verdict without becoming passes** are in
`dev/shapes-verdicts.R` with their reasons: `:891` (Rhat and the two
ESS columns are MCMC diagnostics) becomes a divergence; `:326` becomes
`cannot transfer / fixture` because fixture 1 does not converge, so
`fitted()`'s new `Est.Error` column is `NA` and `all(fi > 0)` is `NA`
where the `Estimate` column is positive; `:595`, `:635`, `:772`,
`:775`, `:314`, `:317` and `:747` keep their outcome and gain the
reason the new behavior gives them.

### What punch round 1 did to the tier

**Nothing, in the totals**: 231 pass, 83 defect, 35 divergence, 145
cannot transfer, before and after. That is a finding and not a
non-result. The five majors and two blockers the review found sit on
surfaces the ported tier does not reach: `dev/brmsport-ledger.tsv` has
no `fixef()` or `vcov()` row on an ordinal fixture, no
`insight::get_residuals()` row, no `coef()`/`vcov()` pairing row and
no multivariate `predict()` row. The tier is the acceptance for the
brms API, and it says nothing about the seams a third-party package
reaches, which is why `dev/shapes-interop.R` had to become a harness
with contracts of its own.

Three rows moved WITHIN the ledger:

- `brmsfit-methods:775` was a defect because `predict()` refused
  `sample_new_levels`. That refusal is gone, the call runs, and the
  row is STILL a defect, for the fixture's own pre-existing reason:
  its `newdata` is `fit5$data[1:5, ]`, `fit$data` partial-matches
  `fit$data2` (the fit has no `data` element), so `newdata` has 2 rows
  and the answer is 2 x 4 where the assertion wants 5 x 4. That is
  exactly what `:764` already records. The verdict's class moved from
  `argument` to `fit-data`.
- `brmsfit-methods:772` keeps its defect, with the reason it actually
  has now: `"old_levels"` resamples the posterior draws of the levels
  the fit saw, and a maximum-likelihood fit has none.
- `brmsfit-methods:340` gained a NOTE and a new ledger class,
  `weak-pass`. It holds, so no verdict can be given for it, but it is
  WEAKER than the assertion brms wrote: `fitted()` is brms's n x 4
  matrix now, so `fi[1, ]` and `fi[2, ]` each carry `Estimate` plus
  three `NA` cells on this non-converging fixture and
  `expect_equal(fi[1, ], fi[2, ])` compares `c(19.122, NA, NA, NA)`
  with itself. Before the shapes item it compared two numbers and
  errored on `fi[1, ]`. It moved from defect to pass on that
  weakening, not on a fix. `dev/brmsport-notes.tsv` carries the note
  and `dev/brmsport-ledger.R` reads it, so the generated tier file
  prints it beside the assertion.

**The mixture fixture found a real defect in `predict()`'s own
construction**, which only the new `sample_new_levels` path could
reach. `predict(fit5, newdata, allow_new_levels = TRUE)` died with
"NA in probability vector", a message about the family simulator's
arguments rather than about the draw. The cause: fixture 5's
`(1 | patient)` log standard deviation has variance **7.97e6** in
`vcov(fit5, full = TRUE)`, so a draw of 1365 makes the unseen level's
variance `exp(2730)` and `mu1` is `Inf`. A replicate whose
distributional parameters are not finite is now DROPPED, counted, and
reported in a warning that names the cause and points at
`vcov(object, full = TRUE)` and `param_uncertainty = FALSE`; the
summary is taken over the replicates that are there. Every replicate
non-finite is an error, not an empty answer.

### The ledger's own generated summary

Pasted verbatim from `dev/brmsport-log/ledger-summary.md`:

<!-- BEGIN GENERATED: dev/brmsport-ledger.R -->

### Per file

| file | assertions | pass | defect | divergence | cannot transfer |
|---|---|---|---|---|---|
| `tests.brm.R` | 23 | 16 | 5 | 0 | 2 |
| `tests.brmsfit-helpers.R` | 3 | 0 | 0 | 0 | 3 |
| `tests.brmsfit-methods.R` | 223 | 100 | 63 | 30 | 30 |
| `tests.brmsformula.R` | 16 | 7 | 0 | 0 | 9 |
| `tests.brmsterms.R` | 5 | 0 | 0 | 0 | 5 |
| `tests.data-helpers.R` | 6 | 0 | 3 | 0 | 3 |
| `tests.emmeans.R` | 11 | 4 | 2 | 0 | 5 |
| `tests.families.R` | 84 | 45 | 1 | 1 | 37 |
| `tests.priors.R` | 36 | 17 | 4 | 4 | 11 |
| `tests.standata.R` | 87 | 42 | 5 | 0 | 40 |
| **total** | **494** | **231** | **83** | **35** | **145** |

### Outcome by class

| outcome | class | assertions |
|---|---|---|
| pass | - | 203 |
| pass | own-words | 27 |
| pass | weak-pass | 1 |
| defect | accepts-refused | 3 |
| defect | argument | 24 |
| defect | different-error | 1 |
| defect | filed | 1 |
| defect | fit-data | 14 |
| defect | internal-error | 5 |
| defect | misparse | 1 |
| defect | naming | 3 |
| defect | output | 3 |
| defect | refuses-accepted | 22 |
| defect | silent | 2 |
| defect | spelling | 4 |
| divergence | class | 18 |
| divergence | hollow | 1 |
| divergence | no-draws | 9 |
| divergence | policy | 7 |
| cannot transfer | absent | 106 |
| cannot transfer | brms-internal | 10 |
| cannot transfer | fixture | 2 |
| cannot transfer | mcmc | 4 |
| cannot transfer | no-draws | 1 |
| cannot transfer | stan | 22 |

Bin 1 passes: 231 of 494 (46.8%).
Against bins 1 and 2: 231 of 823 (28.1%); bin 2 was not ported.
Runs testthat alone would count as a pass and the harness does not (missing function or object, stale object, argument-name refusal, a NULL read through a partial $ match), over both packages: 28; hollow passes marked by hand: 1.

### The frmtmb.sample half

| file (tier) | assertions | cannot transfer | defect | divergence | pass |
|---|---|---|---|---|---|
| `tests.brmsfit-methods.R (sample)` | 61 | 13 | 6 | 0 | 42 |
| `tests.brmsformula.R (both)` | 16 | 9 | 0 | 0 | 7 |
| `tests.brmsterms.R (both)` | 4 | 4 | 0 | 0 | 0 |
| `tests.families.R (both)` | 84 | 37 | 1 | 1 | 45 |
| **total** | **165** | **63** | **7** | **1** | **94** |

frmtmb.sample passes 94 of its 165 runs (57.0%): 42 of the 61 sample-tier and 52 of the 104 both-tier assertions.

<!-- END GENERATED -->

## 7. Suites and checks

**The whole ungated suite, one R process per file, all eight packages**
(`dev/shapes-suite.ps1`, which is `dev/release/run-suite.ps1` pointed
at the lane library; log `dev/shapes-log/punch1-final.log`, and
`dev/shapes-log/final.log` is round 0's):

```
RAN 265 of 265 files   (15737 expectations)
== frmtmb 154 files ==      == frmtmb.eam 26 ==    == frmtmb.sample 29 ==
== frmtmb.coupling 9 ==     == frmtmb.spline 14 == == frmtmb.learn 14 ==
== frmtmb.latent 9 ==       == frmtmb.ode 10 ==
```

**Zero files with a failure, an error or a load error**, after one
was fixed and its file re-run. The run itself reported ONE:
`frmtmb.spline/test-difference.R pass=58 fail=1`, at
`as.numeric(D %*% stats::coef(fit)[keep])` with
`keep = c("(Intercept)", "facB", "x", "facB:x")`. `coef()` takes
brms's names since punch round 1, so `[keep]` gave `NA` for the
intercept. The test reads `fixef_by_dpar(fit)$mu[keep]` now, which is
the design-column vocabulary `keep` is written in, and passes 59 of
59. That and the frmtmb.sample `brmshypothesis` assertion the check
found are the only two call sites in the eight packages the class and
naming changes reached.

An earlier run of the same 264 files had ONE,
`frmtmb.sample/test-draws-spellings.R pass=86 fail=1`, and it is worth
naming because only that guard would have caught it. The guard
compares each of ten brms-facing methods against brms's own POSITIONAL
signature; adding `point_estimate` to `posterior_predict()` had pushed
brms's `ntrys` out of position 11, so a positional brms call would have
landed `ntrys` on `point_estimate` and been ANSWERED, silently. `ntrys`
and `cores` are now in brms's slots and refused by name.

**The gated ported tier, verdicts asserted**
(`sh dev/brmsport-tier.sh`, log
`dev/shapes-log/tier-punch1.txt`), rerun after punch round 1 with
the regenerated ledger and verdicts; guards 60 of 60:

```
frmtmb       brm 23   brmsfit-helpers 3   brmsformula 16   brmsterms 5
             data-helpers 6   emmeans 11   families 84   methods 162
             priors 36   standata 87
frmtmb.sample brmsformula 16   brmsterms 4   families 84
             helper-copy 33   methods 61
TIER ran 15 of 15 files
```

631 expectations over 15 files; 0 failures, 0 errors, 0 skips.

**`R CMD check --as-cran`**, built WITH vignettes and WITHOUT
`--no-manual`, on the two packages whose API this lane changed
(`dev/shapes-check.ps1`, log `dev/shapes-log/check.log`). The base
build's recorded statuses, for comparison, are in
`dev/release/check.log`: frmtmb `1 NOTE`, frmtmb.sample `OK`.

```
frmtmb          Status: 1 NOTE
frmtmb.sample   Status: OK
```

Each package's run is one `R CMD build` plus one `R CMD check`, and
the driver clears its log per invocation, so
`dev/shapes-log/check.log` is the LAST run of the pair. Punch round 1
ran it twice: the first run is kept as
`dev/shapes-log/check-punch1-run1.log` and is where the two findings
below were measured (`frmtmb` `2 NOTEs`, `frmtmb.sample` `1 ERROR`);
the second is the one quoted above, after the knitr leftover was
removed and the extension test was fixed.

The one NOTE is `checking HTML version of manual ... NOTE`,
"Skipping checking math rendering: package 'V8' unavailable", which is
this machine's and is what the base build reports for frmtmb too.
Neither package has a WARNING or an ERROR, and neither status differs
from the base build's.

**Two real findings came out of the first check and are fixed.** They
are recorded because neither would have been found by the suite:

- **Three vignettes failed to BUILD**, which failed `R CMD build` and
  left frmtmb with no tarball to check at all: `case-studies.Rmd` at
  `fixef(frm(bf(ls ~ mo(income) * z) + gaussian(), data = dmo))$mu`,
  `habit.Rmd` at `lapply(fixef(f), unname)` inside a bootstrap
  statistic, and `spectral.Rmd` at
  `fixef(frm(...))[["mu"]][[2]]`. Each is a `$` or `[[` read of the
  per-dpar list on the new matrix, and each was missed by the
  mechanical sweep for the same reason: the sweep's pattern requires
  the object to be one token, and these three wrap a whole `frm()`
  call. `case-studies.Rmd`'s second site was never reached by the
  build, because a vignette stops at its first error.
- **frmtmb.sample checked `1 WARNING, 1 NOTE` where the base build is
  `OK`**: `arg_desc()` is a frmtmb internal that frmtmb.sample does not
  import, and two Rd pages documented arguments no function on them has
  any more. Both are fixed.

**A third came out of punch round 1's check, and it was not a code
defect at all.** frmtmb checked `2 NOTEs`, the second of them:

```
* checking files in 'vignettes' ... NOTE
The following directory looks like a leftover from 'knitr':
  'figure'
```

`vignettes/figure/` holds 16 PNGs that a knitr run in the worktree
left behind: the review's `dev/shapes-rev-vig.R` knitted all 21
vignettes in the tree. `.gitignore` already has `figure/`, so git
never showed it, and `.Rbuildignore` did NOT, so `R CMD build` copied
it into the tarball. The directory is deleted, and
`^vignettes/figure$` is in `.Rbuildignore` now, so a knit in the tree
cannot put the NOTE back. **Six of the seven extensions have the same
leftover from the same run** (`frmtmb.coupling`, `frmtmb.eam`,
`frmtmb.latent`, `frmtmb.learn`, `frmtmb.ode`, `frmtmb.spline`); those
directories are deleted too, but their `.Rbuildignore` files are
another lane's and were not edited, so the same NOTE will appear the
next time anyone knits in one of those trees and then builds it. Worth
one line in each at release time. `frmtmb.sample` never had one, which
is why its check was clean before and after.

**A fourth came out of the same check and was a real miss.**
frmtmb.sample checked `1 ERROR`:

```
Failure ('test-reparam.R:543:3'):
  Expected hypothesis(ds, "sd_g__Intercept > 0", class = NULL)
  to inherit from "brmshypothesis".
  Actual class: "frmtmb_hypothesis".
[ FAIL 1 | WARN 10 | SKIP 11 | PASS 1711 ]
```

Dropping the `brmshypothesis` class was applied to the four core tests
that asserted it and to none in the EXTENSIONS, because the sweep that
found them read `tests/testthat/` and not
`extensions/*/tests/testthat/`. The file now asserts the class frmtmb
gives, asserts that brms's is NOT on it, and asserts the five element
names, so it guards the shape as well as the class.
`extensions/frmtmb.sample/vignettes/brms-posterior.Rmd` said the
generics "return brms's objects ... and brms's `brmshypothesis` list";
it says brms's SHAPES under frmtmb's own class now. The 10 WARNs are
testthat's count of the deprecation warnings the draws accessors are
supposed to raise and are not new.

This is what `R CMD check` is for on a lane that changes a class: the
one-per-file suite run had not reached that extension when the check
did.

### Files, by what they are

**New in core:** `R/brms-shapes.R` (the summary-matrix machinery, the
brms row order and the outer-parameter seam) and `R/predict-brms.R`
(brms's `predict()` and the simulator behind it). **New in
frmtmb.sample:** the three summarizing methods, the three deprecated
accessors and `draws_at_point_estimate()`, all in
`R/methods-draws.R`. **New exports:** `frm_linpred()`,
`fixef_by_dpar()` and, as extension API, `vcov_estimated()`,
`fam_is_category_valued()`, `predict_category_props()`,
`brms_summary_matrix()`, `brms_summary_array()`,
`brms_summarize_draws()` and `brms_prob_cols()`.

**New tests:** `tests/testthat/test-brms-shapes.R` (21 blocks),
`extensions/frmtmb.sample/tests/testthat/test-brms-shapes-draws.R`
(6 blocks) and, from punch round 1,
`tests/testthat/test-brms-shapes-punch.R` (14 blocks, 55
expectations), one per item the review found.

**What punch round 1 changed, file by file.**

| file | what |
|---|---|
| `R/insight.R` | `get_residuals.frmtmb_fit`, insight's `brmsfit` body; `get_varcov()` and `get_parameters()` read `interop_vcov()` and `interop_coef_vector()` |
| `R/interop.R` | `interop_coef_names()`, `interop_coef_vector()`, `interop_vcov()`; `get_coef()`, `set_coef()` and `get_vcov()` carry an ordinal fit's thresholds and `cs()` coefficients |
| `R/brms-shapes.R` | `brms_extra_fixef()`, `brms_fixef_values()`, `brms_fixef_extra_vcov()`, `fd_gradient_row()`; `brms_fixef_rows()` returns the extra rows too; `brms_summarize_draws()` is NA-aware and leaves row dimnames NULL |
| `R/brms-names.R` | `brms_lp_rows()` and `brms_lp_coef_names()`, plus the `n_beta`/`keep_d` attributes the map needs |
| `R/methods-fit.R` | `coef()` takes brms's names and the threshold rows; `vcov_brms_block()` gains the delta path and the cluster-robust warning; `summary_fixed_frame()` carries the thresholds; `summary_empty_block()` and brms's empty shapes; `print()` tests rows rather than NULL |
| `R/predict-brms.R` | the multivariate array, the unseen-level draw (`predict_new_level_spec()`, `predict_new_level_draw()`, `mvn_draw_cov()`), the non-finite-replicate guard (`dpv_all_finite()`), `sample_new_levels = "gaussian"`, `ndraws` checked before coercion |
| `R/predict.R` | ordinal `fitted()` quantiles clamped to [0, 1]; row dimnames NULL on `fitted()` and `residuals()`; the RE1.3 claim rewritten to what the code does |
| `R/confint.R` | the `brmshypothesis` class dropped |
| `R/influence.R` | `dfbetas()` and `cooks.distance()` read `interop_vcov()`, so they cannot recycle against a shorter vector |

**New and rebuilt dev scripts.** `dev/shapes-interop.R` is a HARNESS
now, not a log: class, length, dim and dimnames for every accessor,
contracts that are checked, `HARNESS FAILED` and exit 1 when one
breaks. `dev/shapes-p1-brmsord.R` fits the ordinal fixtures in brms so
that where the thresholds and `cs()` rows go is measured rather than
assumed. `dev/shapes-p1-ordjudge.R` judges them against `MASS::polr`.
`dev/shapes-coverage.R` gained the `newlevel` design, the `exact`
positive control and the `condvar` attribution test.
`dev/brmsport-notes.tsv` and the `weak-pass` class carry a note on a
pass that is weaker than brms wrote it.

**The test suite and the vignettes moved with the API**, mechanically
and then by hand. Counted from the tree as it stands: the test files
now hold 631 `frm_linpred(` calls, 641 `fixef_by_dpar(` calls and 261
`[, "Estimate"]` reads (or `[, "Estimate", ]` on an ordinal,
categorical or matrix-valued response), and the vignettes hold 74 and
64. Every sweep needed its reverts and its hand corrections, and those
are the interesting part:

- 10 `predict(` sites were renamed and put BACK, because the object was
  a reference fit from glmmTMB, lme4, mgcv, `nls()` or brms, where
  `predict()` still means what it always did. `frm_linpred()` is not
  generic, so each of those errored loudly rather than answering.
- 23 `glmmTMB::fixef()` and `lme4::fixef()` calls were renamed across
  the namespace qualifier and put back; those packages have no
  `fixef_by_dpar()`.
- A parse check of all 288 test files after each sweep caught 14 files
  where an insertion had landed INSIDE a string literal, and one
  regexp that had gained an unescaped `\(`.
- 122 lines that went past 80 columns were re-wrapped
  (`dev/shapes-reflow.py`, which breaks only at a comma outside every
  string and reports what it cannot take); 16 more were done by hand.

## 7a. The port harness, made lane-agnostic

`dev/brmsport-run.R`, `dev/brmsport-record.sh`, `dev/brmsport-tier.sh`,
`dev/brmsport-guards.R` and `dev/brmsport-globcheck.ps1` each hard-coded
the `wt-brmsport` worktree, and two of them asserted the StanHeaders
2.32.10 pin, which is gone. They now take the library from
`FRMTMB_PORT_LIB` (defaulting to this lane's), resolve the tree from
the working directory and refuse to run outside a worktree root, and
assert `tmbstan >= 1.2.1` plus the `-std=gnu++17` flag in the user
Makevars in place of the pin. A populated Stan cache hides a compile
failure, which is why those two assertions and not a green run are the
evidence.

Guards after the change: **60 of 60**
(`dev/brmsport-log/guards-shapes.txt`), including the helper copy
identical to core's while a mutated copy is not, the helper-copy test
passing on the real tree (33 expectations) and failing on the mutant,
the 27 own-words patterns each matching their own row's message and
none of the other 26, and the gate skipping the whole tier with
`FRMTMB_BRMS_FIT_TESTS` unset. The gated glob throws on an empty
directory and finds 15 jobs on the tree
(`dev/brmsport-globcheck.ps1`).

## 8. What was not done, and why

**Left as defects in the ported tier, with the reason on the row.**

- `fitted(ndraws =)`, `fitted(draw_ids =)`, `residuals(ndraws =)` and
  their relatives stay refused. There are no draws to thin, and item
  2.6f's own rule is that a draws-only argument is refused BY NAME. The
  rows that read a value PAST such a refusal (`brmsfit-methods:314`,
  `:317`) therefore stay defects: the setup line dies, so the
  assertion reads a stale object and the harness rejects it.
- `predict(sample_new_levels = "gaussian")` is ANSWERED since punch
  round 1, and the other two values are refused. An unseen level's
  effect is drawn from the block's own estimated covariance, once per
  replicate, at that replicate's parameters, which is exactly what
  brms's `"gaussian"` does. `"uncertainty"` and `"old_levels"`
  resample the posterior draws of the levels that WERE seen, and a
  maximum-likelihood fit has none of those. The design pieces come
  from `lp_extra_var()` and `extra_var_blocks()`, the same two
  functions `frm_linpred(se.fit = TRUE)` already used for the same
  term, so the variance the draw uses and the variance the standard
  error reports cannot drift apart.
- `residuals(newdata =)` stays refused, and `residuals()` on a
  MULTIVARIATE fit stays refused. Both are feature gaps that predate
  this item and neither is a shape.

  **CORRECTION (punch round 1).** The sentence that stood here said
  multivariate `fitted()` stays refused as well. That is not what the
  code does: `fitted(mv)` refuses with "not supported yet for
  multivariate fits" and `fitted(mv, resp = "y2")` ANSWERS with a
  150 x 4 matrix, in the base build too. The claim was wrong either
  way. `predict()` no longer has the gap in either direction: it
  returns brms's `nrow x 4 x nresp` array with the responses named,
  and `resp =` is a filter on that.

  That leaves an ASYMMETRY inside the package, stated rather than
  hidden: `predict(mv)` answers for every response and `fitted(mv)`
  still refuses while `fitted(mv, resp = )` answers. Closing it means
  stacking `fitted_point()` and `fitted_point_se()` over the
  responses, and the second of those is the joint delta method, which
  is where a multivariate fit's cross-response covariance would have
  to be decided. That is a feature question, not a shape one, and it
  was not in this item's brief.
- `summary()$fixed` has no `Rhat`, `Bulk_ESS` or `Tail_ESS`. Those are
  MCMC diagnostics and a maximum-likelihood fit has none, so
  `brmsfit-methods:891` moves from defect to a cited divergence rather
  than to a pass. The four brms columns come first and then the Wald
  test this package reports, which brms has no counterpart for.
- `residuals()` on an ordinal fit still answers, where brms refuses
  ("Predictive errors are not defined for ordinal models"). That is
  defect S7 of `dev/brmsport-findings.md` and it is a BEHAVIOR question,
  not a shape one; it was not in this lane's brief and the row stays a
  defect.

**Deliberately not changed.**

- `residuals()` keeps `type` in its second positional slot, where brms
  has `newdata`. Moving it would turn `residuals(fit, "pearson")`, which
  this package has always accepted, into a `newdata` that is not a data
  frame. brms's `"ordinary"` is accepted as a spelling of
  `"response"`, so a ported call reads the same.
- `frm_linpred()` keeps the response's own LEVEL names on a category
  matrix, where `fitted()` takes brms's `P(Y = k)`. `fitted()` is
  brms's method and `frm_linpred()` is this package's own.
- `variables()` keeps frmtmb's ORDER. brms lists `b_Intercept`,
  `b_sigma_Intercept`, then the non-intercepts; frmtmb lists each
  predictor's coefficients together. `fixef()`, `vcov()` and
  `summary()$fixed` take brms's order, because a ported script indexes
  those by position. Reordering `variables()` is a separate decision
  and would move `posterior_samples(pars = "^b_")`'s column order,
  which is why `brmsfit-methods:635` stays a defect.
- `hypothesis()` no longer carries `brmshypothesis`. The stated
  reason for keeping it, that removing it changes which `print()` and
  `plot()` methods dispatch, was MEASURED FALSE: frmtmb exports both
  generics for `frmtmb_hypothesis`, which is first in the class
  vector, so stripping brms's class leaves the printed output
  identical and `plot()` working
  (`dev/reviews/20260918-shapes.md`). Keeping it also contradicted the
  eighteen rule-2 divergences the ledger already cites. The class is
  `"frmtmb_hypothesis"` alone, and `is(x, "brmshypothesis")` in a
  ported script is a rule-2 divergence like the other eighteen.

**A cost that is real and is not a defect.** `fitted()` and
`residuals()` now compute a standard error, so both trigger
`sdreport()` where they used to be free. brms's own `fitted()` is not
free either, and `frm_linpred(type = "response")` is the free route to
the same estimates. The lazy-sdreport test still passes because
`fixef_by_dpar()` does not trigger one.

**CORRECTION (punch round 1).** Counted rather than timed
(`dev/shapes-rev-cost.R`, `-cost2.R`, which count
`autoscale_sdreport()` calls; control: `confint()` on a fresh fit is 1
and on the same fit again is 0, in both builds), the change is wider
than the paragraph above said. **`fixef()` costs 1 sdreport where base
cost 0, and `print(fit)` costs 1 where base cost 0**; `print()` is the
most common interactive call in the package and it went from free to
one. On a large mixed fit (n = 4000, 200 groups) `fitted()` on a FRESH
fit costs 2 where base cost 0, 0.58 s against 0.20 s, and 0 on the
second call. What stays free is what matters for the refitting paths:
`plot()`, `fixef(flatten = TRUE)` and `frm_allfit()` are all 0 in both
builds, and the ordinal finite-difference delta method costs 1
sdreport in total, not one per perturbed parameter.

## 8a. Anything that needs the user

Nothing blocks this item. Both questions this section put to the user
were ANSWERED in punch round 1 and are recorded here as decisions:

1. **`variables()` keeps frmtmb's ORDER.** brms lists every
   predictor's intercept first (`b_Intercept`, `b_sigma_Intercept`,
   then the rest); frmtmb lists each predictor's coefficients
   together. `fixef()`, `vcov()` and `summary()$fixed` take brms's
   order because a ported script indexes them by position;
   `variables()` keeps its own, because reordering it moves
   `posterior_samples(pars = )`'s column order and the draws matrix's,
   and brms's `variables()` also carries `Intercept`, `lprior`,
   `lp__` and the `r_g[...]` rows, so it will differ in CONTENT
   whatever its order. `brmsfit-methods:635` is the one ledger row
   that turns on it and a ported script fixes it with one reindex.
   Revisit is DONE: the thresholds became `fixef()` rows in punch
   round 1 and `variables()` already listed them, in frmtmb's order.
2. **`hypothesis()` no longer returns a `brmshypothesis` class.**
   Dropped in punch round 1. The stated reason for keeping it, that
   dropping it changes which `print()` and `plot()` methods dispatch,
   was measured false: frmtmb exports both generics for
   `frmtmb_hypothesis`, which comes first, so the printed output is
   identical and `plot()` works
   (`dev/reviews/20260918-shapes.md`). What a ported script loses is
   `is(x, "brmshypothesis")`, which the ledger already records as a
   rule-2 divergence eighteen times.

One thing is FILED rather than fixed, because it is not this item's
and the fix belongs with another lane's:

- **`re_formula = ~(1 | nosuch)` is silently `NULL`.**
  `predict()`, `fitted()`, `frm_linpred()` and
  `conditional_effects()` all return the CONDITIONAL prediction, bit
  for bit what `re_formula = NULL` returns, with no warning, on a
  grouping factor the model never had. `check_re_form()` accepts any
  formula and `re_form_keeps()` only asks whether it is `NA` or `~0`.
  It is the same family as the partial-`re_formula` defect another
  lane filed and wants one fix for both. Recorded under
  `## Open - high priority` of `dev/test-backlog.md`. Pre-existing;
  this lane did not touch `check_re_form`.

## 9. Every test seen failing first

Two new files, each run against the code that did not have the change
before a line of it was written. The logs are in the tree.

**Core, `tests/testthat/test-brms-shapes.R` against frmtmb 0.60.0**
(`dev/shapes-log/seen-failing-core.txt`, `dev/shapes-run.R`): 19 of 21
blocks fail, `pass=7 fail=25 err=14`. The two that already pass are the
two that assert a refusal the base build already makes
(`fixef(summary = FALSE)`, `fitted(ndraws =)`); they are in the file so
that the refusals are not lost while the shapes move around them.

```
  [ERROR] fixef() is brms's summary matrix, in brms's row order
  [ERROR] fixef() takes brms's pars and probs
  [FAIL]  vcov() covers brms's population-level coefficients only
  [ERROR] vcov() takes brms's correlation and pars
  [FAIL]  ngrps() is brms's named list
  [ERROR] fitted() is brms's four-column summary
  [ERROR] fitted() on an ordinal fit is brms's n x 4 x K array
  [ERROR] fitted() takes nlpar and allow_new_levels
  [ERROR] residuals() is brms's four-column summary
  [FAIL]  summary() carries brms's fixed and random slots
  [ERROR] summary(priors = TRUE) prints the priors, as brms does
  [FAIL]  print() of a fit uses brms's section headings
  [ERROR] variables() lists an ordinal fit's thresholds and cs terms
  [ERROR] predict() is brms's predictive summary
  [ERROR] predict(summary = FALSE) gives the simulated draws
  [ERROR] predict() on an ordinal fit gives brms's P(Y = k) columns
  [FAIL]  predict() refuses the retired glmmTMB type vocabulary
  [ERROR] frm_linpred() is the old predict(), unchanged
  [ERROR] predict() covers at close to the nominal rate
RESULT test-brms-shapes.R pass=7 fail=25 err=14 skip=0
```

**frmtmb.sample,
`extensions/frmtmb.sample/tests/testthat/test-brms-shapes-draws.R`
against the shared base build `rellib-r3`**
(`dev/shapes-log/seen-failing-sample.txt`, `dev/shapes-run-base.R`,
which never installs anything): 6 of 6 blocks fail,
`pass=0 fail=9 err=5`.

```
  [FAIL]  fitted() on draws summarizes posterior_epred()
  [ERROR] predict() and residuals() on draws are brms's summaries
  [ERROR] point_estimate collapses the draws, as in brms
  [ERROR] nsamples() is brms's, with brms's deprecation warning
  [ERROR] posterior_samples() is brms's, with brms's warning
  [ERROR] parnames() is variables(), with brms's warning
RESULT test-brms-shapes-draws.R pass=0 fail=9 err=5 skip=0
```

The `[FAIL]` on `fitted()` there is the S4 defect itself: the base
build does not ERROR, it returns `NULL` through
`stats::fitted.default()`, which is why the defect was silent.

**Punch round 1, `tests/testthat/test-brms-shapes-punch.R` against
the build the review read** (`dev/shapes-log/seen-failing-punch1.txt`,
`dev/shapes-run1.R`). It stops at ten failures, which it reaches:

```
  [FAIL]  predict(allow_new_levels = TRUE) widens at an unseen level
          fresh 1.00 against sqrt(sigma^2 + tau^2) = 1.38
  [ERROR] predict() takes brms's sample_new_levels = "gaussian"
          predict() cannot honor `sample_new_levels`
  [FAIL]  coef() uses brms's names, so coef()/vcov() pair up
          names(coef(ordf)) "x" against
          rownames(vcov(ordf)) "Intercept[1]" "Intercept[2]" "x"
  [FAIL]  predict() on a multivariate fit is brms's n x 4 x nresp
          length(dim(pr)) 2 against 3; dim(pr)[2:3] 4 NA against 4 2
  [ERROR] predict() on a multivariate fit: dimnames(pr)[[3]] out of bounds
  [FAIL]  an ordinal fitted() cannot report a probability outside [0, 1]
          min -0.0055, max 1.0058
  [FAIL]  empty summary blocks take brms's empty shapes
          s$random is a list where brms gives NULL
  ... and 9 more (the reporter caps at ten)
```

The four items whose failure only a third-party package can show are
recorded by the interop harness instead
(`dev/shapes-log/interop-lane-seen-failing.txt`, 13 contract
failures): `insight::get_residuals()` returning a 150 x 4 matrix and
carrying no `insight_residuals` class on four fits,
`lmtest::coeftest()` losing a `vcov()` row on two fits, and the
ordinal `avg_slopes()` standard error against `MASS::polr`. The base
build reports 7 of the same 13, which is the pre-existing half.

**The ported tier is the third body of evidence** and is not a new
test: section 6 gives the rows that moved, each of which was a recorded
`defect` or `pending 2.6d` verdict measured on the base code by lane
`wt-brmsport` before this lane started.

## 10. Punch round 2

The recheck confirmed both round-1 blockers fixed and no standard error
moved on any non-ordinal fit, and found two new silent wrong answers
that the round-1 fixes introduced, three majors and five minors. Every
number below is from `dev/shapes-p2-measure.R`
(`dev/shapes-log/p2-measure.txt`, seed 20260921) unless another script
is named. The tests are `tests/testthat/test-brms-shapes-punch2.R`
(22 expectations) and two new blocks in
`extensions/frmtmb.sample/tests/testthat/test-brms-shapes-draws.R`;
both were run against the round-1 build first and recorded failing
(`dev/shapes-log/seen-failing-punch2.txt`, 14 failing expectations,
the reporter capping its detail at ten, and
`dev/shapes-log/seen-failing-punch2-sample.txt`, 5).

**BLOCKER A. Every unseen level shared one effect draw.** Round 1 drew
the unseen level's effect through `extra_var_blocks()`, which groups
the design rows of one independent draw by a key, and single-membership
terms carried no key: every unseen level of a block fell under one
constant key and loaded ONE draw per replicate. A per-row variance
cannot see that, so the coverage arm could not; `summary = FALSE` could,
and three distinct new levels correlated at 0.215 to 0.221 in the
recheck, the same-level value. The fix is at the source: `pred_design()`
puts the level label on each random-effect part as `new_key`, which is
what the multi-membership parts already did. Measured with
`param_uncertainty = FALSE` and 4000 draws: three distinct unseen
levels correlate at -0.0070, 0.0358 and 0.0042; two rows in one unseen
level correlate at 0.3798 against `tau^2 / (tau^2 + sigma^2)` =
0.3645. Per-row variances, and every standard error `frm_linpred()`
reports at an unseen level, are unchanged: each row sits in exactly one
level's group either way. One follow-on the full suite caught: a
`conditional_effects()` grid writes `NA` for "a group the fit did not
see", and an `NA` label matched no key, so those rows lost the
between-group variance and `test-ce-bands.R` and `test-effects.R`
each failed one expectation (a band that no longer widened). An `NA`
label is now its own unseen level per row, and both files pass.

**BLOCKER D. `coef()` on a mixed ordinal fit.** `ord ~ x + (1 | g)`
gave thresholds repeated unchanged in every group beside a stray
`(Intercept)` column holding the mode alone. brms's `coef.brmsfit()`
moves each threshold by the group's random intercept with the sign the
family's specials give: `thres_minus_eta` families (cumulative, sratio)
report `threshold - r_g`, `eta_minus_thres` families (cratio, acat)
report `r_g - threshold`. frmtmb's four ordinal likelihoods use the
same two forms (`ord_cat_probs()`), so the rule transfers as written.
`dev/shapes-p2-brmsref.R` fits a cumulative and an acat model in brms
and checks which formula reproduces brms's own `coef()` draw by draw:
the cumulative one reproduces `threshold - r` exactly (maximum absolute
difference 0) and misses `r - threshold` by 4.65; the acat one the
other way round (0, against 5.34). brms's column names are
`Intercept[1]`, `Intercept[2]`, `x` and there is no intercept column;
frmtmb's are the same now.

**MAJOR B. The non-finite drop was per replicate.** A replicate with
one non-finite cell was dropped whole, so every other row was
summarized over a SELECTED subset of parameter draws. The mask is per
cell now: that cell is NA, the replicate stands for every other row,
and the warning names the rows and how many replicates each lost. The
parameter draws and one seed per replicate are also taken up front
from the caller's stream, and each replicate simulates from its own
seed. A poisson row at x = 4 predicted alone and beside a row where
`exp(eta)` overflows in 202 of 400 draws: Estimate 12.37, Est.Error
4.343993, interval (5, 21.025) both ways, `identical()` TRUE.

**CORRECTION (punch round 3).** This paragraph said that with those
changes "a row's draws [do] not depend on its neighbours at all", and
NEWS and a code comment said the same. That is too strong. The rows
of one replicate still share one random stream, so a row's
Monte Carlo draws SHIFT with the number of rows simulated before it;
the bitwise identity above holds because the tested row comes FIRST.
What is gone is the SELECTION bias: no row is summarized over a subset
of parameter draws chosen by another row's overflow, so a row's
summary is the same in distribution whatever else is in `newdata`,
and bitwise the same when no row precedes it. In round 1 the same comparison
gave 15.264 [7, 26] alone and 12.947 [6, 21] beside it. The near-flat
case the recheck also found, a draw that is FINITE but absurd
(Estimate 1.4e146), is filed in `dev/test-backlog.md`, not fixed.

**MAJOR C. Multivariate draws ignored the residual correlation.**
Each response was drawn from its own simulator, so on a fit with
rescor estimated at 0.7 the draws correlated at 0.009. Under
`set_rescor(TRUE)` the responses of one replicate are now drawn from
their joint multivariate normal law, with the correlation matrix read
at that replicate's parameters; `rescor = TRUE` is refused for every
family but gaussian, so multivariate t never arises. Joint drawing was
feasible and `summary = FALSE` is not refused. Measured: estimated
rescor 0.5882, mean within-row draw correlation 0.5896 with
`param_uncertainty = FALSE` and 0.5852 with it. brms on the same kind
of fit draws at 0.724 against its own rescor of 0.725
(`dev/shapes-p2-brmsref.R`).

**MAJOR E. `fixef()` on ordinal draws** reported `x` alone, because the
sampler stores the thresholds as `tau_raw`. It maps each draw through
the family's threshold map now and reports `Intercept[1]`,
`Intercept[2]`, `x` in brms's order, the fit's rows. `brms_fixef_rows()`
joined the extension API so frmtmb.sample reads the same rows the fit
does.

**Minors.**
- `ranef()` names its columns as brms does (`Intercept`,
  `sigma_Intercept`), as `coef()` already did; brms says `Intercept` in
  both (`dev/shapes-p2-brmsref.R`).
- `insight::get_parameters()` and `get_varcov()` report an ordinal
  fit's thresholds on the model's scale under `fixef()`'s names, with
  the delta-method covariance. marginaleffects keeps the internal
  vector, because it writes the vector back through `set_coef()`.
- On draws, `allow_new_levels` and `sample_new_levels` are refused by
  name with a reason that is true. They used to vanish into `...`,
  after which the design builder refused the new level with "Use
  allow_new_levels = TRUE". Supporting them on draws means drawing a
  new level's effect per posterior draw, or resampling the seen levels
  as brms's default does; neither is built.
- `predict(summary = FALSE)` names neither draws nor rows. brms returns
  `list(NULL, NULL)` from `predict(summary = FALSE)` and
  `posterior_predict()`, and `list(NULL, NULL, c("y1", "y2"))` on a
  multivariate fit (`dev/shapes-p2-brmsref.R`).
- NEWS says that under `na.omit` the row names of `fitted()`,
  `residuals()` and `predict()`, formerly the way to line results up
  with the surviving rows, are gone, as in brms.

**Records.** Section 3 carries two CORRECTIONS: brms's draws DO carry
group-effect uncertainty at a known level (the posterior of `r_g`), so
leaving `Var(b | y)` out diverges from brms rather than matching it,
with the decision filed in `dev/test-backlog.md` and behavior
unchanged; and the review's round-1 target of 1.02794 double-counted
`tau^2`, the plug-in value being 0.936 and 0.944 on the review's two
fits.

**CORRECTION (punch round 3): round 2's "gated tier" was the ported
suite only.** It ran the fifteen generated `test-brms-suite-*` files
with `FRMTMB_BRMS_FIT_TESTS=true`, and it did NOT set the gate for the
four brms-agreement files that share it, `test-brms-methods.R`,
`test-brms-likelihood.R`, `test-brms-agreement.R` and
`test-brms-priors.R`. In the ungated suite those files skip most of
their blocks (methods pass=16 skip=45, likelihood pass=20 skip=32 in
round 2's log), so their breakage was never seen. Run gated, the
round-2 build gave methods 368 pass and 18 errors and likelihood 329
pass and 7 errors (the recheck's numbers). Section 11 has the fix and
the gated counts.

**Suites and checks, punch round 2.**

- Gated ported tier (`dev/shapes-log/tier-punch2.txt`): 15 of 15 files,
  631 expectations, 0 failures, rerun after the last code change.
  Guards 60 of 60. The ledger rebuilt without a problem and its totals
  are unchanged: 231 pass, 83 defect, 35 divergence, 145 cannot
  transfer. No bin-1 assertion reads anything round 2 changed.
- Full suite (`dev/shapes-log/punch2-final.log`): 266 of 266 files,
  15,766 expectations. Seven files reported a failure; four were this
  round's and are fixed, three are not this lane's:
  - `test-ce-bands.R` and `test-effects.R`: the `NA` level label
    (BLOCKER A above). Fixed; 167 and 54 pass.
  - `test-id-kron.R`, `test-mm.R` and `test-v14.R` read `ranef()`
    columns as `(Intercept)`, `mmc(c1, c2)` and `y1.mu:(Intercept)`;
    `ranef()` gives brms's names now, `Intercept`, `mmcc1c2` and
    `y1_Intercept`. Tests updated; 60, 116 and 81 pass.
  - **Not this lane's**: `frmtmb/test-nl-rtmb-scope.R` (the shadow set
    lacks `atan2`), `frmtmb.sample/test-sample-direct.R` (2) and
    `frmtmb.sample/test-stan-control.R` (1). RTMB 2.0 was installed into
    the shared user library on 2026-09-21 at 10:36, between round 1 and
    round 2, and it newly exports `atan2`. All three files fail the SAME
    way on the untouched base build `rellib-r3` with the base commit's
    own test files (`dev/shapes-run-base.R`: `nl-rtmb-scope` 1 fail and
    1 error, `sample-direct` 2 fails, `stan-control` 1 fail), so the
    change of RTMB is the cause and not this lane. The main tree has
    uncommitted edits to `R/ad-env.R` and `test-nl-rtmb-scope.R` that
    look like the fix for the first; this lane did not duplicate them.
- `R CMD check --as-cran` (`dev/shapes-log/check.log`; round 2's first
  run is kept as `check-punch2-run1.log`): frmtmb `1 ERROR, 1 NOTE`,
  the ERROR being the one RTMB 2.0 test above (`FAIL 1 | PASS 10323`)
  and the NOTE the V8 math-rendering one base has; frmtmb.sample
  `1 ERROR`, the two RTMB 2.0 sampler tests. Nothing else is flagged
  on either package.

## 11. Punch round 3

The recheck confirmed A, B, C, D, E and every minor on the reviewer's
harness, with no standard error moved. What blocked the merge was test
breakage that round 2's gated runs never saw (section 10's CORRECTION),
plus one regression in the new-level refusal on draws.

**BLOCKER. The brms-agreement files, run gated.** Three test-side
faults, fixed in the tests and not in product code:
- `brms_frm_coef()` (`tests/testthat/helper-brms.R`) mapped brms's
  coefficient names onto the OLD `ranef()` columns (`(Intercept)`,
  `y.sigma:(Intercept)`), so round 2's rename stopped it with "no
  unique group-level column". It tries brms's own name first now,
  `<resp>_<dpar or nlpar>_<coef>`, which is what `ranef()` gives; the
  nlpar prefix was a second miss (`a_Intercept` on the nonlinear row 6
  of check C) that the first gated rerun found.
- `test-brms-methods.R` read `fitted()[, "Estimate"]` and then took
  `[, "Estimate"]` of that again, and read the Estimate of a shape whose
  `fitted()` is `n x 4 x K` as if it were `n x 4`. Both sides now go
  through one reader that takes `[, "Estimate", ]` on a 3-d array.
- Three stale names: `colnames(coef()$Subject)` is `Intercept`, `Days`,
  and `re[lvl1, "Intercept"]` in the conditional-effects block.

Run with `FRMTMB_BRMS_FIT_TESTS=true` and `NOT_CRAN=true`, one process
per file (`sh dev/shapes-gated.sh`, logs `dev/shapes-log/gated-*.txt`):

| file | pass | fail | error | skip |
|---|---|---|---|---|
| `test-brms-methods.R` | 954 | 0 | 0 | 0 |
| `test-brms-likelihood.R` | 404 | 0 | 0 | 0 |
| `test-brms-agreement.R` | 187 | 0 | 0 | 0 |
| `test-brms-priors.R` | 103 | 0 | 0 | 0 |

The likelihood file ran to completion every time here, each run
writing its RESULT line; the empty log the reviewer saw did not recur,
so there was no death to diagnose on this machine. The same counts
come out of the gated full suite below.

**MAJOR. The refusal on draws was too wide.** Round 2 refused whenever
`allow_new_levels`, `allow.new.levels` or `sample_new_levels` was
PRESENT, so `allow_new_levels = FALSE`, brms's default, was refused, and
so was `TRUE` with no `newdata` or with known levels; all of those
answered on base. The refusal fires now only for `TRUE` together with a
`newdata` that holds a level the fit did not see, decided by asking the
fit whether the rows build without the flag and build only with it.
`fitted()`, `predict()` and `residuals()` on draws check before they
delegate, so the message names the function called, not
`posterior_epred()` or `posterior_predict()`. The test runs all five
methods: `FALSE` equals the call without it, `TRUE` with no `newdata`
and with known levels equals the call without it, and `TRUE` with an
unseen level refuses naming the function
(`dev/shapes-log/seen-failing-punch3-sample.txt` is the round-2 build
failing it; the file passes 61 of 61 now).

**MINOR. `?royston_parmar`** told users to read
`ranef(fit)[["centre"]][, "time.gamma1:(Intercept)"]`. On the paired
block it recommends the columns are `Intercept` and
`gamma1_Intercept`, and the snippet reads `"gamma1_Intercept"` now.
Run once by hand on a paired df = 1 fit (`dev/shapes-p3-rpsnippet.R`):
the columns are `Intercept` and `gamma1_Intercept`, the slopes run
1.354 to 1.436, and the check comes back empty. frmtmb.spline is
re-roxygenised.

**Suites and checks, punch round 3.**
- Full suite, GATED this time (`dev/shapes-suite.ps1 -Gated`, log
  `dev/shapes-log/punch3-gated.log`): 266 of 266 files, 18,155
  expectations, 17 skips. Every gated file runs: the four agreement
  files at the counts above and the fifteen ported-tier files at 631.
  Three files fail, all RTMB 2.0 fallout that the coordinator has fixed
  on main and this worktree does not carry: `test-nl-rtmb-scope.R`
  (`atan2`), `frmtmb.sample/test-sample-direct.R` (2) and
  `test-stan-control.R` (1).
- `R CMD check --as-cran` (`dev/shapes-log/check.log`): frmtmb
  `1 ERROR, 1 NOTE`, the ERROR the `atan2` test (`FAIL 1 | PASS
  10323`) and the NOTE the V8 one base has; frmtmb.sample `1 ERROR`,
  the two sampler tests (`FAIL 3 | PASS 1745`); frmtmb.spline
  `1 NOTE`, the V8 one. Nothing else is flagged on any of the three.
- Ledger unchanged: 231 pass, 83 defect, 35 divergence, 145 cannot
  transfer.

## 12. Recheck of punch round 3: MERGEABLE

The reviewer reinstalled this worktree into its own library
(`shapesrev-lib`, RTMB 2.0): all eight packages match the source with 0
differing function bodies. The reviewer did not edit this file; this
section records its report.

**Gated brms-agreement files**, one process each with
`FRMTMB_BRMS_FIT_TESTS=true` and `NOT_CRAN=true`:

| file | pass | fail | error | skip |
|---|---|---|---|---|
| `test-brms-methods.R` | 954 | 0 | 0 | 0 |
| `test-brms-likelihood.R` | 404 | 0 | 0 | 0 |
| `test-brms-agreement.R` | 187 | 0 | 0 | 0 |
| `test-brms-priors.R` | 103 | 0 | 0 | 0 |

Base gives 953 on the methods file. The one extra expectation is the
`brmshypothesis` check. The round-2 likelihood runs that died with an
empty log were stopped from outside, not killed by the code.

**The test edits weaken nothing.** Against f8b45ef the expectation
counts are the same except for that one addition, the `skip` counts
are the same in all six files, and no tolerance changed. The edits
rename columns (`(Intercept)` to `Intercept`), read the 3-D `fitted()`
through one reader applied to BOTH sides, and move `fixef()` and
`predict(type =)` reads to `fixef_by_dpar()` and `frm_linpred()`.
`brms_frm_coef()` tries brms's name first and keeps the old mapping as
a fallback.

**Draws and `allow_new_levels`.** `FALSE`, and `TRUE` with no
`newdata` or with known levels only, answer `identical()` to the call
without the argument, in all five methods and for both spellings. `TRUE`
with an unseen level is refused, classed `frmtmb_sample_error`, and
names the function called. A `newdata` with a missing column keeps its
original error in all five methods, with or without the flag and with
an unseen level present, so the detection cannot misfire.

**Nothing from rounds 1 and 2 regressed.** New-level independence,
per-cell masking, joint rescor draws, ordinal `coef()` against brms
and ordinal `fixef()` on draws all reproduce. No standard error moved
on a non-ordinal fit: `vcov()` and the cluster-robust matrix are
bitwise identical on all 11 fit-and-matrix pairs, and the downstream
quantities show 0.000e+00 change. Rd examples run with 0 errors in
frmtmb, frmtmb.sample and frmtmb.spline.

**Left open, pre-existing.** On draws, an unseen level passed without
the flag still gets core's hint "Use allow_new_levels = TRUE", which
then leads to the refusal. Filed in `dev/test-backlog.md`.

## 13. Found at consolidation: two REML defects the lane shipped

The 0.61.0 release run found two defects on fits whose fixed effects are
integrated out (`REML = TRUE`, or `control(profile = TRUE)`), where the
fixed effects are NOT outer parameters. Base 0.60.0 has neither. The
fuzz tier found the first (three `summary_prints` findings); the second
was found by asking whether the first's cause had siblings.

**1. `summary()` stopped on every REML and profile fit**, with "row names
contain missing values", including the plain `y ~ x + (1 | g)`.
`summary_spec_frame()` read standard errors from `vcov(full = TRUE)` by
position. That is the OUTER covariance, which under REML has no
fixed-effect rows, so sigma's standard error was an NA-named NA and
`data.frame()` refused the name. It now reads `vcov_estimated()`, the
covariance over exactly the coefficients `fixef_estimated()` returns,
under every mode.

**2. `predict()`'s interval dropped the fixed-effect uncertainty under
REML**, with nothing said. The parameter draw perturbed the outer vector
only, so `beta` stayed at its estimate. At an extrapolated point
(x = 3, 40 rows):

| mode | analytic sqrt(se^2 + sigma^2) | Est.Error before | after | plug-in |
|---|---|---|---|---|
| ML | 1.1419 | 1.1544 | 1.1604 | 0.9867 |
| REML | 1.1611 | 1.0147 | 1.1802 | 1.0007 |
| profile | 1.1419 | not measured | 1.1604 | 0.9867 |

`fit_draw_space()` adds `beta` to the draw under REML and profile, with
its covariance taken jointly with the outer parameters from the joint
precision, the source `vcov_estimated()` reads. `fit_fd_se()` uses it
too. Under ML it returns `outer_par_map()` and `vcov(full = TRUE)`
unchanged: summary, spec_pars, predict with parameter draws and
fitted() are `identical()` before and after on a gaussian, a
distributional and a cumulative fit.

**Why three reviews missed it.** No test in any tier calls `summary()` or
`predict()` on a REML fit, and every harness in the lane and the review
used ML fits. `tests/testthat/test-reml-shapes.R` now pins both, and was
seen failing on the lane build (`shapes-lib`).

**Also at consolidation, not defects:** the fuzz tier's `vcov_dim`
invariant compared `vcov()` with the estimated-coefficient count, which
2.6f made false by design (`vcov()` is brms's population-level block);
it now compares `vcov_estimated()`. And four ported rows (:345, :350,
:764, :775) HOLD on the merged build, because this lane's shapes and
lane adefects' `fit$data` fix together remove each row's recorded
reason.
