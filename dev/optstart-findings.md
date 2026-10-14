# Lane optstart: starting values that do not put a fit at a false stop

Round of 2026-10-13, base `d2ef6312` (frmtmb 0.70.0), worktree
`frmtmb-wt-optstart`, private library `wt-optstart-lib`, before arm
`rellib-r8`. Nothing is committed. Logs are in `dev/optstart-log/`
(gitignored); every number below is pasted from the log named beside
it. Final build: install 13 (`dev/optstart-log/install-13.txt`).

STATUS: done, with fix round 1 (the review
`dev/reviews/2026-10-13-optstart.md`) applied in build 13.

## Fix round 1 (the review, 2026-10-13)

1. exgaussian starts sigma and beta at rms / 2, the initializer's
   divisor; the comment that claimed an even split of the variance is
   gone. Battery rebig: 5 errors of 10, as 0.70.0 (was 10); the other
   four exgaussian designs keep 0 short, 0 nonzero codes, 0 warned.
2. test-optrap.R is 0.70.0's again (restored from HEAD with `git
   checkout -- tests/testthat/test-optrap.R`, the one write-side git
   command of the lane); it passes, 145 of 145.
3. `restart_result()` relabels only a restart that also did not move
   (every coordinate within `x.tol` of its size or natural unit), reads
   the tolerance of the optimizer that ran (`restart_rel_tol()`), and
   keeps the restart's code, message and iterations in `opt$restart`,
   which diagnose() prints. The review's constructed case (`rel.tol =
   1e-4`, objective 1e4, a gain of 0.9 and a move of 5) is a test and is
   not relabeled. re1000 and R1 still end at code 0.
4. The condition number is one computation, `lp_design_kappa()`, read
   by the start rule and by newton_polish()'s gate.
5. NEWS rewritten: shorter sentences, the family list as a list, no
   "solved once", the "bit for bit" scope qualified, the four
   user-visible changes added, R1's figure the review's independent
   2.1e-6 and 6.0e-7.
6. Filed: the pre-existing short stops at kappa 2e4 to 1e5 (below).
7. For consolidation: lane illse's test-se-conv.R (below).

Rerun on build 13: the 372-file scan (0 fail, 0 error, 0 skip, 0
escaped warnings), the fit diff, the battery, the 24 moving files
under OpenBLAS 0.3.26 with LAPACK (24 of 24 pass), the relations.
No Rd changed in this round, so R CMD check was not rerun.

## Fix round 2 (build 14)

1. NEWS: every sentence of the development section is 25 words or
   fewer (`dev/optstart-rev2-news.R`: 71 sentences, 0 over 25).
2. NEWS on exgaussian: from the true values with the group effects at
   0 the density fails on 10 of 10; from the true group effects, on 3
   of 10 (the review's `optstart-rev-exg.R`).
3. A relabeled restart now carries the first run's `iterations` and
   `evaluations` too, so the whole of `fit$opt` describes that run; the
   restart's code, message, iterations and evaluations are in
   `fit$opt$restart`. NEWS says the step tolerance is `x.tol` (1.5e-8 by
   default) for nlminb and 1.5e-8 for `optimizer = "optim"`.

Rerun on build 14 (`dev/optstart-log/r2/`), test-optstart.R and every
test file that reads `iterations` or `evaluations`: test-optstart 69
pass (two assertions of the restart test merged into one), test-optrap
145, test-se-check 127, test-se-tier 91, test-diagnostics-ux 130,
test-edgecases 64, test-fd-chain 11, test-nl-rtmb-scope 105,
test-predict-re-uncertainty 78, test-brms-port 21, frmtmb.eam
test-defects 60 and test-gddm-gradients 22; 0 fail, 0 error, 0 skip, 0
warning in each. No git command was used in this round.

## What changed

- `R/start-ls.R` (new). `start_ls_plan()` and `start_ls_one()`: the
  least-squares start, per response; `start_ls_rules`, one entry per
  family (the link of `mu`, the scale parameters it sets, their relation
  to the residuals); `ls_rms()`, `trigamma_inverse()`, `lp_key_of()`;
  `response_observed()` (an `mi()` row that is imputed is not data).
- `R/fit.R`, `make_start()`: where the plan applies, `mu`'s
  coefficients are the least-squares ones and each scale's intercept is
  the link of its relation; the residuals feed the initializers that
  read them (skew_normal's sigma and alpha); the response unit reuses the
  same solve. Elsewhere the old code runs unchanged.
- `R/fit.R`, `optimize_obj()`, `restart_result()` and
  `restart_rel_tol()` (new): a restart from a converged run that ends
  at a nonzero code, gains no more than the optimizer's relative
  function tolerance and moves no coordinate by more than its relative
  step tolerance keeps the run's code; its own code is kept in
  `opt$restart`. `R/confint.R`, `diagnose()`: shows `opt$restart`.
- `R/autoscale.R`: `predictor_ls()` (new) is the least-squares solve
  `predictor_resid_sd()` made, now with the coefficients and a weighted
  form; `lp_design_kappa()` (new) is one predictor's condition number,
  which `predictor_design_kappa()` (newton_polish()'s gate) and the
start rule now share; `predictor_unit_input()` and
  `predictor_scale_response()` (new)
  are the first half of `predictor_unit_of()`; `predictor_units()` and
  `predictor_unit_of()` take the fit the start made. The unit of
  zero_inflated_asym_laplace() reads its nonzero rows only, a weighted
  response's unit reads the weighted fit, and an imputed `mi()` row is
  left out.
- `R/autoscale.R`, `predictor_col_spread()` (item R1): a `mi()` or
  `me()` placeholder column spreads as the values it stands for
  (`spread_mi_values()`, `spread_me_values()`, `placeholder_values()`,
  new).
- `R/se-check.R`, `fd_hessian()` (item R1): a matrix with no Cholesky
  factor is retaken along the coordinates whose step raised the
  objective by more than 1e-2, at the step that raises it by 1e-4; the
  retake is kept when it factors and its smallest eigenvalue clears its
  own noise.
- `R/par-template.R`, `R/fit.R`: the section "The default starting
  values" of `?par_template` and the `start` argument of `?frm`.
- `NEWS.md`: three breaking-change bullets and two bug fixes under
  `# frmtmb (development version)`.
- Tests: `tests/testthat/test-optstart.R` (new, 9 tests). Changed:
  `test-se-check.R` and
  `test-se-tier.R` (the gr-by fixture), `test-grad-verdict.R` (three
  fits), `test-structure.R`, `test-skew-normal-start.R` (two tests),
  `test-car-spde.R` (one assertion), and frmtmb.sample's
  `test-gr-by-draws.R` (an allowed warning). Each gives the fit 0.70.0's
  start where the test's premise is the stop from that start, or says
  why the assertion moved; "Tests changed" below has the measurement
  for each.

## The rule

`start_ls_plan()` applies to a response when all of these hold.
Anywhere else `mu` and the scale start as on 0.70.0, bit for bit. The
group sds and the units can differ there too, because the unit of an
`mi()` response and of zero_inflated_asym_laplace() changed: on `ym |
mi() ~ s(x) + (1 | g)` and an intercept-only zero_inflated_asym_laplace
the unit is 8 against 16 and 32 against 16 (the review's
`optstart-rev-absent.R`). See "Where the rule does not apply".

- The family is one of: gaussian, student, skew_normal, exgaussian,
  asym_laplace, zero_inflated_asym_laplace, huber, lognormal,
  shifted_lognormal, hurdle_lognormal (identity link of `mu`); Gamma,
  weibull, negbinomial (log link; Gamma's default link is the inverse
  and is left out); beta (logit link).
- `mu` is linear (not nonlinear, constant or tied), has a covariate or
  an offset (not intercept only), has no `s()`, `t2()` or `gp()` term,
  and its design's condition number (columns at unit length) is at most
  1e5.
- Every scale it sets is estimated; the response's addition terms are
  among `weights()`, `subset()`, `index()`, `mi()`.
- `start` and priors do not turn it off: the entries `start` names
  are written over the result, and a prior on a coefficient moves no
  start (test-par-template.R pins both). A bound from a prior is met
  by the optimizer's projection of the start (nlminb: start 10, upper
  5, reaches 3).

The least-squares fit is on `mu`'s fixed design, the offset taken out,
weighted by the likelihood weights, on evenly spaced rows when the
design has more than `frmtmb.unit_ls_cells` (5e6) entries: the solve
the response unit makes, so make_start() solves it once for the
identity-link families. The log and logit families solve on their own
response, and `fit_assembled()` solves the unit's once more (filed). The
response on the link scale: y; log y (lognormal, Gamma, weibull, with
the shifted lognormal's shift and the hurdle's positive rows); the
logit of y (beta); glm()'s first IRLS step from `y + 0.1` with its
working weights (negbinomial).

### Each family's relation, measured

`dev/optstart-relations.R` (20 seeds, n 2000, data from the family's own
simulator with a covariate that explains most of the response): the
cold start of the scale over its true value (median and range), and of
the slope, on the final build (`relations-lane11.txt`) and on rellib-r8
(`relations-r8.txt`).

| family | scale | relation | lane: start / true | 0.70.0 | slope, lane |
|---|---|---|---|---|---|
| gaussian | sigma | rms of residuals | 0.996 [0.970, 1.034] | 13.363 | 1.000 |
| student | sigma | MAD / qt(0.75, 4) | 1.010 [0.946, 1.064] | 13.416 | 1.000 |
| skew_normal | sigma | sd(resid) (family) | 0.998 [0.965, 1.044] | 0.998 | 1.000 |
| exgaussian | sigma, beta | each rms / 2 | 0.707 [0.678, 0.731] (true sigma = beta: the relation is 1 / sqrt(2)) | 6.711 | 1.000 |
| asym_laplace | sigma | mean abs. dev. / 2 | 1.010 [0.972, 1.044] | 10.106 | 1.002 |
| zi_asym_laplace | sigma | the same, nonzero rows | 1.000 [0.960, 1.048] | 10.121 | 1.000 |
| huber | sigma | MAD | 1.071 [1.024, 1.117] | 13.283 | 1.000 |
| lognormal | sigma | rms on log y | 0.996 [0.970, 1.034] | 2.243 | 1.000 |
| shifted_lognormal | sigma | rms on log(y - min/2) | 0.957 [0.913, 0.997] | 2.133 | 0.956 |
| hurdle_lognormal | sigma | rms, positive rows | 1.004 [0.969, 1.035] | 2.252 | 1.002 |
| Gamma (log) | shape | trigamma(shape) = var(log y) | 1.003 [0.947, 1.045] | 0.226 | 0.992 |
| weibull | shape | pi / (sqrt(6) rms of log y) | 0.998 [0.971, 1.043] | 0.600 | 1.004 |
| beta | phi | mean mu(1-mu) / mean (y-mu)^2 - 1 | 0.935 [0.872, 0.992] | 0.186 | 1.094 |
| negbinomial | shape | sum mu^2 / sum((y-mu)^2 - mu) | 0.875 [0.696, 1.020] | 0.227 | 0.888 |

Gamma and weibull move the intercept by `log(k) - digamma(k)` and
`lgamma(1 + 1/k) + euler/k` (the mean of `log y` is not `log mu`);
negbinomial moves it so the fitted mean is the mean. skew_normal's
sigma and alpha come from its own initializers, which already read the
residuals; only the residuals' source changed. exgaussian's rule went
through three versions: a third-moment split, then rms / sqrt(2) to
each, then (fix round 1, from the review's measurement) rms / 2 to each,
the divisor of the family's own initializer (`relations-lane13.txt`;
"Battery").

## The two backlog cases

**Wide high-R^2 design** (`dev/optstart-cand.R wide`,
`dev/optstart-cases.R`; n 101000, 50 columns, seeds 1 to 3). rellib-r8,
seed 1: code 0, logLik - lmer -9108 (the backlog's -9107.66), 3
iterations (`cand-wide-r8.txt`). Build 12 (`cases-wide-lane12.txt`; fix
round 1 does not reach this design):
seed 1 code 1 at -3.912e-08, seed 2 code 0 at -1.397e-09, seed 3 code 0
at -2.733e-08. The silent 9107.66 stop is gone.

Seed 1 ends at "false convergence (8)" at lme4's maximum. That code is
not the start's: with `verbose = 2` the first run takes 254 iterations
from 146180.05 to 144801.20 and stops at code 1 itself
(`wide-trace2-lane2.txt`), and across three starts given through
`start =` (`dev/optstart-wide-arms.R`, `wide-arms-s{1,2,3}.txt`; first
run alone, `restarts = 0`) the code-1 stops land on different seeds for
each start:

| start | seed 1 | seed 2 | seed 3 | iterations |
|---|---|---|---|---|
| rms (the lane's) | 0 | 1 | 0 | 255, 265, 263 |
| within-group sigma | 1 | 0 | 0 | 177, 221, 183 |
| within sigma, anova group sd | 1 | 0 | 0 | 83, 69, 48 |

All 18 fits of those runs (first run alone and default) are within
1.1e-7 of lme4. At n 1e5 nlminb's choice between "relative convergence"
and "false convergence" at the optimum follows rounding whatever the
start; the anova start halves the time (83 to 100 s against 147 to 235)
but not the codes, so it was not built (filed). The warning on seed 1
is the convergence-warning policy's, which this lane does not change.

**Group sd 1000 times the residual sd** (`dev/optstart-cand.R re1000`,
20 seeds). rellib-r8: 18 code 0, 2 code 1 (seeds 6 and 14,
`cand-re1000-r8.txt`). The backlog's guess, that the sigma start would
remove it, is wrong: the least-squares start alone (build 1) gave 19
code 0 and 1 code 1 (seed 8, `cand-re1000-lane1.txt`). The cause,
traced on seeds 6 and 14 (`dev/optstart-re1000-trace.R`,
`re1000-trace-r8.txt`): the first run ends at code 0 at lme4's maximum,
its gradient in natural units (0.0021, 0.00104) is above `grad_tol`
(1e-3), the restart from there stops in one iteration at "false
convergence (8)" at the same objective, and `optimize_obj()` took its
code because its objective was no higher. Final build: 20 of 20 at
code 0 (`cases-re1000-lane12.txt`, and the same on build 13,
`cases-re1000-lane13.txt`), logLik - lmer from -8.282e-11 to
7.958e-13.

The restart rule went through three versions. The first kept the
converged run's point; under OpenBLAS 0.3.26 with LAPACK that left
test-bcm-simple.R's SIMPLE_1 with a gradient of 0.0025 against the
test's bound of 1e-3 (`optstart-ob26-lane-suite`), where 0.70.0 (which
took the restart's point) passed, and with 0.70.0's rule patched back
in (`dev/optstart-toggle.R ... norestart`) the lane passed it too. The
second kept the restart's point and the converged run's code whenever
the gain was within `rel.tol`; under OpenBLAS the file then passed 13 of
13 (`optstart-ob26-lane12.sum`). The review showed that a caller's
loose `rel.tol` would relabel a restart that moved 5 units and gained
0.9 at an objective of 1e4, that optim's tolerance was read as
nlminb's, and that the relabeled `opt` hid the restart's own code. The
final rule (fix round 1) relabels only a restart that also moved no
coordinate by more than `x.tol` (1.5e-8) of its size or natural unit,
reads the tolerance of the optimizer that ran (`restart_rel_tol()`:
nlminb's `rel.tol`, L-BFGS-B's `factr` times epsilon, optim's
`reltol`), and keeps the restart's code, message and iterations in
`opt$restart`, which diagnose() prints. Every relabel the review saw on
real fits moved by exactly 0, and the re1000 and R1 cases still end at
code 0 (test-optstart.R, 70 pass). The review's constructed case is a
test. test-v14.R's asym_laplace fit, relabeled by the second version,
now keeps its own code 1: its restart moved (both builds warn there).

## Item R1 (from the coordinator: lane remeasure's R1)

Repro copied as `dev/optstart-r1-repro.R`.

Mechanism, measured:

1. The `mi(chl)` and `mi(chl):age` columns are zero placeholders in
   `X`, so `predictor_col_spread()` gave them spread 1 and their
   coefficients the response's unit itself.
2. The imputed rows of an `mi()` response are stored as 0 in
   `frame$y`, and the unit read them as data: bmi (9 of 25 rows
   imputed) and chl (10 of 25). On rellib-r8 the units are 8 8 8 8 128
   128 1 1 (`dev/optstart-r1-peek.R`); on the lane 2 2 0.0625
   0.0078125 64 64 1 1 (`dev/optstart-r1-peek3.R`).
3. The SE route takes `fd_hessian()` in `par / u` with a step of 1e-3:
   with the units of 1 and 2 the steps of `bmi_michl` and
   `bmi_michl:age` raise the objective by more than 1e-2 (columns 3 and
   4, `dev/optstart-retake.R`, `retake-lane9.txt`), and the matrix is
   indefinite (smallest eigenvalue -0.0115 in `par / u`), where a
   Richardson Jacobian is positive definite (0.0016 to 5.45e5).
4. The false convergence on fit_imp2 is the restart (the re1000
   mechanism): on rellib-r8 with `restarts = 0` the fit ends at code 0;
   with the default restart at "false convergence (8)" at the same
   objective (`dev/optstart-r1-attrib.R`, `r1-attrib-r8.txt`). The
   first run's gradient in the units of 1 and 2 is what starts the
   restart. On the lane the fit ends at code 0 from 0.70.0's start and
   from the lane's, with and without the restart (`r1-attrib-lane3.txt`).
   Which change alone would suffice was not measured.

Standard errors against a Richardson Jacobian of the AD gradient
(`dev/optstart-r1-hess.R`; `r1-hess-r8.txt`, `r1-hess-lane12.txt`;
build 13 the same, `r1-hess-lane13.txt`),
largest relative error over the six coefficients:

| fit | rellib-r8 | lane (final) | optimHess in raw units (0.67.0's route) |
|---|---|---|---|
| fit_imp2 | all NaN, code 1 | 3.3e-7, code 0, no warning | 1.5e-3 |
| fit_imp3 | up to +35%, sigmas withheld | 3.3e-7, code 0, no warning | 1.5e-3 |

The figure NEWS quotes is the reviewer's independent one
(`dev/optstart-rev-r1ref.R`, a closed-form marginal likelihood with the
missing values integrated out, `r1ref-lane.txt`): within 2.1e-6 on
fit_imp2 and 6.0e-7 on fit_imp3. The Richardson figure above is from the
same AD gradient frmtmb uses, so the reviewer's is the stronger one.

So 0.67.0's standard errors, which lane remeasure took as the
reference, are themselves off by up to 1.5e-3 (bmi_age).

The step rule. The coordinator asked for a rule that does not depend on
the units. `fd_hessian()` keeps its step of 1e-3 in the coordinates it
is given, and judges each column by the rise of the objective over the
step, `tau = 0.5 H_ii 1e-6`, the same in every coordinate system.
Measured on fit_imp2 and fit_imp3 under 0.70.0's units, the lane's and
none (`dev/optstart-r1-steps.R`, `r1-steps-lane4.txt`): a second pass
at `tau = 1e-2, 1e-3, 1e-4` gives SE errors 2.0e-4, 2.0e-5 and 2.0e-6
under all three unit choices alike. Retaking every column above the
band (build 7) moved the standard-error check's verdict on
test-se-tier.R's curved ridge, a direction flat by construction; so did
retaking only where the first matrix does not factor (build 8), because
the retake factored there by its noise alone (smallest eigenvalue
1.5e-6 against a noise of 1.1e-3; R1's retake 1.03 against 0.02;
`retake-lane8.txt`, `dev/optstart-sediag.R`). The final rule retakes
only where the first matrix does not factor, and keeps the retake only
when it factors and its smallest eigenvalue clears its noise. Under
0.70.0's units it gives SE errors of 4.6e-7 on both fits
(`r1-steps-lane9.txt`). What it does not do: a matrix that factors
keeps its step, so a large step that leaves the matrix positive
definite is not caught (raw units on fit_imp2: 1.5e-3, as 0.67.0).

### R1 across the suites

`dev/optstart-profile.R` logs, for every fit whose `par_units` are not
all 1 (at most 40 outer parameters, gradient under 0.2 s), the route's
Hessian (`fit_outer_hessian()`) against a Richardson Jacobian of the
gradient. `dev/optstart-hess-sum.R optstart-base2 optstart-lane12`
(`hess-sum-lane12.txt`; build 13 gives the same counts,
`hess-sum-lane13.txt`):

| | rellib-r8 | lane |
|---|---|---|
| fits on the path (test files) | 2424 (131) | 2423 (131) |
| Richardson positive definite, code 0, no warning | 2307 | 2308 |
| of those, route not positive definite | 1 | 1 |
| of those, SE error above 1e-3 / 1e-2 / 0.1 | 22 / 11 / 8 | 20 / 7 / 5 |
| a first-pass step rising above 1e-2 | 38 | 37 |

R1's mechanism, a step that rises above 1e-2 and leaves the route
indefinite where Richardson's matrix is positive definite, occurs on
one suite fit: test-optstart.R's `mi()` fit on rellib-r8 (rise 14.7,
code 1, warned); on the lane its steps rise less than 1e-2. Every other
fit with a rise above 1e-2 has a positive definite route with SE error
at most 1.7e-4, or a Richardson matrix that is not positive definite
either (flat directions the SE check judges), and that list is the
same on both builds. The clean fits with SE errors above 1e-2 have
rises of 4e-4 or less, so their errors are not oversized steps; they
are, on both builds, binomial `diag(1 + x | g)` fuzz fits (0.346), a
weibull `trunc()` fuzz fit (0.011, 0.023), a negbinomial `(fac || g)`
fuzz fit (0.028 on rellib-r8, 0.126 on the lane) and frmtmb.spline's
test-difference.R (7.08), and on rellib-r8 only the ported fixture
`(1 + mmc(Age, volume) | mm(patient, visit))` (0.342), whose group sds
have Richardson SEs of 5e3 to 1.3e5 on both builds and NaN from
`vcov()` on both (`dev/optstart-mm.R`, `mm-r8.txt`, `mm-lane11.txt`):
flat directions where neither matrix determines an SE. No other `mi()`
or multivariate fit of the suites shows the R1 symptom.

New tests that fail on rellib-r8 (`t-optstart-r8-final.txt`): "a mi()
coefficient is stepped in its variable's units" (code 1, units 1.0 not
32, SEs NaN), "an imputed value is not data for the start or the unit"
(intercept 34.4 against 49.8, unit 16).

## Where the rule does not apply (the absent case)

`dev/optstart-guard.R` prints the cold start in hex on rellib-r8 and on
the lane; `dev/optstart-guard-cmp.R` compares (`guard-cmp-lane11.txt`;
build 12 changed no start). Every value of every absent case is
identical in hex: binomial (2 of 2), cens (3), Gamma inverse link (3),
Gamma without intercept (2), gaussian log link (3), hurdle_gamma (4),
intercept only with a group term (13), a degree-8 raw polynomial (10),
mixture (7), nonlinear (3), cumulative (5), poisson (2), se() (3),
sigma fixed (3), a smooth (12), student with se() (4), trunc (3),
tweedie (4). With `start = list(theta = 0)` and with a prior on `b`,
11 of 14 values match: theta (the named entry) and the group effects;
beta and sigma are the rule's, as without a start or prior.

The test "where the rule does not apply the start is the initializer's"
pins censoring, poisson, intercept only, nonlinear, a smooth, the
ill-conditioned design and the Gamma inverse link with
`expect_identical()` against 0.70.0's formulas; it passes on rellib-r8
too (`t-optstart-r8-final.txt`), as an absent-case test should.

## Battery

`dev/optstart-battery.R` (11 families, 5 designs, 10 seeds, n 400 in 20
groups; references lme4 for gaussian and lognormal, glmmTMB for student,
Gamma, beta and negbinomial, survival::survreg for weibull, otherwise
the best of both builds and a fit started at the true values),
summarized by `dev/optstart-battery-sum.R r8 lane13`
(`battery-sum-lane13.txt`, build 13, after fix round 1). Short means
more than 1e-4 below the best;
an error counts as short.

| family | design | rellib-r8: short (worst), code != 0, warned | lane |
|---|---|---|---|
| student | hir2 | 10 (-31.4), 0, 10 | 0, 0, 0 |
| student | offterm | 10 (-23.6), 0, 10 | 0, 0, 0 |
| skew_normal | hir2 | 8 (-62.6), 0, 8 | 0, 0, 0 |
| exgaussian | hir2 | 10 (-20.9), 0, 9 | 0, 0, 0 |
| exgaussian | offterm | 10 (8 errors, -1370), 10, 2 | 0, 0, 0 |
| exgaussian | rebig | 5 (5 errors, seeds 2 5 8 9 10) | 5 (5 errors, seeds 2 3 8 9 10) |
| gaussian | offset | 0, 0, 0 | 0, 1, 1 (at the optimum) |
| asym_laplace | all five | 31 of 50, 49 code != 0 | 34 of 50, 49 |
| huber | all five | 24 of 50, 39 code != 0, 4 degenerate | 30 of 50, 40 |

Every other row (gaussian, lognormal, Gamma, weibull, beta, negbinomial
in all designs; skew_normal and student in the rest) is 0 short, 0
nonzero codes, 0 warned on both builds. Totals: short 108 to 69,
nonzero code 103 to 95, warned 129 to 91; errors 13 to 5; degenerate
stops (more than 50 above the true-start fit) 4 to 0. Build 12 (rms /
sqrt(2) for exgaussian) had 74, 100, 91 and 10: the only rows that
differ between builds 12 and 13 are exgaussian rebig's.

asym_laplace and huber have kinked objectives, and nlminb stops at a
nonzero code on most fits from either start. How far short
(`dev/optstart-battery-kink.R`, `battery-kink-lane13.txt`):

| | at best | to 0.1 | 0.1 to 1 | 1 to 10 | beyond 10 | sum of gaps |
|---|---|---|---|---|---|---|
| asym_laplace, rellib-r8 | 19 | 10 | 9 | 4 | 8 | -5530.0 |
| asym_laplace, lane | 16 | 8 | 16 | 9 | 1 | -40.3 |
| huber, rellib-r8 (4 degenerate) | 22 | 19 | 0 | 0 | 5 | -466.9 |
| huber, lane | 20 | 20 | 10 | 0 | 0 | -3.3 |

The lane removes the large stops and the degenerate ones (huber at
+1.4e5 on 400 rows, rellib-r8 seeds 2, 3, 4, 7 of hir2), and stops 0.1
to 1 short more often. exgaussian rebig is no worse than 0.70.0 after
fix round 1 (5 errors of 10 on both, 4 of the 5 seeds shared). The
density (RTMBdist::dexgauss) fails there from the true values with the
group effects at 0 too (`ll_truth` NA on all 10), and on 3 of 10 from
the true group effects (the review's `optstart-rev-exg.R`). Build 12's
divisor, rms / sqrt(2), raised the errors to 10 of 10: the review's
factorial shows the sigma divisor alone did it (`exg-lane.txt`), and
rms / 2 restores 5 of 10 with every other exgaussian gain kept
(`exgfix-lane.txt`, and this battery).

## Cost of the start

`dev/optstart-cost.R` (minimum of 3, elapsed seconds; `cost.txt`,
builds 1 to 2, whose solve is the final one): `make_start()` takes 0.24
s on the wide design (rellib-r8 0.20 s, which made the unit's solve
too), 0.34 s on 1e6 by 5 (0.29 s) and 1.12 s on 1e5 rows with a
500-level factor (1.11 s); the fits take 140 to 250 s (wide), 140 s and
93 s (`cost-fit.txt`). On the factor design rellib-r8 stopped at
objective 153976.08, code 0, and the lane at 141661.42 (not checked
against lme4). The one solve serves the start and the unit;
`fit_assembled()` solves the unit once more, as 0.70.0 did (filed).

## The suites, and the fit diff

Every test file of core and the seven extensions, one file per R
process, every gate on (`dev/optstart-scan.sh`, `dev/ciharden-scan.sh`
with the fit logger `dev/optstart-profile.R`), on rellib-r8 (run
`optstart-base2`) and on the final build (run `optstart-lane13`, build
13, after fix round 1). Generated by `dev/optstart-counts.R`
(`counts-lane13.txt`,
`counts-base2.txt`):

<!-- BEGIN GENERATED: dev/optstart-counts.R optstart-lane13 wt-optstart-lib -->
372 files; 372 with a RESULT line; 372 whose lib: line names wt-optstart-lib

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 230 | 19016 | 0 | 0 | 0 | 0 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1765 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 51 | 2806 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **372** | **26159** | **0** | **0** | **0** | **0** |

files with a failure, an error, an escaped warning or no result: 0 
<!-- END GENERATED -->

<!-- BEGIN GENERATED: dev/optstart-counts.R optstart-base2 rellib-r8 -->
372 files; 372 with a RESULT line; 372 whose lib: line names rellib-r8

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 230 | 18967 | 31 | 1 | 0 | 7 |
| frmtmb.coupling | 11 | 548 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 29 | 1765 | 0 | 0 | 0 | 0 |
| frmtmb.latent | 10 | 366 | 0 | 0 | 0 | 0 |
| frmtmb.learn | 15 | 515 | 0 | 0 | 0 | 0 |
| frmtmb.ode | 11 | 551 | 0 | 0 | 0 | 0 |
| frmtmb.sample | 51 | 2806 | 0 | 0 | 0 | 0 |
| frmtmb.spline | 15 | 592 | 0 | 0 | 0 | 0 |
| **total** | **372** | **26110** | **31** | **1** | **0** | **7** |

files with a failure, an error, an escaped warning or no result: 1
- frmtmb/test-optstart.txt: pass 21 fail 31 err 1 skip 0 warn 7
<!-- END GENERATED -->

Every file passes on the final build, with no skip and no escaped
warning. rellib-r8's failures are test-optstart.R, the new file (an
earlier version of it in that run; the final version on rellib-r8 is
pass 25, fail 33, error 1 before fix round 1's extra tests,
`t-optstart-r8-final.txt`). The ported brms-suite files pass on both
builds with the same counts, so no ledger row flips.

Fits matched by test file, a digest of the response and the call, and
occurrence (`dev/optstart-fitdiff.R optstart-base2 optstart-lane13`,
`fitdiff13.txt`): 10007 fits logged on each build in 315 files; 7761 at
the same optimum (within 1e-8 relative), 13 higher on the lane, 14
lower, 125 errors on both, 186 unfitted or non-finite. 1908 and 1908 fits do
not match (their call or data differ between the runs, for example a
test whose call now passes `start =`); they are covered by their files'
results above.

Fits whose optimum moved:

| file | fit | lane - base (objective) | codes base, lane | why |
|---|---|---|---|---|
| test-optstart.R | the wide design, seeds 1, 14, 18 | -61.59, -33.28, -32.65 | 0, 0 | the case this lane fixes |
| test-fuzz.R (4) | `y \| weights(w) ~ 1 + x + ar1(...)`, beta | +8.6e14 to +1.1e15 | 1, 1 | rellib-r8 reached objectives near -1e15, a degenerate end of the likelihood; the lane stops at a finite point; both warned |
| test-fuzz.R (2) | the same with `mo(xo) * z`, Gamma | +3937, +1557 | 1, 1 | rellib-r8 at -4074 and -1740 (log likelihood +4074 on fuzz data), the lane at -137 and -183; both warned |
| test-fuzz.R (2) | the same, weibull | +31.53, -18.64 | 1, 1 | two non-converged stops each way |
| test-fuzz.R (4) | `mo(xo) * z + (1 \| ga/gb)`, student | -0.0020 | 0, 0 | a higher maximum on the lane, both at code 0 |
| test-grad-verdict.R | "a start far away" (start beta only, 1 iteration) | +1.45e5 | 1, 1 | sigma, which `start` leaves out, now starts at the residual spread; the test is about the cap |
| test-se-check.R | `ls ~ mo(inc) + age`, `iter.max = 3` | +4.95 | 1, 1 | a stop at the cap from another start; the test is about the cap |
| test-gr-by-draws.R (sample) | the gr(g, by = f) fixture | -1.44e-4 | 0, 1 | the lane reaches the boundary (the correlation at its end), "singular convergence (7)"; 0.70.0 stopped short (`dev/optstart-grby.R`); test-se-check.R and test-se-tier.R now give that fit 0.70.0's start, so they match |
| test-v14.R | asym_laplace | -1.84e-4 | 1, 1 | nearer the optimum; both warned (build 12 relabeled its restart to code 0; the restart moved, so build 13 keeps its code) |
| test-v14.R | zero_inflated_asym_laplace | +5.7e-4 | 1, 1 | further from the maximum: 1.25e-3 short of the exact one (rq and glm) against 6.85e-4 on 0.70.0 (the review's `verify-*.txt`); both stop at "false convergence (8)" and warn twice |
| test-importance.R (4) | `y ~ x + (1 \| g)` | -2.0e-5, -6.2e-6 | 0, 0 | nearer the optimum |
| test-brms-pins.R (sample) | `(1 + x + z + w \| g)` | +8.9e-6 | 1, 1 | further from the maximum: 1.72e-5 short of lme4 (singular) against 8.2e-6 on 0.70.0; both warn |
| test-gp-by.R | `ysig ~ x` | +5.4e-6 | 0, 0 | the optimizer's precision (5.4e-6 on 53.9) |
| test-review-v29.R | mvbf with `(1 \| q \| id)` | +8.9e-7 | 1, 1 | both warned |

Same optimum, convergence code changed: 13 to code 0 and 6 to a
nonzero code (`fitdiff13.txt`). To code 0: six frmtmb.spline
test-rp-floored.R fits,
frmtmb.ode test-ode.R (blowup) and test-scale.R, test-bcm-simple.R,
test-xbeta-zibb-hurdle-cum.R, test-optstart.R's three restart cases
(re1000 and the `mi()` fit): the restart's "false convergence (8)" at
the same objective (`restart_result()`); their warnings go. To a
nonzero code, each at the same optimum: three test-fuzz.R
`mo(xo) + ar1(...)` Gamma fits (objective 4e-7 lower), test-optrap.R's
negbinomial at the Poisson limit (`c ~ x`, false convergence; the shape's
maximum is at infinity), two test-sparsex.R fits (`y ~ f + x + z + (1 |
g)`, singular convergence): the first run from the least-squares start
stops at that code where 0.70.0's did not. They warn; the tests pass.

Warnings that appear or disappear on matched fits: 34 fits
(`fitdiff13.txt`, "warnings that appear or disappear"; 35 on build
12, whose extra one was test-se-tier.R's gr-by fit). 16 are the
code changes above at the same optimum (the three fuzz Gamma fits'
code changes carry no changed warning); 2 change a warning's details
between two warned stops at the same optimum (test-se-tier.R's
large-gradient warning names a theta where it named sigma, test-v14.R's
zero_inflated_asym_laplace names `(Intercept)` where it named `x`); 16
are fits whose optimum moved, in the table above, among them the
SE-check warnings test-optstart.R's wide design no longer gives and the
gr-by fixture's singular convergence.

Under OpenBLAS 0.3.26 with LAPACK, 4 threads (`dev/ciharden-openblas.sh
0.3.26 lapack`; runs `optstart-ob26-base12`, `optstart-ob26-base13`
(the three files fix round 1 changed) and `optstart-ob26-lane13` over
the 24 files above that move or that this lane edited,
`dev/optstart-log/ob-list12.txt`): every file passes on the final
build (24 of 24, 0 fail, 0 error, 0 warning), test-bcm-simple.R with
the Stan tier among them; on rellib-r8 test-optstart.R fails (28
pass, 30 fail, 1 error) and test-skew-normal-start.R fails the one
assertion this lane changed (69 pass, 1 fail), as on the reference BLAS.
The first restart rule failed test-bcm-simple.R there (see "The two
backlog cases").

## Tests changed, and why

- test-optrap.R: reverted in fix round 1. Its two designs have
  condition numbers 5.77e6 and 7.52e7, so the rule never applies there
  and 0.70.0's file passes (145 pass; the edits dated from before the
  condition-number bound).
- test-se-check.R, "a small loading on a downward direction keeps its
  SE", and test-se-tier.R, "ranef(condVar = TRUE) reads the repaired
  joint covariance": the gr(g, by = f) fixture. From the least-squares
  start the fit reaches 1.44e-4 lower, at the boundary, with "singular
  convergence (7)", where the SE check does not run
  (`dev/optstart-grby.R`); both tests are about the SE verdicts at
  0.70.0's stop and give its start. That stop is 1.44e-4 short of the
  maximum (filed).
- frmtmb.sample test-gr-by-draws.R: the same fixture; the fit is not
  what the test tests, so the singular-convergence warning is allowed
  beside the SE warning.
- test-grad-verdict.R, three fits ("a gaussian GLM at n = 50000 trips
  the trip-wire", "a fit stopped at its iteration cap still warns",
  "grad_tol tightens"): from least squares a gaussian GLM starts at its
  closed form, so there is no gradient to trip and no headroom; the
  tests give 0.70.0's start.
- test-structure.R, "the family applies the row weights": struct_wgauss()
  is not in the rule and keeps its start; the gaussian reference gets
  the same start, so the two stop at one point (the deviance ratio's
  spread was 4.7e-6 against the test's 1e-6).
- test-skew-normal-start.R: the sparse() case now reads the
  least-squares intercept (the first level's mean); the stationary-point
  case gives mu its old start too, since from least squares the fit
  leaves the stationary point without the refit.
- test-car-spde.R: the sd's relative gap at con_sd 1e-5 was asked to be
  a fiftieth of the default's; it stops falling near 5e-7 from con_sd
  1e-4 on (5.4e-7 and 6.7e-7 on the lane, 6.7e-7 and 3.9e-7 on
  rellib-r8, `dev/optstart-car.R`), and `rel.tol = 1e-14` does not move
  it (`dev/optstart-car-tight.R`); the assertion asks a tenth.
- test-par-template.R is unchanged: its tests of a named `start` and of
  a prior that moves no start are why the rule does not switch off for
  either.

## Filed (found, not fixed)

- exgaussian's density (RTMBdist::dexgauss) ends at "NA/NaN gradient"
  on a group sd 30 times sigma on 5 of 10 fits on both builds, and from
  the true values too. Pre-existing; no longer a regression.
- Pre-existing, from the review (`optstart-rev-collin.R`,
  `collin-*.txt`): at a condition number of 2e4 to 1e5 with `(1 | g)`,
  where newton_polish() does not run, 7 of 40 fits stop short of lme4
  at code 0 with no warning, up to 0.143 on the lane (seeds 7, 8, 9)
  and 0.038 on rellib-r8 (seeds 2, 5, 10). The start moves which seeds
  stop short; it neither causes the defect nor fixes it.
- huber() reached a log likelihood of +1.39e5 on 400 rows on rellib-r8
  (battery hir2 seeds 2, 3, 4, 7): a degenerate stop of its density.
- 0.70.0's own initializers read an `mi()` response's imputed zeros:
  `bf(y | mi() ~ 1)`, outside the rule, starts mu at the mean with the
  zeros in (34.4 against 49.8 in the new test's data with `~ x`).
- The wide design ends at "false convergence (8)" at lme4's maximum on
  1 of 3 seeds from any start; an anova start halves its time.
- `fit_assembled()` solves the unit's least squares a second time.
- The gr(g, by = f) fixture of test-se-check.R and test-se-tier.R pins
  standard-error verdicts at a point 1.44e-4 below the maximum.
- `fd_hessian()` keeps a step that leaves the matrix positive definite
  however large its rise; fit_imp2 in raw units keeps 0.67.0's 1.5e-3.
- Fuzz specs with `weights(w)` and `ar1()` on beta and Gamma reach
  degenerate objectives near -1e15 on rellib-r8.

## For consolidation

- Lane illse's new `test-se-conv.R` fails 2 assertions on a merge with
  this lane (the review's `merged-test-se-conv.txt`, 15 pass, 2 fail).
  It needs a gaussian `y ~ x` with `iter.max = 1` to stop short, and
  from the least-squares start that stop is the maximum. The
  coordinator has illse change its test.
- The review's other semantic interactions: illse's `newton_polish()`
  steps from coded stops, and a relabeled restart is now code 0; illse's
  at-the-maximum SE judgment runs on the fits that now stop at
  "singular convergence (7)" (the gr-by fixture, two test-sparsex.R
  fits, test-brms-pins.R).

## Not done

- The within-group or anova start for sigma and the group sd: measured
  on the wide design, not built.
- Families outside the list (tweedie, compois, nbinom1, beta_binomial,
  the zero-inflated and hurdle count and beta families, von_mises,
  inverse.gaussian, Gamma's inverse link, cox) keep their starts.
- Nonlinear formulas keep their start (lane postfix's).

## Overlap

- Lane postfix (nonlinear start path, ledger row brmsfit-methods:955):
  nonlinear formulas keep their start here. In R/fit.R postfix edits
  `fit_assembled()` and `fit_recovery_starts()`; this lane edits
  `optimize_obj()` (adjacent to `fit_recovery_starts()`),
  `restart_result()`, `make_start()` and the `start` documentation of
  `frm()`.
- Lane illse (SE check on fits with a convergence warning): this lane
  edits `fd_hessian()` in R/se-check.R, which illse's diff does not
  touch at the time of reading. Its change would make the SE check run
  on the fits that now stop at "singular convergence (7)" (the gr-by
  fixture, two test-sparsex.R fits).
- R/autoscale.R (`predictor_unit_*`, `predictor_resid_sd()`,
  `predictor_col_spread()`), R/par-template.R, NEWS.md,
  frmtmb.sample's test-gr-by-draws.R.

## R CMD check

`R CMD check --as-cran`, built and checked in `dev/optstart-check/` and
`dev/optstart-check-sample/` (`dev/optstart-check.sh`,
`dev/optstart-check-sample.sh`; `R_LIBS` the lane library, rellib-r8,
the user library; logs kept as `dev/optstart-log/check-core.log`,
`check-sample.log` and their tests' `.Rout`, the check trees removed):

- frmtmb: Status: 1 NOTE (`checking HTML version of manual`, V8
  unavailable, the expected one); tests `[ FAIL 0 | WARN 0 | SKIP 388 |
  PASS 13104 ]`.
- frmtmb.sample (its test file changed): Status: OK; tests
  `[ FAIL 0 | WARN 0 | SKIP 80 | PASS 669 ]`.

## Version

A minor bump: breaking changes to every location-scale fit's start and
to the reported convergence code of restarts. No floor changes: the
frmtmb.sample test change accepts the warning of either core.
