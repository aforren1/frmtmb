# Lane rescorpp: joint predictive draws of `set_rescor(TRUE)` models

Round opened 2026-10-28, base ac4906bd (frmtmb 0.77.0, frmtmb.sample
0.24.0). Library `C:/Users/adf44/source/r/wt-rescorpp-lib` (core and
frmtmb.sample from this tree), behind it rellib-r15. brms 2.23.0 from
the user library. Every number below is pasted from the log named
beside it; logs are local (`dev/.gitignore`).

## 1. The defect, seen to fail on rellib-r15

`dev/rescorpp-repro.R` is `dev/review-mointerop-rescorpp.R` with its
library lines pointed at this round (`dev/rescorpp-common.R`). Seed 3,
300 rows, residual correlation 0.8 (`dev/rescorpp-repro-base.log`):

    == posterior_predict(ds) dim ==
    [1] 300 300
    == cor of draws across responses, rows 1..5 ==
    [1] -0.052 -0.040 -0.005  0.082  0.026
    == core predict cor rows 1..5 ==
    [1] 0.795 0.809 0.861 0.831 0.818
    == log_lik(ds) vs bivariate normal: max rel diff ==
    [1] 1.621908e-15

What the record needed to say more exactly: the five correlations
come from two calls, `resp = "y1"` and `resp = "y2"`, and two calls are
independent in brms as well. The defect is that NO call returned the
joint draw: `posterior_predict(ds)` gave the first response alone (300
x 300), and `resp = c("y1", "y2")` was an error. On the lane the same
repro prints `300 300 2` for the shape; its two-call correlations stay
near 0 (0.003 to 0.102, `dev/rescorpp-repro-lane.log`), as brms's
would. The joint call is measured in section 3.

The new test file, run against rellib-r15 with the lane's tests
(`dev/rescorpp-t-rescor-base.log`):

    [ FAIL 13 | WARN 0 | SKIP 0 | PASS 1 ]
    RESULT test-rescor-draws.R pass 1 fail 5 error 8 skip 0 warn 0

The behavioral failures among them, beside the shape failures:
`log_lik(ds, resp = "y1")` of a rescor fit against y1's own normal
density, max abs difference 7.73 (base returned the joint density of
all three responses); the Student-t margin, 9.66;
`posterior_predict(resp = c("y1", "y3"))` stopped with "family '' has
no simulator yet"; `log_lik(resp = c(...))` stopped with "'length = 2'
in coercion to 'logical(1)'".

## 2. What brms 2.23.0 does, measured

Fits: `dev/rescorpp-brmsfit.R` (2 chains, 2000 iterations, seed 21) on
the four data sets of `dev/rescorpp-data.R`: g3 (three gaussian
responses, constant sigma, seed 11, 200 rows), gs (two, `sigma ~ x` on
y1, seed 12), st (two, student, nu 5, seed 13), g2f (y1 and y2 of g3,
`set_rescor(FALSE)`). Measurements `dev/rescorpp-brmsbehav.R` and
`dev/rescorpp-brmsstudent.R`.

- Without `resp`, `posterior_predict()`, `posterior_epred()`,
  `posterior_linpred()` and `predictive_error()` return draws x
  observations x responses, third margin named by the responses;
  `fitted()`, `predict()`, `residuals()` return observations x 4 x
  responses. The same with `set_rescor(FALSE)`.
- `predictive_interval()` without `resp`: "For the default method
  'object' should be a matrix." `pp_check()` without `resp`: "Argument
  'resp' must be a single variable name when applying this method to a
  multivariate model." `conditional_effects()` returns every response
  (`y1.y1_x`, `y2.y2_x`, `y3.y3_x`).
- `log_lik()` is draws x observations. `combine = FALSE` gives the same
  2000 x 200 under rescor, 2000 x 200 x 2 without it.
- `log_lik(resp = "y1")` is y1's own normal density: max abs difference
  0 against `dnorm()` at the draws (g3 and g2f).
- `resp = c("y1", "y2")` of three: 2 layers, drawn jointly. Mean over
  rows of the residual correlation 0.7488.
- brms's own predictive residual correlation (prediction minus
  expectation, mean over rows) against its rescor posterior mean: g3
  0.7492 / 0.7491, 0.5508 / 0.5508, 0.3110 / 0.3118; gs 0.7226 /
  0.7224; st 0.6901 / 0.6883; g2f 0.0009 / 0.
- brms is WRONG on one response of a student rescor model. Its
  one-response preparation has `nu` 1 in every draw (`str(nu)`: `num
  [1:2000] 1 1 1 ...`, where the draws of `nu` are 8.64, 5.41, ...).
  `log_lik(resp = "y1")` equals `dt(z, 1)` to max abs difference 0 and
  differs from `dt(z, nu)` by 6.44. `posterior_predict(resp = "y1")` is
  a Cauchy draw: |prediction minus expectation| at the 50, 90 and 99
  percent quantiles is 1.09, 6.95 and 68.8 against 0.78, 2.15 and 4.19
  for y1 of the joint draw. frmtmb uses the shared `nu`, the margin of
  the multivariate t (project tiebreaker: brms is clearly wrong here).
  Not posted upstream.
- brms refuses a permuted `resp` under rescor. `log_lik(resp = c("y3",
  "y1"))` stops with "The following variables are missing in the draws
  object: {'rescor__y3__y1'}" (`dev/rescorpp-compare-g3.log`). Corrected
  after the review (its finding 8): not on log_lik alone. On g3, gs, st
  and re, brms 2.23.0 stops the same way ("{'rescor__y2__y1'}") in
  every method that prepares several responses: `posterior_predict`,
  `posterior_epred` (which reads no rescor), `posterior_linpred`,
  `predictive_error`, `fitted`, `predict`, `residuals`, `log_lik` and
  `loo_predict` (`dev/review-rescorpp-layout-*.log`). Under
  `set_rescor(FALSE)` it answers. The cause is `get_cornames(resp, ...)`
  in `prepare_predictions.mvbrmsframe()`, which builds the correlation
  names in the order given. frmtmb answers a permuted `resp` in every
  method, with the correlations in the order named; frmtmb.sample's
  NEWS now says so.

## 3. The fix

frmtmb.sample `posterior_predict()` (`R/methods-draws.R`):

- `resp` is a set: `NULL` is every response, several names give their
  layers in that order, either spelling (`y_a` or brms's `ya`).
- One response: the path it had, unchanged (`draws_pp_setup()` +
  `draws_pp_one()` are the old body split in two).
- Several, `set_rescor(TRUE)`: `draws_pp_joint()`. Per draw, one
  `draws_fit_at()`, one `eval_dpars()` for every response, then core's
  `predict_rescor_draw()`, which `predict()` on the fit already used:
  `Z %*% chol(R)` with `Z` an n x K standard normal matrix, each column
  scaled by that row's `sigma` and shifted by `mu`; for student, every
  row divided by `sqrt(W / nu)`, one `W ~ chi^2(nu)` per row shared by
  the responses. `chol(D R D) = chol(R) D` for a positive diagonal `D`,
  so this IS brms's `rmulti_normal(mu, Sigma_i)` (`draws %*%
  chol(Sigma)`) row by row, with ONE Cholesky factor per draw even when
  `sigma ~ x` makes `Sigma` vary by row. brms factors `Sigma` once per
  row and draw.
- Several, without rescor: each response's own draw, stacked. The
  first layer is bit for bit the draw `resp = "y1"` gives from the same
  stream (test below), which is what 0.24.0 returned without `resp`.

The other readers:

- `posterior_epred()`, `posterior_linpred()`: brms's multivariate
  layout without `resp` or with several (each response computed as
  before and stacked; the stack refuses a response that is not one
  value per row, in the words `loo_draws_expect()` already used, which
  now calls the same `draws_resp_stack()`).
- `predictive_error()`: several responses through ONE
  `posterior_predict()` call, so the errors of a rescor model carry its
  correlation; `y - yrep` per layer.
- `fitted()`, `predict()`, `residuals()`: their summaries of those, n x
  4 x K. `residuals(type = "pearson")` divides by the predictive sd per
  row and response.
- `predictive_interval()` without `resp` on a multivariate model:
  refused before any draw, as brms stops.
- `log_lik()` (`R/loo.R`): under rescor, `resp` of length 1 is that
  response's own density (the family's, through `row_lpdf()`, with the
  shared `nu` for student); several are their joint density, through
  core's `rescor_row_loglik(fit, dpv, resp)`. `combine = FALSE` returns
  the joint matrix under rescor. `resp` may name several responses.
- Core: `rescor_row_loglik()` takes `resp` and reorders or subsets the
  correlation matrix by it (`R/predict-brms.R`); `predict_rescor_draw()`
  joins the sampling API (`R/sampling-api.R`, one line in the export
  list, its own line).
- Also fixed, because the roxygen run reported it on a file this lane
  changes: `draws_refuse_newdata()`'s doc block had been separated from
  its function by `draws_pp_structured()`, so `@noRd` was followed by
  text ("methods-draws.R:2880: @noRd must not be followed by any text",
  a pkgcheck failure). The two blocks are swapped back; no text changed.

## 4. Evidence on the lane

### 4.1 Against brms, residual correlation of the predictive draws

`dev/rescorpp-compare.R` per data set (`dev/rescorpp-compare-<w>.log`):
frmtmb draws from `frm_sample(chains = 2, iter = 2000, seed = 21)`
(2000 draws, 200 rows); predictive seed 101, brms seed 102. Three
quantities per pair, each the mean over rows of the across-draw
correlation of prediction minus expectation:

- frmtmb: `posterior_predict(ds)`.
- brms code on frmtmb's draws: `brms:::prepare_predictions()` built
  from the brms fit with frmtmb's draws substituted column by column
  (every column brms reads has a frmtmb counterpart: "brms variables
  read, absent from frmtmb's draws: character(0)"), then
  `rstantools::posterior_predict()` of that prep, so brms's own
  `posterior_predict_gaussian_mv()` / `_student_mv()` runs at the SAME
  parameter draws.
- target: given the draws, `mean(w r s_j s_k) / sqrt(mean(w s_j^2)
  mean(w s_k^2))` per row, `w = nu / (nu - 2)` for student and 1
  otherwise.

z is the mean over rows of the difference in units of its own standard
error (sd over rows / sqrt(200)); ratio is mean_a / mean_b.

    g3 pair 12: frmtmb vs target      0.7495 / 0.7505  z -1.289
    g3 pair 12: brms code vs target   0.7499 / 0.7505  z -0.846
    g3 pair 12: frmtmb vs brms code   0.7495 / 0.7499  z -0.380
    g3 pair 13: frmtmb vs target      0.5542 / 0.5542  z  0.001
    g3 pair 13: brms code vs target   0.5532 / 0.5542  z -0.980
    g3 pair 13: frmtmb vs brms code   0.5542 / 0.5532  z  0.631
    g3 pair 23: frmtmb vs target      0.3147 / 0.3149  z -0.133
    g3 pair 23: brms code vs target   0.3147 / 0.3149  z -0.144
    g3 pair 23: frmtmb vs brms code   0.3147 / 0.3147  z  0.007
    gs pair 12: frmtmb vs target      0.7212 / 0.7204  z  1.053
    gs pair 12: brms code vs target   0.7207 / 0.7204  z  0.379
    gs pair 12: frmtmb vs brms code   0.7212 / 0.7207  z  0.470
    st pair 12: frmtmb vs target      0.6909 / 0.6891  z  1.361
    st pair 12: brms code vs target   0.6896 / 0.6891  z  0.382
    st pair 12: frmtmb vs brms code   0.6909 / 0.6896  z  0.761
    st |e|^2 Spearman, frmtmb vs brms 0.4010 / 0.4024  z -0.741
    g2f pair 12: frmtmb vs target     -0.0011 / 0      z -0.728
    g2f pair 12: brms code vs target   0.0010 / 0      z  0.645

The student line is the shared chi-square: the Spearman correlation of
the squared residuals is 0.40 in both, where independent chi-squares
would leave only what the 0.69 correlation gives.

The spread of each response and of a sum, against their expectations
given the draws (`sd(e_j)` per row against `sqrt(mean(w s_j^2))`, and
`Var(e1 + e2)` against `mean(w (s1^2 + s2^2 + 2 r s1 s2))`):

    g3  resp 1 0.9767 / 0.9782 z -1.446   resp 2 1.9987 / 2.0022 z -1.481
        resp 3 0.4948 / 0.4959 z -1.836   Var(e1+e2) 7.877 / 7.905 z -1.528
    gs  resp 1 1.0241 / 1.0246 z -0.449   resp 2 0.9834 / 0.9843 z -0.848
        Var(e1+e2) 3.633 / 3.636 z -0.403
    st  resp 1 1.3664 / 1.3695 z -1.087   resp 2 2.1335 / 2.1347 z -0.298
        Var(e1+e2) 10.458 / 10.462 z -0.087
    g2f resp 1 0.9803 / 0.9828 z -2.346   resp 2 2.0011 / 2.0019 z -0.357
        Var(e1+e2) 4.962 / 4.974 z -0.995

Every sd z above is negative, the g2f one (the per-response simulator
this lane did not change) included, so I checked that it is one
predictive seed and not a bias: `dev/rescorpp-sdrep.R`
(`dev/rescorpp-sdrep.log`) holds the g3 draws fixed and repeats the
draw over predictive seeds 1 to 30. Mean z over the 30: 0.158, 0.091
and 0.210 for the three responses of the joint draw, 0.008 for
`resp = "y1"`; sd 0.89 to 1.22. A pure-R control (`rnorm()` at the
draws' sigma) gives mean -0.039, sd 0.914. A null simulation of the
statistic with no frmtmb code (`dev/rescorpp-sdstat.R`, seed 7, 300
replicates) gives mean -0.10, sd 1.02. The single-seed negatives were
chance.

Core `predict()` on the ML fit, for the record (seed 103, 2000
replicates): mean residual correlation 0.7562 against
`rescor_matrix()` 0.7572 (g3), 0.7267 / 0.7282 (gs), 0.6908 / 0.6921
(st). With `resp = c("y1", "y3")` (`dev/rescorpp-audit-lane.log`, seed
7, 3000 replicates): 0.5574 against 0.5601.

### 4.2 `log_lik()` against brms at the same draws

Same transplant, `rstantools::log_lik()` of the prep. Control first:
the route with brms's OWN draws against `log_lik(b)`, max abs
difference 0 on all four. Then frmtmb's `log_lik(ds)` against brms's
code at frmtmb's draws, max abs difference over max |ref|:

    g3  joint 8.32e-16   resp = "y1" 0   resp = c("y1", "y3") 4.81e-16
    gs  joint 8.78e-16   resp = "y1" 0
    st  joint 5.11e-16   resp = "y1" 1.33 (brms's nu = 1; section 2)
    g2f joint 0          resp = "y1" 0

and for st, frmtmb's `resp = "y1"` against `dt(z, nu) - log(sigma)` at
the draws: 2.97e-16.

### 4.3 The audit of every other reader

`dev/rescorpp-audit.R`, both arms (`dev/rescorpp-audit-<arm>.log`; g3
and st, 2 chains of 1000, seed 21). "Joint as brms" means: reads the
joint law where brms does, and a margin where brms does.

Each reader, with rellib-r15 first and the lane second:

- `posterior_predict()`: the first response only, and no call returned
  the joint draw. Lane: joint, in brms's layout.
- `posterior_epred()`, `posterior_linpred()`: the first response only.
  Lane: every response; each value unchanged.
- `predictive_error()`: the first response. Lane: joint, from one
  predictive call.
- `fitted()`, `predict()`, `residuals()`: the first response. Lane: n x
  4 x K.
- `predictive_interval()`: the first response. Lane: refused without
  `resp`, as brms stops.
- `log_lik()`: the joint density, for every `resp`. Lane: the joint
  density; with `resp`, a margin or the joint density of a subset.
- `log_lik(combine = FALSE)`: refused. Lane: the joint matrix, as brms.
- `loo(resp = "y1")` elpd: -724.69, the joint. Lane: -280.79, y1's own.
- `loo()` and `waic()` without `resp`: -724.69 and -724.65 on both.
- `loo_predict()`, `loo_epred()`: 3D, with weights from the joint
  density for every `resp`. Lane: 3D, with weights by `resp`, as brms's
  `loo(resp =)` computes them.
- `loo_R2()`: weights from the joint density. Lane: weights from
  `log_lik(resp = r)`, as brms's `loo_R2()` computes them.
- `bayes_R2()`: per response, from the expectation only. Unchanged.
- `pp_check()`: needs `resp`, reads the margin. Unchanged.
- `conditional_effects(method = "posterior_predict")`: the margin, with
  the shared `nu`; the first response only. Unchanged.
- Core `predict()`: joint, and subsets by `resp`. Unchanged.
- Core `simulate()`: refused, "simulate() is not supported yet for
  multivariate fits". Unchanged.
- Core `frm_simulate()`: refused, "frm_simulate() supports univariate
  models". Unchanged.
- `mi()` and `subset()` with rescor: refused at the fit ("mi() cannot be
  combined with rescor = TRUE"; "subset() cannot be combined with
  rescor = TRUE: ..."). Unchanged.

`loo_R2(ds3)` moved by up to 2.6e-4 (R2y2 0.03143 to 0.03117), from the
weights. The conditional-effects band at x = 0 against
`posterior_predict()`'s own quantiles there: g3 (-0.93, 2.99) against
(-0.77, 2.95); st (-1.60, 3.89) against (-1.65, 3.85). These read the
margin and did not change between arms (identical numbers).

### 4.4 Guard absent: `set_rescor(FALSE)`

Without rescor the joint path is not taken: `draws_pp_one()` runs per
response, exactly the old loop. The test "without rescor each layer is
the response's own draw" asserts `identical(pp[, , 1], posterior_predict(
ds, resp = "y1"))` from the same seed, and that the cross-response
residual correlation is within 4.5 standard errors of 0 (g2f in 4.1:
z -0.73). The one-response path is the old code split in two
functions with no change of call or order. The rescor = FALSE outputs
that changed are shapes only: without `resp` they now carry every
response.

## 5. Cost

`dev/rescorpp-timing.R`: 3 gaussian responses, 500 rows, 4000 draws
(`frm_sample(chains = 2, iter = 4000, seed = 5)`). Arms interleaved in
one process, 3 rounds, minimum; a fixed `crossprod()` of a 2000 x 2000
matrix is the load control. Run alone on the machine except for other
lanes' processes (`dev/rescorpp-timing2-<arm>.log`):

    lane:  control 2.22  independent 4.25  joint 2.08   (seconds, min)
           joint / independent 0.489
    base:  control 2.22  independent 4.25

"independent" is the three one-response calls, the only route on
rellib-r15; it costs the same on both arms, so the split of the old body
into `draws_pp_setup()` and `draws_pp_one()` added nothing measurable.
The joint call costs half of three one-response calls because it runs
one `draws_fit_at()` and one `eval_dpars()` per draw for all three
responses. A first version called `eval_dpars()` once per response
(it evaluates every response each time): its profile put 67 percent of
the joint call in `eval_dpars()`, and the measured joint was 0.88 of
independent under load (`dev/rescorpp-timing-lane.log`, a busier run).

Profile of the joint call (`Rprof`, 5 ms): `eval_dpars()` 52 percent,
`predict_rescor_draw()` 25 percent, `draws_fit_at()` 10 percent. Per
call: `predict_rescor_draw()` 0.13 ms (one 3 x 3 Cholesky factor, one
500 x 3 normal matrix and its product), `eval_dpars()` 0.22 ms,
`draws_fit_at()` 0.07 ms. The draw allocates the n x K normal matrix
and K result vectors per draw, which is the size of the output; there
is no factorization per row. The profiled call took 0.61 s against the
timed 2.08 s minimum; the difference is the timed loop's garbage
collection of the arms before it, which the profile does not see.

## 6. Tests

New: `extensions/frmtmb.sample/tests/testthat/test-rescor-draws.R`, 9
blocks, 51 expectations; one block added to core's
`tests/testthat/test-mv-gaps.R`. Seen to fail on rellib-r15: section 1
(sample file, 13 failures) and `dev/rescorpp-t-mvgaps-base.log` (core
block: "unused argument (c("y1", "y3"))", the weak form; its
behavioral counterpart is the sample file's `log_lik(resp =)` block).
On the lane: `RESULT test-rescor-draws.R pass 51 fail 0 error 0 skip 0
warn 0`.

The z bounds (4.5) are on a Monte Carlo error each test measures from
the spread over rows; the old independent draw sits at z above 20 on
the same statistic, which the test asserts too, so the bound
discriminates. Exactness tests use a multiple of `.Machine$double.eps`
times the largest reference value.

All nine suites, one file per process, `NOT_CRAN=true`,
`FRMTMB_BRMS_FIT_TESTS=true`, every `lib:` line the lane library
(`dev/rescorpp-tests/all.log.lib`, 0 lines otherwise):

dev/rescorpp-tests/all.log (all nine suites, gated, one file per process):
    frmtmb           files 285 pass 22811 fail 0 error 0 skip 15 warn 0
    frmtmb.coupling  files  11 pass   551 fail 0 error 0 skip  5 warn 0
    frmtmb.eam       files  32 pass  1844 fail 0 error 0 skip  3 warn 0
    frmtmb.latent    files  12 pass   681 fail 0 error 0 skip  2 warn 0
    frmtmb.learn     files  15 pass   501 fail 0 error 0 skip  2 warn 0
    frmtmb.ode       files  11 pass   549 fail 0 error 0 skip  1 warn 0
    frmtmb.sample    files  74 pass  3945 fail 2 error 0 skip  1 warn 0
    frmtmb.spatial   files   8 pass   561 fail 0 error 0 skip  0 warn 0
    frmtmb.survival  files  15 pass   552 fail 0 error 0 skip  0 warn 0
    all ran 463 of 463
dev/rescorpp-tests/sample2.log (rerun after the last sample change):
    frmtmb           files  27 pass  2875 fail 0 error 0 skip  0 warn 0
    frmtmb.eam       files   3 pass   136 fail 0 error 0 skip  0 warn 0
    frmtmb.latent    files   1 pass    13 fail 0 error 0 skip  0 warn 0
    frmtmb.ode       files   1 pass    15 fail 0 error 0 skip  0 warn 0
    frmtmb.sample    files  74 pass  3946 fail 0 error 0 skip  1 warn 0
    frmtmb.spatial   files   2 pass   258 fail 0 error 0 skip  0 warn 0
    sample2 ran 108 of 108

The 2 failures of the first run were `test-message-uniqueness.R`: two
refusal templates I had copied ("Invalid argument 'resp'..." and "...
stacks one matrix per response..."). Both now come from one function
each (`draws_refuse_resp()`, `draws_resp_stack()`, which
`loo_draws_expect()` now calls too). The skips are other tiers' gates,
none reached by this change: `FRMTMB_SCALE_TESTS` (the `test-scale.R`
files and `test-curve-scale.R`), `FRMTMB_FUZZ` (`test-fuzz.R`) and
`FRMTMB_DRMTMB_FIT_TESTS` (`test-drmtmb-agreement.R`, 13).

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
`dev/rescorpp-ob-*.txt`), the files that touch the changed numerics:

    RESULT test-mv-gaps.R pass 62 fail 0 error 0 skip 0 warn 0 seconds 1.9 
    RESULT test-draws-methods.R pass 174 fail 0 error 0 skip 0 warn 0 seconds 12.9 
    RESULT test-loo-predict.R pass 153 fail 0 error 0 skip 0 warn 0 seconds 157.4 
    RESULT test-loo.R pass 146 fail 0 error 0 skip 0 warn 0 seconds 142.5 
    RESULT test-rescor-draws.R pass 51 fail 0 error 0 skip 0 warn 0 seconds 18.7 

The ported bin-2 tiers on the lane build (`FRMTMB_PORT_LIB` the lane
library, `R_LIBS_SITE` rellib-r15; `dev/rescorpp-bin2-*.txt`):

    RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=105 fail=0 err=0 skip=0 blocks=17
    RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=65 fail=0 err=0 skip=0 blocks=16

`dev/release/run-sessions.sh` with `FRMTMB_RELLIB` the lane library,
all nine suites in forward, reverse and preload order, gated
(`dev/rescorpp-tests/sessions/rescorpp-sess.log`; core and
frmtmb.sample from the lane library in all 27):

    SESSION frmtmb forward files=285/285 pass=22811 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb reverse files=285/285 pass=22961 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb preload files=285/285 pass=23054 fail=0 err=0 skip=17 warn=0
    SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam forward files=32/32 pass=1844 fail=0 err=0 skip=3 warn=0
    SESSION frmtmb.eam reverse files=32/32 pass=1844 fail=0 err=0 skip=3 warn=0
    SESSION frmtmb.eam preload files=32/32 pass=1844 fail=0 err=0 skip=3 warn=0
    SESSION frmtmb.latent forward files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent preload files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn forward files=15/15 pass=501 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn reverse files=15/15 pass=501 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn preload files=15/15 pass=501 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample forward files=74/74 pass=3946 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample reverse files=74/74 pass=3946 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample preload files=74/74 pass=3946 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=561 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=561 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=561 fail=0 err=0 skip=0 warn=0
    RESCORPP-SESS ran 27 of 27

Every LEAK line is tagged ([with load of frmtmb.sample], [with load
of a Stan model (rstan)]); untagged LEAK lines: 0.

`R CMD check --as-cran` once per changed package, built and checked
in `dev/rescorpp-check/` (`dev/rescorpp-check.sh`;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH,
`R_LIBS` the lane library, rellib-r15, the user library):

    frmtmb         Status: 1 NOTE  (HTML manual: V8 unavailable)
    frmtmb.sample  Status: 1 NOTE  (HTML manual: V8 unavailable)

The first core build failed while the sample check was running in
the same `dev/rescorpp-check/`: `R CMD build` copies the tree before
it applies `.Rbuildignore`, and a file the running check held could
not be copied. Run one after the other, both passed.

## 7. Ledger

No verdict flips. `posterior_predict:69` and `:72` pass through core's
`predict()` (their note) and stay pass; their note's clause "which
draws the responses independently (dev/test-backlog.md, filed at the
0.77.0 release)" is now false and should be reworded at consolidation,
for example "frmtmb.sample's posterior_predict() draws them jointly
since this round". `log_lik:107` and `:109` call
`rescor_row_loglik(fit, dpv)` through the bin-2 adapter; the new `resp`
argument defaults to every response, so their values are unchanged
(section 6 has the bin-2 runs).

## 8. Not done, and why

- `log_lik(combine = FALSE)` without rescor is still refused; brms
  returns a draws x observations x responses array. Not this lane's
  defect; it is one `abind` away.
- The bin-2 ledger is not regenerated (brief).

### Filed, not fixed (for `dev/test-backlog.md` at consolidation)

- **`conditional_effects()` on multivariate draws** (pre-existing, the
  same on rellib-r15 and the lane; review finding 7,
  `dev/review-rescorpp-readers-{base,lane}.log`). Without `resp` it
  returns one effect, for the first response; brms returns one per
  response (`y1.y1_x`, `y2.y2_x`, `y3.y3_x`,
  `dev/rescorpp-brmsbehav.log`). With several responses it fails with
  a message that is not about the call:
  - `conditional_effects(ds, resp = c("y1", "y2"))`: "'length = 2' in
    coercion to 'logical(1)'";
  - the same with `method = "posterior_predict"`: "method = 'predict'
    needs a family with a simulator", false for a gaussian model.
  Fix: one panel set per response, as brms; until then, refuse several
  `resp` with a sentence that says so.
- **A stored `loo_R2` in saved 0.77.0 draws keeps the old weights**
  (review finding 9, `dev/review-rescorpp-readers-lane.log`). Draws
  saved on rellib-r15 after `add_criterion(c("loo", "loo_R2"))` return
  their stored 0.77.0 values on the lane (`R2y2` -0.002062348),
  computed with joint-density weights for `resp`. frmtmb has no users,
  so no seam; mention it in the release record.
- **One-response matrices carry `NULL` dimnames where brms has
  `list(NULL, NULL)`** (pre-existing nit, review finding 12): every
  `resp = first` difference on `posterior_predict`, `predictive_error`,
  `predict(summary = FALSE)` and `log_lik` in
  `dev/review-rescorpp-layout-*.log` is this and nothing else.
- **For `dev/upstream-bugs.md`, two brms 2.23.0 defects** (no issue
  found by `gh search issues --repo paul-buerkner/brms`, 2026-10-28):
  1. One response of a student rescor model reads `nu = 1`.
     `brmsterms.R` lines 184 to 189 set `pfix$nu <- 1` on every response
     when `no_nu()` is true (rescor + student); the shared `nu` reaches
     the prep only when `length(resp) > 1`
     (`prepare_predictions.mvbrmsframe()`). So
     `posterior_predict(resp = "y1")` is a Cauchy draw and
     `log_lik(resp = "y1")` is `dt(z, 1)` (max abs difference 0; 6.44
     from `dt(z, nu)`; `dev/rescorpp-brmsstudent.log`).
  2. A `resp` out of model order under rescor stops every method that
     prepares several responses, with "The following variables are
     missing in the draws object: {'rescor__y2__y1'}", because
     `get_cornames(resp, ...)` builds the names in the order given
     (section 2).

## 9. Overlaps with the other lanes

- Core `R/sampling-api.R`: one line added to the `@rawNamespace
  export(...)` list (`predict_rescor_draw,`, on its own line after
  `rescor_row_loglik,`), one `@aliases` line and one paragraph;
  `NAMESPACE` and `man/frmtmb-sampling-api.Rd` regenerate from it. A
  lane that adds an export there merges line by line.
- Core `R/predict-brms.R`: `rescor_row_loglik()` only.
- Core `NEWS.md`: a new `# frmtmb (development version)` section with
  a "Sampling API" subsection; the meme lane will add its own section
  under the same heading.
- frmtmb.sample `R/methods-draws.R` (posterior_epred, posterior_linpred,
  posterior_predict, predictive_interval, predictive_error, predict,
  residuals), `R/loo.R` (log_lik), `R/loo-predict.R`
  (`loo_draws_expect()` now calls `draws_resp_stack()`), `R/zzz.R`
  (the rescor x frm_sample compat note), `NEWS.md`. No lane of this
  round named these files; meme works in `R/priors.R` and the frame.
- equiv may flag rescor pairs: posterior_predict() of rescor draws now
  returns draws x observations x responses without `resp`, so a
  comparison written against the 0.24.0 shape (first response only)
  must take `[, , 1]` or pass `resp`.

## 10. Version and floor

frmtmb: a patch bump would do for the code (one export, one argument
with a default), but the round's other lanes decide it. frmtmb.sample:
a minor bump (BREAKING changes in output shape and in `log_lik(resp =)`
values). frmtmb.sample's `Imports: frmtmb (>= ...)` must rise to the
core version that carries this round: it calls `predict_rescor_draw()`
and `rescor_row_loglik(resp =)`, which 0.77.0 lacks. I did not edit
DESCRIPTION (no version numbers, brief).

## 11. Files

Changed: `R/predict-brms.R`, `R/sampling-api.R`, `NAMESPACE`,
`man/frmtmb-sampling-api.Rd`, `NEWS.md`, `tests/testthat/test-mv-gaps.R`;
`extensions/frmtmb.sample/R/methods-draws.R`, `R/loo.R`,
`R/loo-predict.R`, `R/zzz.R`, `man/posterior_epred.Rd`,
`man/sample-log_lik.Rd`, `man/sample-loo_predict.Rd` (fix round 1),
`NEWS.md`. New:
`extensions/frmtmb.sample/tests/testthat/test-rescor-draws.R` and the
`dev/rescorpp-*` scripts and this record. Nothing committed.

## 12. Fix round 1 (review `dev/reviews/2026-10-28-rescorpp.md`, FIX)

What changed, by the review's finding:

1. **Major: `predict()` summarized a category response's codes**
   (finding 1). Without `resp`, or with several, `predict()` on an
   ordinal + gaussian model returned the mean of the category codes
   for the ordinal layer, where `predict(resp = "o")` gives category
   shares; brms stops on both. `predict.frmtmb_draws()` now refuses
   before drawing, when `summary = TRUE` and a response it would stack
   is category valued ("Response 'o' is a category: call predict()
   with resp = for it alone, or use summary = FALSE for the draws").
   `summary = FALSE` still returns the draws, as brms does. Measured on
   the ordinal + gaussian model of the tests (`dev/rescorpp-ordcheck.R`,
   `dev/rescorpp-ordcheck.log`): `predict()` refused, `summary = FALSE`
   5 x 120 x 2, `posterior_linpred()` 5 x 120 x 2 (an ordinal `mu` is
   one value per row), `posterior_epred()` and `fitted()` refused by
   `draws_resp_stack()`, `residuals()` refused as for any ordinal
   family.
2. **The Student-t test** (finding 2). It now also checks each margin's
   tails: the residual standardized by its draw's sigma is beyond the
   draw's `qt(0.995, nu)` in 1 percent of the cells, a binomial count
   over 600 x 120 independent cells, |z| < 4.5. A normal draw puts
   about 0.03 percent there at nu near 6.
3. **Permutations of three responses** (finding 3). `log_lik(resp =
   c("y3", "y1", "y2"))` against `mvtnorm::dmvnorm` written in that
   order, `rescor_row_loglik(fit, dp, c("y3", "y1", "y2"))` in core, and
   `posterior_predict(resp = c("y3", "y1", "y2"))` pair by pair against
   the correlation the draws imply. The test comment that claimed a
   two-response permutation checks the order is corrected.
4. **`log_lik(resp =)`** (finding 4) reads `resp` through
   `draws_resp_set()` on a multivariate model, as the draw methods do:
   brms's `ya`, a repeated name refused, an empty `resp` every
   response. `draws_resp_set()` now refuses a repeated name for every
   method (it dropped it silently before), and every `resp` reader maps
   `character(0)` to every response. `loo()`'s `yhash` reads the same
   spelling, so `loo(resp = "ya")` and `loo(resp = "y_a")` are
   identical.
5. **`predict_rescor_draw()` pairs `dps` with `rspecs` by name** when
   `dps` has names, and refuses a missing name or, unnamed, a length
   that differs (finding 5); the sampling-API page says so. The
   reviewer's subset call (`rspecs[c("y3", "y1")]` with the whole
   `eval_dpars()`) is the core test.
6. **`loo_predict()`, `loo_epred()`, `loo_linpred()` and
   `loo_predictive_interval()` take several responses** (finding 6):
   their layers in the order named, weighted by `log_lik(resp = <the
   set>)`, the joint density under rescor. The test checks the layers
   and that the weights equal `draws_psis_of(log_lik(resp = set))`.
7. **NEWS** (findings 8 and 11): one bullet lists every method that
   takes several responses, the two spellings, the repeated-name
   refusal and the empty `resp`, and says that a permuted `resp` is
   answered where brms 2.23.0 stops. Another says that `posterior_epred()`
   and `fitted()` without `resp` refuse a category response where brms
   binds its category columns beside the other responses, and that
   `predict()` refuses to summarize it, as brms stops. The Rd section
   "Multivariate models" says the same.
8. **Filed** (findings 7, 9, 12, and the brms defects): section 8,
   "Filed, not fixed". Section 2 now describes the permuted-`resp`
   brms defect across every method (finding 8).
9. The floor bump is left to consolidation (finding 10, section 10).

### Seen to fail

On rellib-r15, the lane's test files (`dev/rescorpp-f1-*-base.log`),
then on the lane (`dev/rescorpp-f1-*-lane.log`):

    test-mv-gaps.R       base: pass 59 fail 1 error 2   lane: pass 67 fail 0 error 0
    test-rescor-draws.R  base: pass 1 fail 6 error 11   lane: pass 74 fail 0 error 0

On rellib-r15: `log_lik(resp = "ya")` stops with "names no response
of this model"; `predict(ds)` of the ordinal model answers (the first
response alone) where the test wants a refusal; `loo_epred(resp =
c(...))` stops with "`resp` must name one response"; the name-pairing
block gives a y3 layer that differs from the by-position draw.

The defects new in the lane's first build cannot fail on rellib-r15,
so each was reverted in the loaded namespace by `dev/rescorpp-mutate.R`
(after the reviewer's `dev/review-rescorpp-mutate.R`) and the test
file run (`dev/rescorpp-mut-<mutation>-<file>.log`):

    nonu     (student draw as normal)            test-rescor-draws.R pass 72 fail 2 error 0
    sortC    (log_lik C in model order)          test-rescor-draws.R pass 73 fail 1 error 0
    sortC                                         test-mv-gaps.R      pass 66 fail 1 error 0
    sortR    (draw's R in model order)           test-rescor-draws.R pass 71 fail 3 error 0
    pospair  (dps paired by position)            test-mv-gaps.R      pass 66 fail 1 error 0
    nocatref (no category refusal in predict())  test-rescor-draws.R pass 72 fail 2 error 0
    nodup    (repeated name passed through)      test-rescor-draws.R pass 68 fail 2 error 1

Before this round, "nonu" and "sortC" passed the whole file (the
reviewer's `review-rescorpp-mut-nonu.log`: pass 51 fail 0;
`review-rescorpp-mut-sortC-*.log`: pass 62 and 51, fail 0).

### Reruns

Changed and neighboring files, one per process, gated (every core
`test-*.R` that names rescor, the frmtmb.sample suite, and every file
of the other packages that loads frmtmb.sample;
`dev/rescorpp-tests/fix1.log`, every `lib:` line the lane library):

    frmtmb           files  68 pass  7620 fail 0 error 0 skip 13 warn 0
    frmtmb.eam       files   3 pass   136 fail 0 error 0 skip 0 warn 0
    frmtmb.latent    files   1 pass    13 fail 0 error 0 skip 0 warn 0
    frmtmb.ode       files   1 pass    15 fail 0 error 0 skip 0 warn 0
    frmtmb.sample    files  74 pass  3969 fail 0 error 0 skip 1 warn 0
    frmtmb.spatial   files   2 pass   258 fail 0 error 0 skip 0 warn 0
    fix1 ran 149 of 149

The skips are `FRMTMB_DRMTMB_FIT_TESTS` (13) and `FRMTMB_SCALE_TESTS`
(1), other tiers' gates.

OpenBLAS 0.3.26 with LAPACK (`dev/rescorpp-ob2-*.txt`):

    RESULT test-mv-gaps.R pass 67 fail 0 error 0 skip 0 warn 0 seconds 1.8 
    RESULT test-loo-predict.R pass 153 fail 0 error 0 skip 0 warn 0 seconds 141.6 
    RESULT test-loo.R pass 146 fail 0 error 0 skip 0 warn 0 seconds 129.3 
    RESULT test-rescor-draws.R pass 74 fail 0 error 0 skip 0 warn 0 seconds 23.2 

One-session runs, core and frmtmb.sample, three orders
(`dev/rescorpp-tests/sessions2/rescorpp-sess2.log`):

    SESSION frmtmb forward files=285/285 pass=22816 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb reverse files=285/285 pass=22966 fail=0 err=0 skip=15 warn=0
    SESSION frmtmb preload files=285/285 pass=23059 fail=0 err=0 skip=17 warn=0
    SESSION frmtmb.sample forward files=74/74 pass=3969 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample reverse files=74/74 pass=3969 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample preload files=74/74 pass=3969 fail=0 err=0 skip=1 warn=0
    RESCORPP-SESS2 ran 6 of 6

Untagged LEAK lines: 0; core and frmtmb.sample
from the lane library in all six.

The reviewer's other two mutations on this build
(`dev/review-rescorpp-mutate.R`, `dev/rescorpp-mut-{nochol,sepW}-rescor.log`):

    RESULT nochol test-rescor-draws.R pass 61 fail 13 error 0 skip 0 
    RESULT sepW test-rescor-draws.R pass 73 fail 1 error 0 skip 0 

`R CMD check --as-cran`, run one after the other this time
(`dev/rescorpp-check.sh`; `dev/rescorpp-check-{sample,core}-full.txt`):

    frmtmb.sample  Status: 1 NOTE  (HTML manual: V8 unavailable)
    frmtmb         Status: 1 NOTE  (HTML manual: V8 unavailable)
