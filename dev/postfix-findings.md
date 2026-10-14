# Lane postfix: the open post-processing gaps of the brms port

Date 2026-10-13. Worktree `frmtmb-wt-postfix` (branch `wt-postfix`),
base d2ef6312 (frmtmb 0.70.0, frmtmb.sample 0.18.0). Lane library
`C:/Users/adf44/source/r/wt-postfix-lib` (frmtmb and frmtmb.sample from
this worktree), base `C:/Users/adf44/source/r/rellib-r8`. brms 2.23.0,
rstan 2.32.7, mice and emmeans from the user library, R 4.6.1.

STATUS: FINAL after fix round 2 (see "Fix round 2" and "Status" at the
end). Fix round 2 removed the scan of item f and fix round 1's
acceptance rule: read item f and fix round 1's item 1 as history.

Every "before" below is rellib-r8 and every "after" is the lane
library in front of it. The logs are under `dev/postfix-log/` and
`dev/postfix-tests/` (local, gitignored); the lines quoted here are
pasted from them.

## What was open, re-checked against rellib-r8

`dev/postfix-probe.R before` (log `dev/postfix-log/probe-before.txt`)
calls every item of the brief on 0.70.0. Of vigport's ranked list,
items 1 (`plot()` arguments), 3 (`add_criterion()`) and 5
(group-level `mo()`) were done in 0.68.0 to 0.70.0, and so were three
parts of item 2: `fixef()`, `summary()` and `conditional_effects()` of a
`frm_multiple()` result (NEWS 0.68.0). `summary(waic = TRUE)` was
already refused by name (NEWS 0.68.0). What was still open on 0.70.0:

```text
coef                                           OK NULL [0]
fitted                                         OK NULL [0]
residuals                                      OK NULL [0]
ranef                                          ERROR no applicable method for 'ranef' applied to an object of class "frmtmb_multiple"
VarCorr                                        ERROR no applicable method for 'VarCorr' applied to an object of class "frmtmb_multiple"
vcov                                           ERROR no applicable method for 'vcov' applied to an object of class "frmtmb_multiple"
predict                                        ERROR no applicable method for 'predict' applied to an object of class "frmtmb_multiple"
emmeans                                        ERROR Can't handle an object of class  "frmtmb_multiple"   Use help("models", package = "emmeans") for information on supported models.
logLik                                         ERROR no applicable method for 'logLik' applied to an object of class "frmtmb_multiple"
nobs                                           ERROR no 'nobs' method is available
pp_check                                       ERROR no applicable method for 'pp_check' applied to an object of class "frmtmb_multiple"
loo                                            ERROR no applicable method for 'loo' applied to an object of class "frmtmb_multiple"
formula                                        ERROR invalid formula
update                                         ERROR need an object with call component
pp_check(mv, resp = 'y1')                      ERROR pp_check() is not supported yet for multivariate fits
frm(threshold = 'equidistant')                 ERROR unused argument (threshold = "equidistant")
summary(uv, waic = FALSE)                      ERROR summary() cannot honor `waic`: WAIC averages the pointwise likelihood over posterior draws, ...
:314 fitted(old_levels, ndraws = 10)           ERROR fitted() cannot honor `ndraws`: a maximum likelihood fit carries one estimate, not a posterior, ...
:955 update(fit2, bf(count ~ a + b, nl))       ERROR The nonlinear fit failed from its default starting values (NA/NaN gradient evaluation). ...
ce(draws, method = 'predict')                  ERROR conditional_effects() on draws has no method = "predict": the curves ARE posterior expected-response draws. ...
summary(draws, waic = TRUE)                    ERROR summary() has no argument `waic`. It takes: object, priors, prob, robust, mc_se
```

Three of these are silent wrong answers, not refusals: `coef()`,
`fitted()` and `residuals()` of a `frm_multiple()` result returned
`NULL`, because stats' default methods read `$coefficients` and
`$fitted.values`, which the object does not have. `summary(fit, waic =
FALSE)` was a refusal that fired where its condition (asking for WAIC)
was absent.

Lane remeasure's "Ranked gaps at 0.70.0"
(`frmtmb-wt-remeasure/dev/remeasure-findings.md`, read 2026-10-13
16:32) gives this lane ranks 6 (`pp_check(resp = )` on a multivariate
model), 7 (`threshold =`), 8 (`summary(waic = TRUE)`), 11
(`conditional_effects(method = "predict")` on draws) and 23
(`frm_sample()` of a `frm_multiple()` result). Rank 23 was not in the
brief; it is done here (item g).

## brms 2.23.0 on the same calls

`dev/postfix-brms.R` (log `dev/postfix-log/brms.txt`; `brm_multiple()`
on `mice::nhanes`, 2 imputations of 1 chain, 600 iterations, seed 1;
the mixed model adds `(1 | g)` with five arbitrary groups):

```text
brm(threshold = 'equidistant', chains = 0)     ERROR passing unknown arguments: threshold.
brm(cumulative(threshold =), chains = 0)       OK brmsfit [23]
stancode(threshold = ) has delta: FALSE
stancode(cumulative(threshold = )) has delta: TRUE
-- fm  class: brmsfit_multiple brmsfit
coef                                           ERROR No group-level effects detected. Call method 'fixef' to access population-level effects.
VarCorr                                        OK list [1]
vcov                                           OK matrix/array [4x4]
  [warning] Using only the first imputed data set. Please interpret the results with caution until a more principled approach has been implemented.
fitted                                         OK matrix/array [25x4]
fitted(newdata)                                OK matrix/array [3x4]
  [warning] Using only the first imputed data set. ...
predict                                        OK matrix/array [25x4]
  [warning] Using only the first imputed data set. ...
residuals                                      OK matrix/array [25x4]
emmeans                                        OK emmGrid [1]
  [warning] Using only the first imputed data set. ...
logLik                                         OK matrix/array [600x25]
nobs                                           OK integer [1]
  [warning] Using only the first imputed data set. ...
pp_check                                       OK ggplot2::ggplot/ggplot/ggplot2::gg/S7_object/gg [1]
  [warning] Using only the first imputed data set. ...
loo                                            OK psis_loo/importance_sampling_loo/loo [600x25]
summary(waic = TRUE) prints the same as summary(): TRUE
formals(summary.brmsfit): object, priors, prob, robust, mc_se, ...
pp_check(mv, resp = 'y1')                      OK ggplot2::ggplot/...
pp_check(mv)                                   ERROR Argument 'resp' must be a single variable name when applying this method to a multivariate model.
pp_check(mv, resp = c('y1', 'y2'))             ERROR Cannot coerce 'resp' to a single character value.
pp_check(mv, resp = 'nope')                    ERROR Invalid argument 'resp'. Valid response variables are: 'y1', 'y2'
pp_check(mv rescor, resp = 'y1')               OK ggplot2::ggplot/...
method = 'predict' is posterior_predict: TRUE
```

`conditional_effects(pf, method = "posterior_predict")` on a poisson
brmsfit has `estimate__` 1, `se__` 1.4826, `lower__` 0, `upper__` 3 at
the low end of the grid: the median, the MAD and the quantiles of the
predictive draws. brms's NEWS dates `threshold`: deprecated in `brm()`
at 1.8.0 ("Deprecate argument `threshold` in `brm` and instead
recommend passing `threshold` directly to the ordinal family
functions", `dev/brms-suite/brms/NEWS.md` line 1333) and removed with
the deprecated arguments at 2.0.0 (line 1172).

brms answers every in-sample method of a `brm_multiple()` fit from the
FIRST imputed data set and warns that it does ("until a more principled
approach has been implemented").

## a. Pooled post-processing of `frm_multiple()`

**Decided per method.** A method has a Rubin's-rules meaning when each
imputation gives an estimate and a standard error of the same quantity;
then the pooled estimate is their mean and the total variance the mean
within-imputation variance plus `1 + 1/m` times the between variance.

Pooled:

- `fixef()`, `summary()`, `conditional_effects()`: done in 0.68.0.
- `vcov()`: Rubin's total covariance of the `fixef()` rows, from each
  fit's `vcov()` and estimates; `mice::pool()`'s `t` is its diagonal.
- `fitted()`: row by row on the response scale, where one fit's Wald
  interval is symmetric, with a t interval on the Barnard-Rubin df. In
  sample each imputation predicts its own rows (brms uses the first
  data set and warns); with `newdata` all predict the same rows.
- `coef()`, `ranef()`: the mean of the levels' coefficients and modes;
  `ranef(condVar = TRUE)` gives Rubin's total variance, since a
  conditional mode with its conditional variance is an estimate with a
  variance.
- `VarCorr()`: on the log (sd) and Fisher-z (cor) scales, the scales
  `frm_multiple()`'s pooled table already uses, back-transformed, with
  `Est.Error` by the delta method; a covariance on its own scale.
- emmeans: `recover_data()` and `emm_basis()` as emmeans pools a
  `mice::mira`, so the result is emmeans's own on `glm()` fits.
- `nobs()`, `formula()`, `family()`: the value every imputation shares.

Refused by name, each naming what to call:

- `predict()`: a predictive interval is a quantile of simulated
  responses, not an estimate with a standard error; brms mixes the
  imputations' draws. Names `fitted()` and `predict(x$fits[[j]])`.
- `residuals()`, `pp_check()`: an imputed response differs between
  imputations; brms uses the first data set. Names `x$fits[[j]]`.
- `logLik()` (so `AIC()`, `BIC()`): m likelihoods of m data sets
  average to the likelihood of no model. Names `anova(m1, m0)`.
- `loo()`, `waic()`: posterior quantities. Names `anova()` and one
  imputation's draws.
- `update()`: each fit's stored call reads its data under a name local
  to `frm_multiple()`. Names `frm_multiple()` again.
- `model.frame()`: stats' default would rebuild a frame from
  `formula()` and its environment, which holds none of the data sets.

`frm_multiple()` itself now refuses imputations whose fits have
different coefficients and names them; 0.70.0 stopped inside
`vapply()` (`dev/postfix-repro-levels.R`, logs
`dev/postfix-log/repro-levels-{before,after}.txt`):

```text
before: "ERROR values must be length 4,\n but FUN(X[[2]]) result is length 5"
after:  "ERROR frm_multiple(): imputation 2 has other coefficients than the first (fc), so there is nothing to pool: a factor level, say, is missing from one of the imputed data sets"
```

**Validation** (`dev/postfix-validate.R`, log
`dev/postfix-log/validate.txt`). 1: poisson, n = 120, `x1` 25% missing,
data seed 11, `mice(m = 5, seed = 12)`, against `glm()` pooled by
`mice::pool()` (a poisson `glm()` has the maximum-likelihood covariance
`frm()` has, so the numbers must agree to the optimizers' precision).
2: gaussian on `nhanes` (mice seed 1) against `lm()` with its
covariance rescaled to the ML variance, `(n - p) / n`. 3: a
`(1 + x | g)` model with 20 groups, data seed 21, 5 hot-deck
imputations. Largest relative differences:

```text
## poisson, m = 5, n = 120 (seed 11, mice seed 12)
fixef Estimate vs pool()                                   1.351e-05
fixef Est.Error vs pool() std.error                        4.232e-06
pooled df vs pool() df                                     3.101e-05
vcov() vs glm() total covariance (all entries)             8.463e-06
diag(vcov()) vs pool() t                                   8.463e-06
fitted(newdata) Estimate vs Rubin over predict.glm()       1.833e-06
fitted(newdata) Est.Error vs Rubin over predict.glm()      7.396e-06
emmeans(~ x2) emmean vs emmeans(mira)                      4.402e-05
emmeans(~ x2) SE vs emmeans(mira)                          6.664e-06
emmeans(~ x2) df vs emmeans(mira)                          2.709e-04
pairs() estimate vs emmeans(mira)                          1.351e-05
pairs() SE vs emmeans(mira)                                7.004e-06
pairs() df vs emmeans(mira)                                4.281e-04

## gaussian, nhanes, m = 5 (mice seed 1)
fixef Estimate vs lm()                                     8.851e-07
vcov() vs lm() total covariance at the ML variance         2.150e-07
nobs()                                                     0.000e+00

## gaussian with (1 + x | g), 20 groups, hot-deck imputations, m = 5 (seed 21)
ranef() mode vs mean of the modes                          0.000e+00
ranef(condVar) condSD vs Rubin by hand                     0.000e+00
coef() = fixef() + ranef() (mean of sums)                  4.637e-15
VarCorr() sds vs the pooled_varcorr table (log scale)      0.000e+00
VarCorr() sd 2.5% vs pooled_varcorr lwr (same t quantile)  5.345e-11
VarCorr() cor vs the pooled_varcorr table (Fisher z)       0.000e+00
VarCorr() cor 97.5% vs pooled_varcorr upr                  2.416e-12
VarCorr() sd vs Rubin by hand on the log scale             0.000e+00
VarCorr() sd Est.Error vs the delta method by hand         2.196e-16
```

The zeros in block 3 are identities (the method computes the same
expression), stated as a check that the code computes what the
documentation says; the 5e-11 and 2e-12 are the difference between
`VarCorr()`'s finite-difference delta method and the theta covariance
the pooled table reads. The emmeans degrees of freedom agree to 4e-4
rather than 1e-5 because a df near 800 moves by 0.1 for a 1e-5 change
in the between variance. emmeans's pooled df of a `frm_multiple()` uses
the asymptotic complete-data reference that emmeans uses for one fit
(`Inf`), where `fixef()` uses `df.residual()` as `mice::pool()` does;
each matches its reference package.

## b. `pp_check(resp = )` on a multivariate fit

`pp_check.frmtmb_fit()` selects the response with `pp_check_resp()`,
whose three refusals are brms's words, and simulates that response's
marginal replicates with `sim_fit_draws(resp = )`, an internal
argument: the group effects are drawn for the whole model (a block two
responses share is one draw), a residual correlation is not drawn, and
the marginal of one response does not depend on it. A Student-t rescor
model's shared `nu` is added (`with_shared_nu()`), and a `subset()`
response is checked on its own rows, and on the newdata rows its
subset keeps. Public `simulate()` still refuses a multivariate fit.
frmtmb.sample's draws method speaks brms's words for the same three
refusals (it said "not supported yet" for a missing `resp`).

After (`dev/postfix-log/probe-after.txt`):

```text
pp_check(mv, resp = 'y1')                      OK ggplot2::ggplot/ggplot/ggplot2::gg/S7_object/gg [1]
pp_check(mv)                                   ERROR Argument 'resp' must be a single variable name when applying this method to a multivariate model.
pp_check(mv, resp = 'y2', type = 'stat')       OK ggplot2::ggplot/ggplot/ggplot2::gg/S7_object/gg [1]
pp_check(uv)                                   OK ggplot2::ggplot/ggplot/ggplot2::gg/S7_object/gg [1]
pp_check(uv, resp = 'y1')                      OK ggplot2::ggplot/ggplot/ggplot2::gg/S7_object/gg [1]
```

The test (`test-pp-check-types.R`) asserts that the observed series is
the named response's, in row order, and that the replicates' grand
mean sits within 0.1 of that response's own sd of its mean, while the
two responses are more than 2 sds apart, so a swapped response fails.

## c. `threshold =` as an `frm()` argument

A refusal, as the brief allowed: brms 2.23.0 fails the same line. `frm()`
has a trailing formal `threshold = NULL` that only refuses, naming the
family argument with the family of the call:

```text
frm(threshold = 'equidistant')                 ERROR frm() has no `threshold` argument: the threshold structure is an argument of the ordinal family, as in brms 2.23.0, whose brm() refuses it too ("passing unknown arguments: threshold"); brms 1.x took it in brm(). Write family = cumulative(threshold = "equidistant")
frm(family = cumulative(threshold =))          OK frmtmb_fit [21]
```

The absent case: R matches `thres = 1` to the formal by prefix, and the
first version refused it as `threshold` (`probe-after.txt` of the first
run). A prefix written in the call now gets R's own "unused argument
(thres = 1)", as before the formal existed. A prefix that reaches
`frm()` through another function's dots (`frm_multiple(thres = 1)`) is
not visible in `sys.call()` and is refused as `threshold`; that is
named in the code.

## d. `summary(waic = TRUE)`

brms 2.23.0 has no `waic` formal and prints `summary(fit, waic = TRUE)`
exactly as `summary(fit)` (identical `capture.output()`, above). frmtmb
keeps the refusal of `waic = TRUE` (NEWS 0.58.0's rule: a swallowed
argument that asked for something is refused, and old brms printed
WAIC there), now through an explicit `waic = FALSE` formal so that
`waic = FALSE`, which asks for nothing, is accepted. 0.70.0 refused
`waic = FALSE` too. frmtmb.sample's summary of draws refuses `waic =
TRUE` by name, naming `waic(x)`, where it gave the generic "has no
argument `waic`".

## e. `conditional_effects(method = "predict")` on draws

frmtmb.sample's method resolves the method with core's `ce_method()`
and, for `"predict"`, draws one response per grid row and posterior
draw with `ce_sim_one()` (the natural-scale dpars on the grid, the
addition terms held where the grid holds them, the family simulator
through `sim_response()`), and through `ce_plan_eval(one = )` when a
plan draws new levels per draw. Under the draws method's `robust =
TRUE` default `estimate__` is the median and `se__` the MAD, as brms's
are. `transform` applies to the predicted draws, as brms passes it to
`posterior_predict()`; a category-valued response is refused in the fit
method's words. The test compares the band with `posterior_predict()`
of the same draws on the same grid (median within 0.3 and upper bound
within 0.5 of the predictive sd; independent simulations), requires
`se__` above twice the expected-response band's, and checks that a
shift transform moves every quantile by the shift.

## f. The two bin-1 ledger defects

**`brmsfit-methods:314` (and `:317`).** brms's `ndraws` chooses how
many posterior draws brms summarizes: it sets the Monte Carlo error of
the answer, not the quantity answered. A maximum-likelihood fit answers
that quantity at its estimate with no Monte Carlo error, so on a fit
`ndraws` is accepted, changes nothing, and a message says so
(`fitted()` and `residuals()`). Drawing parameters from the asymptotic
normal (as `predict(ndraws = )` does for a predictive interval) was
considered and not done: it would make `fitted()`'s answer a Monte
Carlo estimate of what the delta method gives exactly, and its unseen
levels under `"gaussian"` would need per-draw effects that
`fitted_point()` does not take. `draw_ids`, `sort`, `summary = FALSE`
and `robust = TRUE` ask for particular draws and stay refused. After:

```text
:314 fitted(old_levels, ndraws = 10)           OK matrix/array [100x4]
:317 fitted(gaussian, ndraws = 1)              OK matrix/array [100x4]
```

**`brmsfit-methods:955`** (history: fix round 2 removed this scan,
and the row is a defect again; fix round 1 had already changed the
rule so that a scan restart stood only on a clean maximum). `update(fit2, formula. = bf(count ~ a + b,
nl = TRUE))` without fit2's priors stops at its zero start: `mu = a +
b = 0` on `Gamma("identity")`. The rule is `nl_intercept_scan()`
(`R/fit.R`), offered by `fit_recovery_starts()` after the first attempt
fails, so a fit that converged never meets it. It moves the intercepts
of the nonlinear parameters (every slope at zero) to one common value
over `+-10^(-2:3)`, then twice over each intercept alone over
`{0, +-1, 2, 5 x 10^(-2:3)}`, keeping each lower finite objective, and
restarts from the result with a message naming it. A coefficient that
`start` or a prior placed keeps its place.

The common value alone was the first version, and it was not enough.
On `brms_nonlinear`'s insurance-loss curve without priors the best
common value is 10 (nll 558.11, against 455.46 at the hand start
(5000, 1, 45)), and nlminb from there stopped 100.73 log-likelihood
units below the maximum, converged, with only a standard-error warning.
Coordinate passes over the powers of ten alone stopped 101.30 units
below, because the best power-of-ten point is (1000, 10, 10), nll
508.6. On the 1-2-5 grid the best joint point is (5000, 1, 50), nll
455.60 (`dev/postfix-nlscan-loss.R`, log
`dev/postfix-log/nlscan-loss.txt`), and the coordinate passes over that
grid reach the maximum. The three rules, swapped into the lane build
by `dev/postfix-nlscan-variants.R <lib> common|coord10|built` (logs
`dev/postfix-log/nlscan-{common,coord10,built}.txt`), on model 7 below:
`logLik(default) - logLik(hand)` = -1.007e+02, -1.013e+02 and
4.218e-11; on models 1 to 6 all three reach the hand start's
log-likelihood within 2.1e-9.

`dev/postfix-nlstart.R` on both builds (logs
`dev/postfix-log/nlstart-{before,after}.txt`): each model fitted from
the default start and from a start placed by hand near the truth.

```text
model                              before            after: logLik(default) - logLik(hand)
1 brmsfit-methods:955 update       ERROR (NA/NaN)    5.414e-11  (and the flat-direction warning naming a_Intercept, a_Age, b_Intercept, b_Age)
2 a + b, Gamma identity            ERROR (NA/NaN)    1.762e-11
3 b1 exp(b2 x), Gamma identity     ERROR (NA/NaN)    2.006e-09
4 logistic growth, gaussian        ERROR (NA/NaN)    2.022e-10
5 log(a) + b x, poisson            ERROR (NA/NaN)    1.478e-12
7 brms_nonlinear fit_loss          ERROR (NA/NaN)    4.218e-11
6 b1 exp(b2 x), gaussian (works)   conv 0, estimates 2.1382253772049342 0.47491431361041336 -0.70797714681122481 on both builds
```

Model 6 is the absent case: its default start fits, the scan never
runs, and its estimates are the same to every digit on both builds.
Model 1 is the ledger's: its body reads `a` and `b` only through their
sum with shared terms, so the maximum is a ridge, and the fit reaches
it and warns, as lane fixes recorded from a start placed by hand. A
start the caller gave is not scanned: `start = list(beta = c(0, 0,
0))` on model 4 still stops (`test-nl.R`).

Overlap: lane optstart changes the sigma and mu starts of
location-scale fits. This change is in `fit_recovery_starts()` and a
new `nl_intercept_scan()` in `R/fit.R`, reads `make_start()` but does
not change it, and runs only on a nonlinear model after a failed first
attempt.

## g. `frm_sample()` of a `frm_multiple()` result (remeasure rank 23)

`frm_sample(x)` with a `frm_multiple()` result re-evaluates the
caller's call once per imputation with that imputation's fit, seed
`seed + j - 1`, and combines the stanfits with `rstan::sflist2stanfit()`
and the draws by `rbind()`, as `brm_multiple()` combines its chains.
The draws object keeps the first imputation's fit, so in-sample methods
read the first imputed data set, as brms's do. The test checks that the
combined draws are bitwise the `rbind()` of three separate
`frm_sample()` calls with seeds 5, 6 and 7, that `nchains()` is 3, and
that an imputation with other parameters is refused by name.

## h. Lane remeasure's R2 to R4, and `loo()` of several models on draws

Added by the coordinator mid-lane (remeasure record R2 to R4, near
line 1050 of `frmtmb-wt-remeasure/dev/remeasure-findings.md`).

- **R2, `loo()` and `waic()` of draws swallowed `newdata`.** They
  passed their dots to `loo::loo.matrix()`, which drops what it does
  not know, so `loo(ds, newdata = nd)` returned the in-sample elpd with
  no message. `newdata` cannot be made to work here as brms does:
  `log_lik()` of draws composes the fitted objective's density over the
  assembled frame and refuses `newdata`, and evaluating draws on a frame
  rebuilt from newdata would rebuild a smooth's basis on it, which
  changes what a coefficient means without changing its name. So
  `loo()` and `waic()` refuse `newdata`, `re_formula` and the other
  arguments `log_lik()` refuses, by name, and refuse `moment_match =
  TRUE` and `reloo = TRUE` (brms refits; draws cannot). `pointwise` and
  `k_threshold` change nothing and are accepted; `r_eff`, `save_psis`
  and `cores` reach the loo package. Any other name is refused.
- **R3, `logLik()` of draws** is `log_lik()`, as brms's is.
- **R4, `bayes_R2()` of draws takes `newdata` and `re_formula`**, which
  it hands to `posterior_epred()`, reading the response off newdata's
  rows. The ported row also reads `model.frame(fit2)` of draws, which
  had no method and stopped in stats' default ("object 'count' not
  found"); `model.frame()` of draws is now the fit's.
- **`loo()` and `waic()` of several draws** return brms's comparison:
  `loos`, `diffs` (`loo::loo_compare()`) and brms's pairwise
  `ic_diffs__`, of class `frmtmb_loolist` (not brms's `loolist`: the
  class-name rule of 2026-09-17), printed as brms prints it. The test
  checks `ic_diffs__` against the definition on two different models,
  and that one model twice differs by exactly zero.

The signature of `loo()` and `waic()` of draws is now brms's, `(x, ...,
compare, resp, ndraws, model_names)`, so `ndraws` must be named; it is
marked BREAKING in frmtmb.sample's NEWS.

**Bin 2 of the ported suite** (lane remeasure's
`dev/brmsport-bin2/frmtmb.sample/test-brms-suite-bin2-methods.R`, run
from this worktree with remeasure's `dev/brmsport-bin2-run.R` on both
builds by `dev/postfix-bin2.sh`; records and diff
`dev/postfix-bin2-log/sample-{before,after}.tsv`,
`dev/postfix-bin2-log/sample-diff.txt`; nothing in the remeasure tree
was written):

```text
rows: 102 ; outcome changed: 14
assert brmsfit-methods:1015  verdict defect     held/ok FALSE -> TRUE
assert brmsfit-methods:1016  verdict defect     held/ok FALSE -> TRUE
assert brmsfit-methods:121   verdict defect     held/ok FALSE -> TRUE
assert brmsfit-methods:434   verdict defect     held/ok FALSE -> TRUE
assert brmsfit-methods:446   verdict defect     held/ok FALSE -> TRUE
assert brmsfit-methods:447   verdict defect     held/ok FALSE -> TRUE
assert brmsfit-methods:448   verdict defect     held/ok FALSE -> TRUE
assert brmsfit-methods:449   verdict defect     held/ok FALSE -> TRUE
assert brmsfit-methods:469   verdict hollow     held/ok TRUE -> FALSE
setup brmsfit-methods:1014   verdict            held/ok FALSE -> TRUE
setup brmsfit-methods:120    verdict            held/ok FALSE -> TRUE
setup brmsfit-methods:445    verdict            held/ok FALSE -> TRUE
setup brmsfit-methods:468    verdict            held/ok TRUE -> FALSE
setup brmsfit-methods:470    verdict            held/ok TRUE -> FALSE
```

Eight 2T rows now hold. `:469` held hollowly (on a NULL `newdata`,
which the swallowing hid) and now fails at its setup on the
`newdata` refusal, so its honest verdict is a divergence: refused by
name, as `log_lik()` refuses it. `:470` reads `:469`'s object. The bin-2
files of frmtmb (`methods`, `log_lik`, `posterior_predict`,
`distributions`, `priors`) change no row
(`dev/postfix-bin2-log/core-diff.txt`).

## The vignette port, before and after

`dev/postfix-bv.sh before|after` runs seven hand-translated vignettes
(`brms_missings`, `brms_multivariate`, `brms_distreg`,
`brms_nonlinear`, `brms_overview`, `brms_monotonic`,
`brms_phylogenetics`, 349 rows) and `dev/postfix-bv-diff.R` diffs them
(log `dev/postfix-log/bv-diff.txt`, the run after item h):

```text
rows: before 349  after 349
rows whose outcome changed: 10
  FALSE -> TRUE: 10
  TRUE -> FALSE: 0
brms_distreg | SAMPLE: conditional_effects(fit1, method = 'predict')
brms_missings | SAMPLE: conditional_effects(s_imp2, method = 'predict')
brms_missings | SAMPLE: fit_imp1 (the pooled object)
brms_multivariate | ML: pp_check(fit1, resp = 'back')
brms_multivariate | ML: pp_check(fit1, resp = 'tarsus')
brms_multivariate | SAMPLE: conditional_effects(..., method = 'predict')
brms_multivariate | SAMPLE: loo(fit1, fit2)
brms_nonlinear | SAMPLE: the vignette's exact call (method = 'predict')
brms_overview | SAMPLE: conditional_effects(fit1, method = 'predict')
brms_phylogenetics | SAMPLE: loo(model_pois, model_normal)
```

In remeasure's ranking these are rank 6 (2 of its 3 hand rows; the
third, `pp_check(s_imp2)` without `resp`, brms refuses too), rank 11
(all 5), rank 18 (both) and rank 23 (its one row).

Eight more rows fail on both builds with a new message: the three
`frm_multiple()` draws refusals (the message now names
`frm_sample(x)`), `pp_check(s_imp2)` without `resp` (brms's words now;
brms refuses it too), the two `threshold =` rows, and the two SAMPLE
`summary(waic = TRUE)` rows. The workaround script
(`dev/brms-vignettes/_workarounds.R`, diff `dev/postfix-wa-diff.R`, log
`dev/postfix-log/wa-diff.txt`): the brms spelling runs on 25 of 35
rows before and 27 after, "WA start: brms_nonlinear fit_loss, no
prior, no start" and "WA mv: pp_check(fit1, resp = 'tarsus') on the
ML fit".

## Not done, and why

- **`predict()` and `conditional_effects(method = "predict")` of a
  `frm_multiple()` result** stay refused. A mixture of each
  imputation's simulated replicates is brms's answer and the natural
  extension, but it is not Rubin's rules, which pool an estimate and a
  standard error; the brief asked for the methods with a Rubin's-rules
  meaning. The refusals name `fitted()` and one imputation's fit.
- **`update()` of a `frm_multiple()` result** is refused rather than
  implemented: each fit's stored call reads its data through a variable
  local to `frm_multiple()`'s `lapply()`, so `update()` of a fit cannot
  re-evaluate it, and `frm_multiple()` keeps neither its call nor the
  data list. Storing them is the change that would make it work.
- **brms's "first imputed data set" warning** on in-sample methods of
  pooled draws (item g) is not copied; the documentation says it.
- **The mechanical rows of ranks 7 and 8** still fail: both are
  refusals by name now, where brms fails rank 7 and silently ignores
  rank 8.

## Defects found and not fixed

- The `vapply()` failure of `frm_multiple()` on imputations with
  different factor levels was found here and is fixed (item a).
- `fitted(ndraws = )` stays a Wald summary. Drawing the parameters
  from their asymptotic normal, as `predict()` does, would give
  `summary = FALSE` an answer; not done (item f).
- My first `dev/postfix-run.sh` edit (adding `POSTFIX_RS`) was made
  while that driver ran the first full tier. The per-file outputs are
  written by the R processes and its summary step completed (371 of
  371), but the totals quoted below come from the final tier, run with
  the unchanged script, and are counted by `dev/postfix-summarize.R`
  from the per-file outputs.
- Pre-existing, filed in fix round 1 (item 6):
  `conditional_effects(<frm_multiple()>, int_conditions = )` stops
  with "matched by multiple actual arguments". Also, the review's Emax
  model M6 converges silently 114.85 units below its maximum.

## Ledger rows

Bin 1 (`dev/brmsport-verdicts.tsv`, not edited): the gated
`test-brms-suite-methods.R` of core reports two stale rows, which is
the expected form (`dev/postfix-tests/gated-after/frmtmb--test-brms-suite-methods.txt`):

```text
DETAIL FAIL [fitted has reasonable outputs] brms brmsfit-methods:314 now HOLDS but is recorded as 'defect' (...)
DETAIL FAIL [fitted has reasonable outputs] brms brmsfit-methods:317 now HOLDS but is recorded as 'divergence' (...)
```

`brmsfit-methods:955` stays a defect (starting values). Its setup
stops as on 0.70.0. After fix round 2 there is no scan, and the error
names `a` and `b` and the two ways to give starts.
The assertion `is(up, "brmsfit")` could not hold anyway, under the
class-name rule. frmtmb.sample's bin-1 file changes no row (61 pass, 0
fail). Bin 2: item h's table (8 rows of frmtmb.sample's methods file
now hold, `:469` no longer holds hollowly); frmtmb's five bin-2 files
change no row.

## Tests run

Every count below is from `dev/postfix-summarize.R` over the per-file
outputs (one R process per file, `dev/release/run-tests.R` with
`FRMTMB_RELLIB` = the lane library; every file's `lib:` line reads
`C:/Users/adf44/source/r/wt-postfix-lib`).

Full ungated tier, `NOT_CRAN=true`, final builds
(`dev/postfix-tests/final-after`, by package in
`dev/postfix-log/final-by-pkg.txt`):

```text
frmtmb files=229 pass=16546 fail=0 err=0 skip=179 warn=0
frmtmb.coupling files=11 pass=542 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=29 pass=1747 fail=0 err=0 skip=3 warn=0
frmtmb.latent files=10 pass=360 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=15 pass=430 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=51 pass=2645 fail=0 err=0 skip=4 warn=0
frmtmb.spline files=15 pass=590 fail=0 err=0 skip=1 warn=0
TOTAL final-after: files 371 of 371; pass=23409 fail=0 err=0 skip=208 warn=0
```

The 208 skips are gates: 179 `FRMTMB_BRMS_FIT_TESTS`, 15
`FRMTMB_SCALE_TESTS`, 13 `FRMTMB_DRMTMB_FIT_TESTS`, 1 `FRMTMB_FUZZ`.
Gated tier (`FRMTMB_BRMS_FIT_TESTS`,
`FRMTMB_DRMTMB_FIT_TESTS`, `FRMTMB_FUZZ`; the bcm, brms agreement,
fuzz and ported brms-suite files, `dev/postfix-jobs-gated.txt`, on the
final builds, as the ungated tier; no package file changed after
either ran):

```text
frmtmb files=31 pass=2893 fail=2 err=0 skip=0 warn=0
frmtmb.sample files=5 pass=201 fail=0 err=0 skip=0 warn=0
```

The two failures are the stale ledger rows above. The first full
ungated tier, before item h, read pass=23374 fail=5: the brms#734
test in `test-diagnostics-ux.R` pinned the failure the scan now
recovers (rewritten: it now pins the message and the maximum), and the
two `test-message-uniqueness.R` files found brms's resp refusal
written twice in each package (now one `pp_check_resp()`, exported to
frmtmb.sample, and `brms_validate_resp()` in core).

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`; a 2000 x 2000 `crossprod()` takes 0.16 s there against 4.88 s
on the reference BLAS), the files whose numbers this lane's changes
move (`dev/postfix-jobs-ob.txt`):

```text
TOTAL ob-after: files 9 of 9; pass=1030 fail=0 err=0 skip=2 warn=0
```

Seen to fail on rellib-r8 (the edited test files on the base build;
`dev/postfix-tests/{edited,r2,fix2}-before.log`):

```text
frmtmb test-multiple-methods.R pass=41 fail=2 err=6
frmtmb test-brms-postfit.R pass=76 fail=0 err=1
frmtmb test-sratio-thresholds.R pass=31 fail=4 err=0
frmtmb test-brms-shapes.R pass=72 fail=0 err=1
frmtmb test-arg-refusal.R pass=112 fail=0 err=1
frmtmb test-pp-check-types.R pass=198 fail=0 err=2
frmtmb test-nl.R pass=85 fail=0 err=1
frmtmb test-diagnostics-ux.R pass=127 fail=0 err=1
frmtmb.sample test-draws-methods.R pass=152 fail=0 err=3
frmtmb.sample test-ppcheck-loo.R pass=29 fail=0 err=1
frmtmb.sample test-conditional-effects-draws.R pass=63 fail=2 err=1
frmtmb.sample test-loo.R pass=100 fail=6 err=4 warn=1
frmtmb.sample test-laplace-draws.R pass=89 fail=1 err=1
```

## Fix round 1 (review `dev/reviews/2026-10-13-postfix.md`)

**1. The scan's restart stands only on a clean maximum.** The review
found the first rule turned the default start's error into a wrong fit
on 3 of its 6 new models: a Michaelis-Menten curve (M2, convergence 1,
79.9 below the maximum), an exponential at scale 2e5 (M4, convergence
1, 181.8 below) and a sinusoid (M7, convergence 0 at a straight line,
176.4 below). It also showed that the grid had been tuned on model 7
of this lane. `nl_scan_accept()` (`R/fit.R`) now judges the restart
once the fit is built. It stands only when the fit converged (code 0)
and the standard-error analysis (`se_fit_analysis()`) loses no
parameter. Otherwise the default start's own error stands, as on
0.70.0, and it ends with the scan's point and a candidate `start`.
When the restart itself errors, the error names the point too
(`dev/postfix-fix1-msgs.R`, log `dev/postfix-log/fix1-msgs.txt`). The
message that names an accepted scan start now comes after the check.
The non-convergence hint of `check_convergence()` says "restarted from
a scan" when a scan start was used, where it said "the fit began at
zero".

The 14 models that fail at the default start: the review's 8
(`dev/postfix-rev-nlscan.R` and `dev/postfix-rev-nlscan7.R`, seeds 101
to 109; logs `dev/postfix-log/rev-nlscan{,7}-fix1.txt`) and this
lane's 6 (`dev/postfix-nlstart.R`, log
`dev/postfix-log/nlstart-fix1.txt`). The gap is to the best of a start
placed by hand and a multistart of 2 to 20 points:

```text
model                                         outcome                          default - best
M1  poisson(identity) a exp(b x)              rescued, conv 0, all SEs         -2.8e-13
M1b M1 with normal(0, 5) on b, ub = 0         rescued, conv 0, all SEs         -2.6e-11
M1c M1 with normal(0, 100) on a, lb = 0       error; no finite scan point      -
M2  Gamma(identity) Michaelis-Menten          error, candidate start named     -
M4  Gamma(identity) a exp(-x / tau)           error, candidate start named     -
M7  gaussian sinusoid                         error, candidate start named     -
M8  Gamma(identity) a ~ 1 + (1 | g)           rescued, conv 0, all SEs         -3.6e-12
M9  Gamma(identity), six intercepts (ridge)   error, candidate start named     -
1   brmsfit-methods:955 update (ridge)        error, candidate start named     -
2   a + b, Gamma identity                     rescued                          1.8e-11
3   b1 exp(b2 x), Gamma identity              rescued                          2.0e-9
4   logistic growth, gaussian                 rescued                          2.0e-10
5   log(a) + b x, poisson                     rescued                          1.5e-12
7   brms_nonlinear fit_loss                   rescued                          4.2e-11
```

8 are rescued, each to its maximum within 2.1e-9; 6 keep the error. No
model returns a wrong fit without an error or a warning. The two
ridges (M9 and row :955) are refused because the maximum is not
identified; no start makes them clean. M2 and M7 are pinned in
`test-nl.R` as errors that name a candidate start, with the absent
case (the same model from a start in the right region converges).

M6, the review's Emax model with x = 0 in the data (seed 106,
`y ~ e0 + emax * x / (ed50 + x)`, gaussian, `emax = -8`,
`ed50 = 2`), is pre-existing and outside the scan: the default start
fits, converges (code 0) silently to logLik -155.23 against -40.38,
and the scan never runs. Filed below.

**2. The scan is computed only when every earlier restart has
failed.** `fit_recovery_starts()` now returns it as a function, last in
the list, and the restart loop calls it. On the review's case
(`dev/postfix-rev-flat.R`, a Gamma nonlinear model with `(1 | g)` and
a located prior, which "the placed intercepts with the nonlinear slopes
at zero" rescues), three fresh processes per build give the same
estimates to every digit on rellib-r8 and on the lane build, 0.50 s
on both (logs `dev/postfix-log/flat-{base,fix}-{1,2,3}.txt`):

```text
par: 4.8800347713965602 0.3415737626993442 2.1585475360700679 2.2853221669448156 -1.3023835890559512
proc 1 identical
proc 2 identical
proc 3 identical
```

**3. The diagnostics of combined imputations are brms's.** When the
draws combine more than one imputation (`draws_multiple()`, which reads
`nimp`), `summary()` leaves out `Rhat`, `Bulk_ESS` and `Tail_ESS`, as
brms's `summary.R:109` does. `loo()`, `psis()` and the loo types of
`pp_check()` use `r_eff = 1`, as brms's `loo.R:770` does. On the
review's case (`dev/postfix-rev-sampleC.R`, 3 imputations of 2 chains,
seed 7; log `dev/postfix-log/rev-sampleC-fix1.txt`) the summary has
the four columns, and `elpd_loo` is -129.53, the review's `r_eff = 1`
value (it was -130.66).

**4. `fitted()` of a `frm_multiple()` result pools on the link scale.**
Each imputation's estimate goes through the link of its predictor
(the logit for a category probability), its standard error through the
derivative, and the result comes back through the inverse link, as the
pooled `conditional_effects()` does. On the review's bernoulli case
(`dev/postfix-rev-rubin.R` section 2, seed 43, mice seed 44; log
`dev/postfix-log/rev-rubin-fix1.txt`), `fitted()` on the curve's grid
equals the curve: largest relative difference 2.9e-16 in the bounds,
0 in the estimate. At x = -1.906 both are 0.0833 (0.0040, 0.6741);
`fitted()` gave a lower bound of -0.1668 before. On the poisson case of
`dev/postfix-validate.R`, against `mice::pool.scalar()` of `glm()`'s
linear predictor carried back through `exp()`: Estimate 1.7e-6,
Est.Error 7.5e-6, Q97.5 4.2e-6, and every lower bound is positive.

**5. Minor.** The `threshold` refusal reads the family from
`bf() + cratio()`. Multi-model `loo()` attaches brms's `yhash`
(the observed responses; brms hashes them with digest, which
frmtmb.sample does not import), so `loo::loo_compare()` gives its own
warning, "Not all models have the same y variable", on models of
different responses, and stays silent on one response. `transform =` in
`conditional_effects(method = "predict")` and in `posterior_predict()`
on draws gives brms 2.23.0's warning, "Argument 'transform' is
deprecated and will be removed in the future." (brms
`posterior_predict.R:161`). The NEWS sentences are split to 25 words or
fewer, and the scan's bullet says "of either sign".

**6. Filed as pre-existing (not this lane's).**

- `conditional_effects(frm_multiple_result, "x", int_conditions = )`
  stops with "formal argument "int_conditions" matched by multiple
  actual arguments": `conditional_effects.frmtmb_multiple()` (0.68.0)
  passes its own `int_conditions` to every imputation after the first.
- The review's M6 (above) converges silently 114.85 log-likelihood
  units below the maximum from its default start. Lane optstart may
  reach it.

**Tests run after fix round 1.** The same driver and summary script as
"Tests run" above. These numbers replace the ones there. Full ungated
tier, `NOT_CRAN=true` (`dev/postfix-tests/fix1-final`, by package in
`dev/postfix-log/fix1-by-pkg.txt`):

```text
frmtmb files=229 pass=16561 fail=0 err=0 skip=179 warn=0
frmtmb.coupling files=11 pass=542 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=29 pass=1747 fail=0 err=0 skip=3 warn=0
frmtmb.latent files=10 pass=360 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=15 pass=430 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=51 pass=2652 fail=0 err=0 skip=4 warn=0
frmtmb.spline files=15 pass=590 fail=0 err=0 skip=1 warn=0
TOTAL fix1-final: files 371 of 371; pass=23431 fail=0 err=0 skip=208 warn=0
```

The first full tier of the round (`fix1-all`) had 5 failures and 1
error, all fixed before `fix1-final`:

- The conditions census flagged a bare `stop(e)` in
  `nl_scan_accept()`. `fit_error_context()` now takes `err =`.
- `test-multiple-methods.R` had a names mismatch, and the identity
  link dropped a dimension in `multiple_pool_summaries()`.
- `test-nl.R` used the old wording.
- `test-message-uniqueness.R` of frmtmb.sample found the transform
  warning written twice. It is now `warn_transform_deprecated()`.

The other runs:

```text
TOTAL fix1-gated: files 36 of 36; pass=3094 fail=2 err=0 skip=0 warn=0   (the stale rows :314 and :317)
TOTAL fix1b:      files 20 of 20; pass=1441 fail=0 err=0 skip=2 warn=0   (the affected, nonlinear, loo and multiple files)
TOTAL ob-fix1:    files 9 of 9;   pass=1051 fail=0 err=0 skip=2 warn=0   (OpenBLAS 0.3.26 with LAPACK)
```

On rellib-r8 the files edited in this round fail
(`dev/postfix-tests/fix1-before`, `dev/postfix-jobs-fix1c.txt`):

```text
frmtmb test-nl.R pass=89 fail=2 err=1
frmtmb test-multiple-methods.R pass=42 fail=5 err=7
frmtmb test-sratio-thresholds.R pass=31 fail=5 err=0
frmtmb test-diagnostics-ux.R pass=127 fail=0 err=1
frmtmb.sample test-draws-methods.R pass=152 fail=0 err=3
frmtmb.sample test-loo.R pass=100 fail=6 err=5 warn=1
frmtmb.sample test-conditional-effects-draws.R pass=63 fail=2 err=1
TOTAL fix1-before: files 7 of 7; pass=604 fail=20 err=18 skip=2 warn=1
```

Each of the new tests of this round is among the failures:

- "a scan restart that is no clean maximum keeps the error" (2 FAIL).
- "fitted() and the pooled curve agree, inside the support".
- "frm(threshold = ) reads the family of bf() + family()".
- "loo() of models on different responses warns as brms does".
- The `nimp` checks in "frm_sample() of a frm_multiple() result
  combines the imputations", which errors at `frm_sample()`.

The bin-2 methods file of frmtmb.sample on the final build
(`dev/postfix-bin2-log/sample-diff.txt`) gives item h's table again:
14 of 102 outcomes changed, the same 14 rows.

## Fix round 2 (review `dev/reviews/2026-10-13-postfix.md`, "Round 2")

**1. The scan is removed, and no fit comes from it.** The review's
stress block (`dev/postfix-rev2-nlnew.R`: `y ~ A * sin(2 pi x / P)`,
P in {1.3, 3, 4, 7, 9}, seeds 1001 to 1006, 30 data sets) gave on the
round-1 build: error 4, rescued 4, and 22 wrong fits with code 0,
finite standard errors and no warning, 93 to 118 log-likelihood units
low. A local check cannot tell those maxima from the global one. So
`nl_intercept_scan()`, `nl_scan_accept()`, `nl_scan_values()`,
`nl_scan_note()`, the lazy entry of `fit_recovery_starts()`, the
`err =` argument of `fit_error_context()` and the scan hint of
`check_convergence()` are gone. `R/fit.R` is 0.70.0's except for the
`threshold` refusal (item c) and the message below. The round-1 file
is kept as `dev/postfix-fit-round1.R.txt` and
`dev/postfix-fit-round1.diff`.

**2. The default-start error says how to give starts and offers none.**
`nl_start_message()` (`R/fit.R`) replaces the old text, whose example
was `start = list(beta = c(ult = 5000))`, a parameter of another model.
The new text names this model's nonlinear parameters. It gives two
ways: `start =` with the names from `par_template()`, or a located
prior, shown as `prior(normal(<location>, <scale>), nlpar = "A")`. It
says that frmtmb offers no start of its own, because a searched start
can reach a local maximum with no warning. It names no point. frmtmb's
`restarts` control restarts from the optimum, not from new points, so
the message does not name it. The sinusoid of the stress block (P 3,
seed 1001; `dev/postfix-fix2-msg.R`, log
`dev/postfix-log/fix2-msg.txt`) gives:

```text
The nonlinear fit failed from its default starting values (NA/NaN gradient evaluation). A nonlinear coefficient starts at 0 unless `start` or a located prior places it, and the body is often not defined there. frmtmb offers no start of its own: a start found by search can lead to a local maximum, which converges with no warning. Give the starts of `A`, `P` in one of two ways. (1) `start = list(beta = c(<name> = <value>))`: par_template(formula, data) gives the names, and the names go straight back. (2) A normal() or student_t() prior with a location, e.g. prior(normal(<location>, <scale>), nlpar = "A"): the fit starts from the location
```

`test-nl.R` pins it on that case, with the absent case (the same model
from `start = c(2, 3)` converges). On rellib-r8 the test fails 4 of
its expectations (`dev/postfix-tests/fix2-nl-before`).

**3. What reverts.** Every model the scan rescued stops again, now with
the new message. Logs are
`dev/postfix-log/{nlstart,rev-nlscan,rev-nlscan7}-fix2.txt` and
`dev/postfix-log/fix2-nlnew.txt`
(`dev/postfix-fix2-nlnew.R` is the review's script with its candidate
step skipped on the message's `<name> = <value>` example, which it
otherwise parses and fails on):

```text
dev/postfix-nlstart.R        6 models error (:955, a + b, b1 exp(b2 x), logistic, log(a) + b x, fit_loss); model 6 fits as on 0.70.0
dev/postfix-rev-nlscan.R     M1, M1b, M1c, M2, M4, M8, M9 error; M6 fits silently 114.85 low (pre-existing); M7 stops on `pi` (the script's formula)
dev/postfix-rev-nlscan7.R    M7 errors
dev/postfix-fix2-nlnew.R     N1 to N6 error; stress block: error 30 of 30
```

- Ledger `brmsfit-methods:955` is a defect again (starting values).
  It was never flipped in the final build of round 1, so the
  flipped-rows list does not change.
- The workaround row "WA start: brms_nonlinear fit_loss, no prior, no
  start" reverts: the brms spelling runs on 26 of 35 rows, against 25
  on rellib-r8 and 27 on round 1 (`dev/postfix-log/wa-diff-fix2.txt`).
- The vignette diff is unchanged: 10 rows FALSE to TRUE, 0 TRUE to
  FALSE (`dev/postfix-log/bv-diff-fix2.txt`). None of the 10 came from
  the scan.
- `dev/postfix-rev-flat.R` gives rellib-r8's estimates to every digit
  (`dev/postfix-log/flat-fix2.txt`).
- `test-diagnostics-ux.R` is 0.70.0's file again (the brms#734 test
  expects the error).
- NEWS: the two scan bullets are removed. A bug-fix bullet describes
  the new message.

**4. Pooled category probabilities.** `multiple_pool_summaries()`
(`fitted()`) and `conditional_effects.frmtmb_multiple()` now pool a
category probability's estimate and standard error by Rubin's rules on
the probability scale. Each imputation's categories sum to 1, so their
means do too. The interval ends still come from the logit scale. An end
on the wrong side of the estimate is moved to the estimate; on the
cases below that never happens. Other links (poisson, bernoulli,
gaussian) are as in round 1. `dev/postfix-catpool.R`, log
`dev/postfix-log/catpool.txt`:

```text
cumulative, 3 imputations (the review's case), x = -3, 0, 3 and in sample (60 rows)
  (1) |sum - 1| at most 0 and 1.1e-16; conditional_effects() 2.2e-16
  (2) estimates outside their interval: 0 and 0; conditional_effects() 0; ends moved: 0
  (3) fitted() vs conditional_effects() on the curve's grid: est 0, lower 1.8e-11, upper 1.7e-11
  (4) bounds in [7.2e-05, 0.9965] and [7.7e-04, 0.9702]
  mean of probabilities vs inverse logit of the mean logit: largest 0.0031 and 0.0046
categorical, 3 categories, 5 imputations, x = -3, 0, 3 and in sample (150 rows)
  (1) |sum - 1| at most 1.1e-16
  (2) estimates outside their interval: 0; ends moved: 0
  (4) bounds in [1.6e-04, 0.9922] and [8.9e-05, 0.9927]
  mean of probabilities vs inverse logit of the mean logit: largest 3.9e-04 and 0.052
```

For a categorical fit, property (3) has no curve to compare.
`conditional_effects()` of one categorical fit has no analytic band
("Use band = \"boot\""), and the pooled method refuses the bootstrap
band. That is older than this lane.

The review's `dev/postfix-rev2-rubin.R` on the new build
(`dev/postfix-log/rev2-rubin-fix2.txt`) differs from its round-1 log
in one line only: "probs sum 0.9966863 0.9998307" is now "probs sum 1
1". The poisson, bernoulli and gaussian lines are identical.

The ordinal test of `test-multiple-methods.R` (0.68.0) asserted the
old rule (the inverse logit of the mean logit). It now asserts the
mean, the sum of 1, `fitted()` equal to the curve, estimates inside
their intervals, and bounds in [0, 1]. On rellib-r8 it fails 5 and
errors 1 (`dev/postfix-tests/fix2-mm-before`).

**Tests run after fix round 2.** R/fit.R changed, so the full tiers ran
again, not only the files that reach `fit_recovery_starts()`. Full
ungated tier, `NOT_CRAN=true` (`dev/postfix-tests/fix2-final`):

```text
frmtmb files=229 pass=16558 fail=2 err=0 skip=179 warn=0
frmtmb.coupling files=11 pass=542 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=29 pass=1747 fail=0 err=0 skip=3 warn=0
frmtmb.latent files=10 pass=360 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=15 pass=430 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=51 pass=2652 fail=0 err=0 skip=4 warn=0
frmtmb.spline files=15 pass=590 fail=0 err=0 skip=1 warn=0
TOTAL fix2-final: files 371 of 371; pass=23428 fail=2 err=0 skip=208 warn=0
```

The two failures were in the new tests themselves: `tapply()` keeps a
dimension in `test-multiple-methods.R`, and "ult" matched "default" in
`test-nl.R`. Both are corrected, and the two files were run again:

```text
frmtmb test-multiple-methods.R pass=113 fail=0 err=0   (dev/postfix-tests/fix2-after)
frmtmb test-nl.R               pass=92  fail=0 err=0   (dev/postfix-tests/fix2-nl)
```

With these two, frmtmb has 229 files, 16560 pass, 0 fail, 0 error and
179 skips. The other runs:

```text
TOTAL fix2-gated: files 36 of 36; pass=3094 fail=2 err=0 skip=0 warn=0   (the stale rows :314 and :317)
TOTAL fix2-after: files 6 of 6;   pass=675 fail=1 err=0 skip=2 warn=0    (test-nl.R's "ult", above)
TOTAL ob-fix2:    files 6 of 6;   pass=676 fail=0 err=0 skip=2 warn=0    (OpenBLAS 0.3.26 with LAPACK, after the fix)
```

`fix2-after` and `ob-fix2` are `test-nl.R`, `test-multiple-methods.R`,
`test-multiple-pooling.R` and `test-diagnostics-ux.R` of core, and
`test-draws-methods.R` and `test-loo.R` of frmtmb.sample
(`dev/postfix-jobs-fix2.txt`). The OpenBLAS copy of R was rebuilt by
`dev/ciharden-openblas.sh 0.3.26 lapack` (log
`dev/postfix-log/openblas-build2.txt`; a 2000 x 2000 `crossprod()`
takes 0.03 s on it).

## Status

FINAL. Every item of the brief is done or decided, with remeasure's
ranks 6, 7, 8, 11 and 23 and the coordinator's R2 to R4 and
multi-model `loo()` (item h).

`R CMD check --as-cran` on the final build of fix round 2, built and
checked in `dev/postfix-check/`
(`dev/postfix-check.sh`, `R_LIBS` = lane library, rellib-r8, user
library; logs kept as `dev/postfix-log/check-{core,sample}-full.txt`,
the check trees deleted):

```text
frmtmb:        Status: 1 NOTE   (checking HTML version of manual: "Skipping checking math rendering: package 'V8' unavailable", the expected one)
frmtmb.sample: Status: OK
```

Bump: minor for both packages (new methods and behavior changes,
BREAKING bullets in both NEWS). Floor: frmtmb.sample needs the frmtmb
release that carries this lane (it calls the new exports `ce_sim_one()`,
`ce_method()`, `ce_predict_family_check()`, `pp_check_resp()` and
`ce_plan_eval(one = )`).

Overlaps with other lanes: `R/fit.R` (optstart: the start path;
after fix round 2 this lane adds only `nl_start_message()`, called from
`fit_error_context()`, and the `threshold` formal of `frm()` with
`frm_refuse_threshold()`), `R/multiple.R`,
`R/predict.R` (`fitted_refuse_draws_args()`, `sim_fit_draws(resp =
)`), `R/conditional-effects.R`, `R/methods-fit.R` (`summary()`'s
`waic` formal), `NEWS.md`, frmtmb.sample (`R/loo.R`, `R/sample.R`,
`R/methods-draws.R`, `R/conditional-effects-draws.R`).
