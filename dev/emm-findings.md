# Lane emm: emmeans on nonlinear and multivariate fits

Lane `emm`, 2026-09-25, worktree `/home/user/wt/emm`, built into
`/opt/rlib/lane-emm`. Base for every "before" below is `/opt/rlib/base`.

## What changed

`R/interop.R` replaces `emm_mu_linpred()` and the two emmeans methods
with brms's interface. `recover_data()` and `emm_basis()` take `resp`,
`dpar`, `nlpar`, `epred` and `re_formula` with brms's defaults
(`re_formula = NA`, `epred = FALSE`), and brms's deprecated
`dpar = "mean"` with brms's warning. One resolver, `emm_target()`,
serves both methods, so the grid and the basis always describe the same
predictor. It picks one of two routes:

- **design**: every selected predictor is linear with parametric terms
  only, at the population level. The basis is emmeans's usual one
  (design, coefficients, `vcov_estimated()` block). This is the old
  path, extended to a distributional parameter, a nonlinear parameter
  and to several responses (block diagonal). For the old default call
  (univariate `mu`) it builds the same `X`, `bhat` and `V` as before.
- **grid**: a nonlinear body, `epred = TRUE`, `re_formula` other than
  `NA`, or a predictor with `s()`, `t2()`, `gp()`, `mo()` or `mi()`.
  The basis is brms's shape: `X` is the identity over the grid rows,
  `bhat` the predictions, `V = G Sigma G'` with `Sigma` the joint
  covariance (`get_joint_cov()`) and `G` the Jacobian from
  `frm_lp_basis()`. That function already tapes a nonlinear body, so no
  new differentiator was written. For `epred = TRUE`, `G` is the chain
  rule `predict_mean_se()` uses (`mean_eta_grad()` per dpar times that
  dpar's `frm_lp_basis()` Jacobian), kept as a whole matrix instead of
  its diagonal, and it reaches a nonlinear dpar.

Without `resp`, a multivariate fit stacks its responses in the order of
`names(fit$spec$responses)`, grid rows fastest, and sets
`misc$ylevs = list(rep.meas = <response names>)`. That is exactly
brms's `emm_basis.brmsfit()` layout (it flattens the draws x grid x
resp array column-major and sets the same `misc$ylevs`).

`misc` now carries the link, built by `emmeans::.std.link.labels()` as
brms builds it, so `type = "response"` applies the inverse link once.
A link that `make.link()` does not know (`softplus`, `logm1`, ...) is
passed as the link object itself. `epred = TRUE` carries no link: the
estimates are already on the response scale. The ordinal location
keeps its latent convention and carries no link either.

### How the refusal message was lost, and the fix

emmeans 2.0.4 `ref_grid()` does

    data = try(recover_data(object, data = NULL, ...))
    if (inherits(data, "try-error"))
        stop("Perhaps a 'data' or 'params' argument is needed")
    ...
    if (is.character(data)) stop(data)

so any condition raised in `recover_data()` becomes the "Perhaps"
error (the real text goes to stderr only, through `try()`'s print).
The last line is emmeans's own convention: a `recover_data()` method
refuses by RETURNING a string. `recover_data.frmtmb_fit()` now wraps
its whole body in `tryCatch(error = function(e) conditionMessage(e))`,
so every error there, refusal or not, reaches the user as the error
text. `emm_basis()` is not wrapped in `try()` by emmeans, so its errors
already surfaced. Because all argument checks run in
`emm_target()`, which `recover_data()` calls first, every refusal about
the arguments is raised there and is visible.

## Silent wrong answers found on the way, all fixed

Each was reproduced on base, and each has a test in
`tests/testthat/test-emmeans.R` that fails on base (see Tests).

1. **`dpar` was ignored.** `emmeans(fit, "f", dpar = "sigma")`
   returned mu's marginal means under the requested label. brms port
   `emmeans:18` passed only because it counts rows. Now it is sigma's
   predictor on the log scale.
2. **No link reached emmeans.** `misc` was empty, so on a Poisson fit
   `type = "response"` printed the log-scale means unchanged, with no
   note. Now `rate` is `exp()` of the mean.
3. **Smooths were left out of the means.** The design basis used the
   parametric terms only, so `y ~ f + s(x)` reported means without the
   smooth's value (5.28 against `frm_linpred()` 5.33 on the probe
   model). Such predictors now take the grid route.
4. **`mo()` terms were left out of the means.** On a probe fit of
   `y ~ f + mo(e)` base printed 4.05 for level A where the prediction
   at the grid's mean `e` is 5.01; the lane prints 5.01. (The NaN
   standard errors on that probe come from the fit's own covariance,
   which is NaN for the intercept and `moe`, not from emmeans.) The
   `mo()` variable now joins the grid; brms ports `emmeans:21` and
   `:24` pass as a result.
5. **An addition term was read as a response transformation.** The
   fit's call was handed to emmeans, which read `y | weights(w)` as a
   transformation named `"|.weights"`, labeled the scale with it, and
   made `type = "response"` a no-op with a note. The call is now
   `<response> ~ 1` for `mu` and `epred`, and `call("frm")` otherwise,
   which is brms's rule (`call("brms")` unless `epred`). `log(y) ~ f`
   still gets `tran = "log"` (tested).
6. **`re_formula` was ignored.** It now keeps the group-level effects;
   their grouping factors join the grid, restricted to the terms a
   formula keeps.

## Validation

All numbers are printed by the scripts named; nothing below is typed
from memory. Relative differences are `max(|a - b| / |b|)`.

### Nonlinear and multivariate against linear references

`dev/emm-validate.R`, log `dev/emm-validate-log.txt` (seed 1, n = 200).
The two fits of each pair are separate optimizations, so the first
line is the floor every other line sits on.

| comparison | rel. difference |
|---|---|
| coefficients, `a + b * x` nl fit vs `y ~ f + x` (floor) | 7.453e-07 |
| `nlpar = "a"` vs linear at `x = 0`: estimate / SE | 3.001e-07 / 3.188e-07 |
| whole mu (grid route) vs linear: estimate / SE | 1.925e-07 / 3.188e-07 |
| `pairs()` on whole mu vs linear: estimate / SE | 1.287e-06 / 3.188e-07 |
| `epred = TRUE` vs whole mu, identity link (identity) | 0 / 6.828e-12 |
| `resp = "y1"` vs univariate `y1 ~ f`: estimate / SE | 1.106e-06 / 5.727e-07 |
| stacked `rep.meas = y2` vs univariate: estimate / SE | 9.520e-06 / 2.392e-07 |
| `y1 - y2` contrast SE vs `sqrt(se1^2 + se2^2)` | 2.315e-07 |
| mean over `rep.meas` vs `(y1 + y2) / 2`: SE | 2.315e-07 |

The `epred` row is an identity (the mean of an identity-link gaussian
is its predictor), so its residual is a numerical check of the central
difference in `mean_eta_grad()`. The contrast rows check that the
stacked covariance has no spurious cross-response block when the
responses are independent. With `rescor = TRUE` the cross-response
covariance of the coefficients enters instead: the `y1 - y2` SE was
0.139 against 0.188 for the independent formula on a probe fit.

### Delta method against a numerical Jacobian and marginaleffects

Same log, seed 11, n = 200. `J` is `numDeriv::jacobian()` of
`frm_linpred()` at the grid under `set_coef()`, and the reference
covariance is `J V J'` with `V = interop_vcov()`.
`marginaleffects::predictions()` shares no code with the emmeans path
(it perturbs `get_coef()` and calls `get_predict()`), so it is an
independent numerical Jacobian.

| model and scale | bhat vs predict | V vs numDeriv | SE vs marginaleffects |
|---|---|---|---|
| `y ~ a * exp(-b * x)`, gaussian, link | 0 | 1.261e-09 | 2.403e-08 |
| `cnt ~ a - b * x`, poisson, link | 0 | 8.867e-08 | 1.502e-08 |
| same, `epred = TRUE`, response | 0 | 1.316e-07 | 2.986e-08 |

`regrid()` of the poisson link-scale grid equals the `epred` grid to
0 (estimate) and 3.014e-11 (SE): `type = "response"` applies `exp()`
once, and never on top of `epred`.

### Delta method against a parametric bootstrap

Same script, section 4. Refits with `update()` on responses from
`simulate()`; the quantities are the three marginal means over `f` and
their three `pairs()`. Ratio is delta-method SE over bootstrap SD.

| B | model | ratios | SD relative error |
|---|---|---|---|
| 1500 | exp decay, whole mu | 0.990 0.994 1.029 0.981 1.007 1.007 | 0.018 |
| 1500 | poisson nl, `epred` | 1.000 0.985 1.011 0.974 1.010 1.004 | 0.018 |
| 400 | exp decay, whole mu | 0.913 1.014 1.076 0.973 1.005 1.061 | 0.035 |
| 400 | poisson nl, `epred` | 1.010 0.994 1.033 0.991 1.041 1.009 | 0.035 |

Logs: `dev/emm-validate-log-b1500.txt` (B = 1500) and
`dev/emm-validate-log.txt` (B = 400); 1500 of 1500 and 400 of 400
refits succeeded. A dissolved signal, recorded as the rules ask: at
B = 400 the first ratio read 0.913, 2.5 of its own standard errors
below 1. At B = 1500, which extends the same simulation seed, it reads
0.990. The 400-draw value was an unlucky prefix, not a bias.

### brms's interface

`dev/emm-brms-compare.R`, log `dev/emm-brms-compare-log.txt`. brms's
`recover_data()` needs no draws, so it was run on `empty = TRUE` fits
beside frmtmb's on the same data. The grid variables, the refusal
texts and the `rep.meas` levels (`y1 y2`) agree, except:

- nonlinear `mu` with no `nlpar`: brms's grid holds only the body's
  covariates (`x`); frmtmb's holds `x f`. brms's `.extract_par_terms()`
  takes `allvars` of the `btnl` object, which lists the body's
  covariates only. The nonlinear parameters' own covariates are then
  not in brms's grid, so a call such as `emmeans(fit, "f")` fails there
  and averages over `f` here. Documented under "Divergence from brms"
  in `?frmtmb-emmeans`.
- mixed families without `resp`: brms refuses when the responses'
  `misc` lists differ, which includes the label (`rate` against
  `response`). frmtmb refuses when the LINKS differ. So gaussian and
  student (both identity) stack in both; a log-link gaussian with a
  poisson stacks here and not in brms, where the stacked scale is the
  same log scale for both.

Part 2 of that script would sample both models with brms and compare
the numbers. It is gated off: on this machine rstan cannot compile any
model: `StanHeaders` includes `tbb/tbb_stddef.h`, which the conda
`RcppParallel` build does not ship, and `compileCode()` stops with
"fatal error: tbb/tbb_stddef.h: No such file or directory". A brms
numeric comparison is therefore not in this
lane's evidence. The two estimators differ anyway (posterior mean and
SD against ML estimate and Wald SE), so it would compare meaning and
layout, which the grid comparison and the linear-reference identities
above already pin.

## brms suite ports

`dev/brmsport-verdicts.tsv`: `emmeans:21`, `:24`, `:27`, `:35`, `:38`
and `:50` flip to `pass`. `emmeans:42` stays `cannot transfer` with a
new reason: fit6's `volume` response has an exact `gp(Age)`, predicted
at the grid's mean `Age`, which is not a fitted position, and the
kriging covariance between grid points is not available; frmtmb
refuses by name.

The generator's source tarball `dev/brms-suite/brms_2.23.0.tar.gz` is
gitignored and absent here, CRAN is refused by the proxy, and GitHub's
`paul-buerkner/brms` has no `v2.23.0` tag. The `cran/brms` mirror has
tag `2.23.0`. I ran `dev/brmsport-gen.R` from a scratch copy whose
`dev/brmsport-blocks.R` points at that clone and skips the sha256 check
(`/tmp/lanes/emm/emm-brmsport-blocks.R`). Control: with the verdicts
unchanged, it rewrote all generated files with no diff at all. After
the verdict change the only file that changed is
`tests/testthat/test-brms-suite-emmeans.R`, in exactly the seven rows.

`test-brms-suite-emmeans.R` (gated): lane `pass=11 fail=0`; base
`pass=5 fail=6`, the six flipped rows.

## Refusals added

Each names its reason; texts follow brms where brms has one.

- `'dpar' and 'nlpar' cannot be specified at the same time.` (brms)
- `Non-linear parameter 'zz' is not part of the model...` (brms)
- `Distributional parameter 'nu' is not part of the model...` (brms)
- `Invalid argument 'resp'. Valid response variables are: ...` (brms)
- different links without `resp`: names the links and both ways out.
- no `mu` (categorical, mixture): names the dpars to choose from.
- `epred = TRUE` on an ordinal or categorical family, on an HMM-style
  structured family, and on a family with no mean.
- a grid row the rank-deficient fit cannot estimate, or a non-finite
  prediction (brms: "created NAs").
- an exact `gp()` off its fitted positions (see "Not done").

## Downstream methods

Nothing outside emmeans changed. `predict()`, `fitted()`,
`simulate()`, `summary()`, `confint()` and `frm_sample()` do not reach
this code. The marginaleffects methods do not share it; they were run
as the independent reference above and agree. The objective is not
touched, so fit time is unchanged by construction and was not measured.

## Registry (`R/compat.R`)

`emmeans x nl`: refused -> works. `rescor x emmeans`: refused -> works.
New `mvbf x emmeans` works (it was caught by `mvbf x kind:method`
refused). New rows: `s()` and `t2()` works, `mo()` works, `gp_pred()`
conditional, `mi_pred()` untested. `emmeans x kind:family` note
rewritten.

## Not done, and why

- **Exact `gp()` kriging covariance.** `pred_design()` returns the
  kriging variance as a per-row vector (`extra_var`); a marginal mean
  needs the covariance between grid rows. Building it means carrying an
  n x n matrix out of `pred_design()` for every `predict()` caller, or
  an opt-in path through `frm_lp_basis()`'s exported return value.
  Refused by name instead, with two ways out that work today
  (`gp(..., k = )`, or `at =` an observed position).
- **`epred = TRUE` on ordinal and categorical families.** brms returns
  the category probabilities as a `rep.meas`-like factor. No delta
  method for category probabilities exists in the package
  (`frm_linpred(se.fit = TRUE)` refuses them too), so this is refused.
- **brms's `type = "response"` on an ordinal location.** brms passes
  the ordinal family's link to emmeans and so back-transforms a latent
  mean with no threshold, which is not a probability of anything.
  frmtmb keeps the latent convention and applies nothing.
- **se.fit on a nonlinear predictor in `predict()`/`frm_linpred()`**
  stays refused; the emmeans path reaches `frm_lp_basis()` directly.

## Tests

New `tests/testthat/test-emmeans.R`, 10 blocks, 38 expectations, about
9 s. Lane: `pass=38 fail=0 err=0 skip=0`. Base: `pass=4 fail=10 err=4`,
every block failing: behavioral failures (not missing symbols) for
`dpar`, `type = "response"`, the weights transform, the smooth and
`mo()`, `re_formula` and the `gp()` refusal; errors for the nonlinear,
multivariate and refusal-text blocks, where base raises "Perhaps a
'data' or 'params' argument is needed". Tolerances are ratios to a
floor the test measures (the two optimizers' coefficient difference,
or `sqrt(.Machine$double.eps)` for an identity) and, for the Jacobian,
to the error of a forward difference against Richardson's.

## NEWS entry

- `emmeans()` now takes brms's `dpar`, `nlpar`, `resp`, `epred` and
  `re_formula`. A nonlinear model is supported: `nlpar =` averages
  that parameter's linear predictor, and the whole `mu` or
  `epred = TRUE` uses the delta method through `frm_lp_basis()`. A
  multivariate fit without `resp` stacks its responses as brms's
  `rep.meas` factor. See `?frmtmb-emmeans`.
- emmeans now shows frmtmb's reason when it refuses, instead of
  "Perhaps a 'data' or 'params' argument is needed".
- Fixed silent wrong answers in `emmeans()`: `dpar =` was ignored and
  `mu` returned; `type = "response"` did not apply the inverse link;
  `s()` and `t2()` smooths were left out of the means; `mo()` terms
  gave NaN standard errors; `y | weights(w)` was read as a response
  transformation; `re_formula` was ignored.
