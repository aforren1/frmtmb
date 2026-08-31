# Changelog

## frmtmb 0.13.0

Method-surface (“sugar”) milestone: the conventional S3 methods that
downstream packages and user muscle memory dispatch on.

- New accessors: [`sigma()`](https://rdrr.io/r/stats/sigma.html),
  [`terms()`](https://rdrr.io/r/stats/terms.html),
  [`weights()`](https://rdrr.io/r/stats/weights.html),
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html),
  [`deviance()`](https://rdrr.io/r/stats/deviance.html),
  [`extractAIC()`](https://rdrr.io/r/stats/extractAIC.html) (enables
  [`step()`](https://rdrr.io/r/stats/step.html) /
  [`drop1()`](https://rdrr.io/r/stats/add1.html)),
  [`ngrps()`](../reference/ngrps.md),
  [`prior_summary()`](../reference/prior_summary.md), and a
  [`profile()`](https://rdrr.io/r/stats/profile.html) method wrapping
  [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html).
- **Breaking**: [`coef()`](https://rdrr.io/r/stats/coef.html) now
  follows the lme4 / glmmTMB / brms convention (per-group coefficients:
  fixed effects broadcast over levels plus the conditional modes). Use
  [`fixef()`](../reference/fixef.md) for the fixed effects alone. Fits
  without random effects still return the coefficient vector.
- [`conditional_effects()`](../reference/conditional_effects.md) with a
  base-graphics [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  method: grid predictions per predictor (or `"x:z"` pair) with Wald
  bands, smooths included, matrix covariates held at column means.
- `plot(fit)`: Pearson-residual diagnostics (residuals vs fitted, normal
  QQ).
- `pp_check()` (registered on the bayesplot generic): predictive checks
  from [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws
  through any `ppc_*` function.
- [`hypothesis()`](../reference/hypothesis.md): tests of arbitrary
  expressions of the parameters, brms spelling, with three methods:
  delta-method Wald (default), `"profile"` (profile-likelihood intervals
  via [`TMB::tmbroot()`](https://rdrr.io/pkg/TMB/man/tmbroot.html)
  lincombs; linear hypotheses, ML fits), and `"boot"` (parametric
  bootstrap percentile intervals; any expression). The expression
  environment includes natural-scale random-effect names
  (`sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`, `sigma`), so ICC
  inference is
  `hypothesis(fit, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)", method = "boot")`.
  The result is a `frmtmb_hypothesis` data frame carrying the bootstrap
  draws or profile curves in attributes, with a
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method
  (bootstrap histogram, profile curve, or the implied Wald normal, per
  hypothesis).
- `frm_bootstrap(fit, FUN, nsim)`: parametric bootstrap over
  [`refit()`](../reference/refit.md)s (the `bootMer` analog), with
  [`confint()`](https://rdrr.io/r/stats/confint.html) percentile
  intervals.
- `refit(fit, newresp)`: refit to a new response reusing the assembled
  design with a warm start - one re-tape plus optimization; the engine
  for parametric bootstrap over
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws.
- The accessor set makes insight’s default methods work
  (`find_response()`, `get_data()`, `get_sigma()`, `n_obs()`, …);
  dedicated insight methods (random-effect formula parts,
  `get_parameters()`) are a roadmap item.

## frmtmb 0.12.0

Pre-release consolidation toward CRAN. Highlights by milestone:

### Model grammar (brms-compatible spelling)

- [`bf()`](../reference/bf.md) with distributional-parameter formulas
  carrying the full predictor grammar (random effects and smooths in any
  dpar), constant dpars, and addition terms
  [`weights()`](https://rdrr.io/r/stats/weights.html), `trials()`,
  `cens()` (left / right / interval),
  [`trunc()`](https://rdrr.io/r/base/Round.html) (with the inclusive
  discrete correction).
- Nonlinear formulas (`nl = TRUE`) with a full linear predictor per
  parameter.
- Multivariate models: [`mvbf()`](../reference/mvbf.md) / `bf() + bf()`
  / `mvbind()`, a family per response, gaussian residual correlation
  ([`set_rescor()`](../reference/mvbf.md)), and `|ID|` random-effect
  correlation across formulas.
- mgcv smooths `s()` / `t2()` in any linear predictor, including
  matrix-covariate summation-convention terms: scalar-on-function,
  function-on-scalar, and function-on-function regression validated
  against mgcv.
- Covariance structures: `us`, `diag`, `homdiag`, `cs`, `ar1`, `toep`,
  `ou` (positions via [`num_factor()`](../reference/num_factor.md)), and
  known-structure terms `gr(g, cov = A)` (dense, with correlated slopes
  via a Kronecker product) and `gr(g, prec = Q)` (sparse GMRF).

### Families

gaussian, poisson, binomial, Gamma, lognormal, student, negbinomial,
nbinom1, beta, tweedie, compois, beta_binomial, skew_normal,
inverse.gaussian, exgaussian, zero-inflated and hurdle counts,
cumulative ordinal, matrix-response multinomial, and
[`custom_family()`](../reference/frmtmb_family.md) as a plain R
log-density with
[`check_custom_family()`](../reference/check_custom_family.md) AD
verification.

### Estimation and inference

- ML and REML by Laplace approximation; `quadrature = TRUE` upgrades
  scalar random effects to adaptive Gauss-Kronrod marginalization
  (matches `glmer(nAGQ = 25)` and GLMMadaptive).
- MAP / regularized ML via brms-style
  [`set_prior()`](../reference/set_prior.md); hard bounds (`lower` /
  `upper`, and per-prior `lb` / `ub`).
- `sdreport` is deferred until standard errors are first needed (roughly
  a quarter off fit time; `se = TRUE` restores eager mode).
- Pluggable optimizers (`frmtmb_control(optimizer = )`), including
  arbitrary user functions.
- [`confint()`](https://rdrr.io/r/stats/confint.html) (Wald / profile /
  likelihood-root), natural-scale
  [`confint_varcorr()`](../reference/confint_varcorr.md),
  [`anova()`](https://rdrr.io/r/stats/anova.html) likelihood-ratio
  tests, [`diagnose()`](../reference/diagnose.md).

### Diagnostics and ecosystem

- Residuals: response, pearson, and one-step-ahead (`type = "osa"`).
- DHARMa simulation-based residuals
  ([`dharma_residuals()`](../reference/dharma_residuals.md)).
- NUTS on the fitted objective:
  [`frm_sample()`](../reference/frm_sample.md) (ML-mode initialization,
  named draws, priors, bounds) and
  [`check_laplace()`](../reference/check_laplace.md); bayesplot works on
  the draws.
- emmeans and marginaleffects support, both verified against glmmTMB.

### Verification

Roughly 500 tests compare fits against glmmTMB, lme4, mgcv, MASS,
survival, nnet, GLMMadaptive, and hand-written RTMB references, plus an
edge-case suite mined from the lme4 / glmmTMB / brms issue trackers
(dev/test-backlog.md).
