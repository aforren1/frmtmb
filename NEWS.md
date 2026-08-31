# frmtmb 0.15.0

Tier-2 sweep of dev/feature-gaps.md plus a method surface for sampled
fits.

* `mo()` monotonic effects (brms syntax): scale coefficient times an
  estimated simplex, exact against direct ML; works in prediction,
  `conditional_effects()`, and with ordered factors. No frequentist
  package offers this. Standalone additive terms only for now.
* Sequential ordinal families `sratio`, `cratio`, `acat`, validated
  against direct ML and collapsing to logistic regression at K = 2.
* Covariance structures: `hetar1`, `homcs`, `homtoep` (glmmTMB
  parameterizations, exact where glmmTMB converges) and spatial
  `exp()`, `gau()`, `mat()` over `num_factor(x, y)` coordinates
  (Matern uses bounded internal transforms for stability where both
  we and glmmTMB otherwise diverge).
* Fixed a latent parser bug: `exp(x)` (or any covariance-structure
  name) used as a plain function in a fixed formula was silently
  stripped by the formula splitter - since v0.1. Such calls are now
  protected automatically.
* `influence(fit, groups = )` and `cooks.distance()`: case-deletion
  diagnostics over warm-started refits.
* `frm_multiple()`: fits across multiply-imputed datasets (list or
  `mice::mids`) pooled by Rubin's rules with Barnard-Rubin df.
* `conditional_effects()` gains `method = "predict"` (prediction
  intervals) and data-frame `conditions` (one condition set per row,
  brms style).
* `vint()`/`vreal()` addition terms pass arbitrary data vectors to
  custom families (brms custom-family convention). The test suite
  includes a full Wiener drift-diffusion model written as a
  `custom_family()` in plain R - the workflow brms needs raw Stan
  code for.
* Method surface for `frm_sample()` draws: `summary()` (with
  Rhat/ESS), `fixef()`, `VarCorr()` (natural scale), `hypothesis()`
  (exact posterior version), `posterior_epred()` /
  `posterior_predict()` (each draw runs the full prediction
  machinery), `pp_check()`, and `posterior::as_draws()`.
* Examples added across the reference (run by R CMD check).

# frmtmb 0.14.0

Gap-closing milestone from a sweep of the lme4 / glmmTMB / brms
vignettes (dev/feature-gaps.md).

* `se()` addition term: known per-observation sampling SDs
  (meta-analysis), gaussian and student. `se(x)` fixes the residual SD
  (sigma is mapped out); `se(x, sigma = TRUE)` adds an estimated sigma
  in quadrature. With `(1 | obs)` this is random-effects
  meta-analysis, and with `gr(g, cov = A)` the phylogenetic version.
* Binomial (and beta-binomial) responses may be proportions of
  `trials()`, the glm/glmer idiom, converted to counts internally.
* New families: `bernoulli`, `geometric`, `exponential`, `weibull`
  (mean-parameterized, matches survreg's likelihood),
  `shifted_lognormal`, `hurdle_gamma`, `hurdle_lognormal`,
  `zero_inflated_binomial`, `zero_inflated_beta` (validated against
  glmmTMB), and `asym_laplace` for quantile regression (fixed
  `quantile` dpar reproduces `quantreg::rq` estimates).
* `ranef(condVar = TRUE)`: conditional SDs of the modes (TMB/glmmTMB
  convention: fixed-effect uncertainty propagates), plus tidy
  `as.data.frame()` methods for `ranef` and `VarCorr` output.
* `frm_allfit()`: refit under every available optimizer (nlminb,
  optim, bobyqa, NLopt) and compare, the lme4 `allFit()` analog.
* `frm_simulate()`: de novo simulation from a `bf()` formula and
  supplied parameters with no fitted model (power analysis; the
  glmmTMB `simulate_new()` analog).
* `frmtmb_control(profile = TRUE)` (experimental): profile the
  primary-coefficient vector into the inner Laplace problem, the
  `glmmTMBControl(profile = TRUE)` / `glmer(nAGQ = 0)` speed
  approximation for many-coefficient models; covariance comes from
  the joint precision.
* `anova()` no longer fails on models sharing a primary formula
  (distributional vs plain fits); row labels now include dpar
  formulas.

# frmtmb 0.13.0

Method-surface ("sugar") milestone: the conventional S3 methods that
downstream packages and user muscle memory dispatch on.

* New accessors: `sigma()`, `terms()`, `weights()`, `model.matrix()`,
  `deviance()`, `extractAIC()` (enables `step()` / `drop1()`),
  `ngrps()`, `prior_summary()`, and a `profile()` method wrapping
  `TMB::tmbprofile()`.
* **Breaking**: `coef()` now follows the lme4 / glmmTMB / brms
  convention (per-group coefficients: fixed effects broadcast over
  levels plus the conditional modes). Use `fixef()` for the fixed
  effects alone. Fits without random effects still return the
  coefficient vector.
* `conditional_effects()` with a base-graphics `plot()` method: grid
  predictions per predictor (or `"x:z"` pair) with Wald bands, smooths
  included, matrix covariates held at column means.
* `plot(fit)`: Pearson-residual diagnostics (residuals vs fitted,
  normal QQ).
* `pp_check()` (registered on the bayesplot generic): predictive checks
  from `simulate()` draws through any `ppc_*` function.
* `hypothesis()`: tests of arbitrary expressions of the parameters,
  brms spelling, with three methods: delta-method Wald (default),
  `"profile"` (profile-likelihood intervals via `TMB::tmbroot()`
  lincombs; linear hypotheses, ML fits), and `"boot"` (parametric
  bootstrap percentile intervals; any expression). The expression
  environment includes natural-scale random-effect names
  (`sd_<group>__<term>`, `cor_<group>__<t1>__<t2>`, `sigma`), so ICC
  inference is `hypothesis(fit, "sd_g__Intercept^2 /
  (sd_g__Intercept^2 + sigma^2)", method = "boot")`. The result is a
  `frmtmb_hypothesis` data frame carrying the bootstrap draws or
  profile curves in attributes, with a `plot()` method (bootstrap
  histogram, profile curve, or the implied Wald normal, per
  hypothesis).
* `frm_bootstrap(fit, FUN, nsim)`: parametric bootstrap over
  `refit()`s (the `bootMer` analog), with `confint()` percentile
  intervals.
* `refit(fit, newresp)`: refit to a new response reusing the assembled
  design with a warm start - one re-tape plus optimization; the engine
  for parametric bootstrap over `simulate()` draws.
* The accessor set makes insight's default methods work
  (`find_response()`, `get_data()`, `get_sigma()`, `n_obs()`, ...);
  dedicated insight methods (random-effect formula parts,
  `get_parameters()`) are a roadmap item.

# frmtmb 0.12.0

Pre-release consolidation toward CRAN. Highlights by milestone:

## Model grammar (brms-compatible spelling)

* `bf()` with distributional-parameter formulas carrying the full
  predictor grammar (random effects and smooths in any dpar), constant
  dpars, and addition terms `weights()`, `trials()`, `cens()` (left /
  right / interval), `trunc()` (with the inclusive discrete correction).
* Nonlinear formulas (`nl = TRUE`) with a full linear predictor per
  parameter.
* Multivariate models: `mvbf()` / `bf() + bf()` / `mvbind()`, a family
  per response, gaussian residual correlation (`set_rescor()`), and
  `|ID|` random-effect correlation across formulas.
* mgcv smooths `s()` / `t2()` in any linear predictor, including
  matrix-covariate summation-convention terms: scalar-on-function,
  function-on-scalar, and function-on-function regression validated
  against mgcv.
* Covariance structures: `us`, `diag`, `homdiag`, `cs`, `ar1`, `toep`,
  `ou` (positions via `num_factor()`), and known-structure terms
  `gr(g, cov = A)` (dense, with correlated slopes via a Kronecker
  product) and `gr(g, prec = Q)` (sparse GMRF).

## Families

gaussian, poisson, binomial, Gamma, lognormal, student, negbinomial,
nbinom1, beta, tweedie, compois, beta_binomial, skew_normal,
inverse.gaussian, exgaussian, zero-inflated and hurdle counts,
cumulative ordinal, matrix-response multinomial, and `custom_family()`
as a plain R log-density with `check_custom_family()` AD verification.

## Estimation and inference

* ML and REML by Laplace approximation; `quadrature = TRUE` upgrades
  scalar random effects to adaptive Gauss-Kronrod marginalization
  (matches `glmer(nAGQ = 25)` and GLMMadaptive).
* MAP / regularized ML via brms-style `set_prior()`; hard bounds
  (`lower` / `upper`, and per-prior `lb` / `ub`).
* `sdreport` is deferred until standard errors are first needed
  (roughly a quarter off fit time; `se = TRUE` restores eager mode).
* Pluggable optimizers (`frmtmb_control(optimizer = )`), including
  arbitrary user functions.
* `confint()` (Wald / profile / likelihood-root), natural-scale
  `confint_varcorr()`, `anova()` likelihood-ratio tests, `diagnose()`.

## Diagnostics and ecosystem

* Residuals: response, pearson, and one-step-ahead (`type = "osa"`).
* DHARMa simulation-based residuals (`dharma_residuals()`).
* NUTS on the fitted objective: `frm_sample()` (ML-mode
  initialization, named draws, priors, bounds) and `check_laplace()`;
  bayesplot works on the draws.
* emmeans and marginaleffects support, both verified against glmmTMB.

## Verification

Roughly 500 tests compare fits against glmmTMB, lme4, mgcv, MASS,
survival, nnet, GLMMadaptive, and hand-written RTMB references, plus
an edge-case suite mined from the lme4 / glmmTMB / brms issue
trackers (dev/test-backlog.md).
