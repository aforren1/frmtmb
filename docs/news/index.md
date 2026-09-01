# Changelog

## frmtmb 0.20.0

- `mixture_mvn(K, D)`: multivariate gaussian mixture components for
  model-based clustering of an n x D matrix response (the mclust /
  clustTMB use case). Every class mean is a full linear predictor, so
  cluster means may depend on covariates and random effects; mixing
  weights are multinomial-logit dpars with their own formulas. Validated
  against a hand-rolled ML fit to 1.3e-10 and against faithful-data
  cluster recovery. Limitations (documented in
  [`?mixture_mvn`](../reference/mixture_mvn.md)): class covariances are
  unstructured (`us`) and covariate-free - no mclust-style constrained
  covariance taxonomy (EII..VEV);
  `cens()`/[`trunc()`](https://rdrr.io/r/base/Round.html),
  [`mvbf()`](../reference/mvbf.md), and
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) are not
  supported.
- Multi-dimensional Gaussian processes: `gp(x1, x2, ...)` takes up to
  three variables, with a separate lengthscale per dimension by default
  (the brms convention) or one shared lengthscale via `iso = TRUE`, in
  both the exact and the Hilbert-space (`k =`) form. The tensor HSGP
  basis is capped at `k^D <= 1000` columns.
  [`confint_varcorr()`](../reference/confint_varcorr.md) reports one
  range row per dimension. Validated against a direct 2-D ML fit to
  3.8e-11.
- Kriging prediction for the exact `gp()`: `predict(newdata = )` at
  positions not seen at fit time now returns the GP conditional mean,
  and `se.fit = TRUE` adds the conditional variance, instead of
  erroring. Kriging mean validated to 1.7e-15 and the standard error to
  the exact decomposition identity. The Hilbert-space form already
  interpolated through its basis.
- `frmtmb_control(sparse_x = TRUE)`: sparse fixed-effect design matrices
  ([`Matrix::sparse.model.matrix`](https://rdrr.io/pkg/Matrix/man/sparse.model.matrix.html)),
  the analog of `glmmTMB(sparseX = )`, for models where a many-level
  fixed factor dominates memory. Estimates are identical to the dense
  path (gradient agreement 3e-14; 13.8% lower tape memory in the
  validation model);
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) on the
  fit returns a sparse matrix. Prediction falls back to the dense
  builder for newdata containing NA factor rows, which
  `sparse.model.matrix` would silently zero.
- `frmtmb_control(autoscale = TRUE)`: internal standardization of badly
  scaled continuous predictors (the lme4 \>= 1.1.37 feature): a scaled
  pre-fit is mapped back exactly and warm-starts the reported fit,
  parameter magnitudes feed `nlminb`’s scale hook, and convergence and
  standard-error Hessian checks run in scaled units. Turns 12 scaling
  warnings into 0 on the validation problem with a 0.0 log-likelihood
  difference against a manually standardized reference. Reported results
  are always on the original scale.
- [`frm_multiple()`](../reference/frm_multiple.md) pooling extensions:
  `$pooled_varcorr` pools random-effect SDs and correlations across
  imputations on transformed scales (log for SDs and GP ranges, Fisher z
  for correlations) with Barnard-Rubin degrees of freedom, and
  [`hypothesis()`](../reference/hypothesis.md) on a `frmtmb_multiple`
  object pools delta-method hypothesis estimates by Rubin’s rules with
  t-based intervals. Fixed-effect pooling matches
  [`mice::pool()`](https://amices.org/mice/reference/pool.html) to
  2.1e-6. Likelihood-ratio pooling across imputations (mice’s D3) is not
  implemented; [`anova()`](https://rdrr.io/r/stats/anova.html) on pooled
  fits remains per-fit.

## frmtmb 0.19.0

- Latent classes now combine with continuous random effects
  (growth-mixture models): `mixture(..., groups = ~g)` accepts random
  effects, smooths, and gp() terms in the component formulas. The class
  sum is computed conditional on the latent effects, so one Laplace
  approximation integrates them - the sum-and-integral swap makes this
  exact as an integrand identity, and even crossed random effects are
  structurally fine. Random effects written in a component formula are
  class-specific by construction. Accuracy: the Laplace approximation of
  the class-mixture integrand carries a small documented bias (about 0.1
  log-likelihood units in the validation problem, with parameter
  estimates matching the exact closed-form marginal to 0.01);
  `quadrature = TRUE` is numerically exact when the per-group integrand
  is univariate (validated to 5e-3) and approximate when class-specific
  intercepts couple.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) now draws one
  class per group and then simulates each observation from its group’s
  component; [`mixture_probs()`](../reference/mixture_probs.md) gives
  empirical-Bayes classification conditional on the modes.
- [`variables()`](../reference/variables.md) (brms spelling): lists
  every parameter name usable in
  [`hypothesis()`](../reference/hypothesis.md) expressions
  (coefficients, `sd_`/`cor_` summaries, `sigma`); on
  [`frm_sample()`](../reference/frm_sample.md) output it lists the draw
  columns.
- [`get_prior()`](../reference/get_prior.md) (brms spelling): enumerates
  every slot [`set_prior()`](../reference/set_prior.md) can target
  (class/coef/dpar/group rows), from a formula plus data or from a fit;
  the default in every slot is flat.
- README rewritten for the current scope, with related-work pointers
  (glmmTMB, BayesRTMB); getting-started vignette extended to the full
  grammar; SPEC.md status and deviations brought current.

## frmtmb 0.18.0

The last three roadmap features.

- `gp(x)` Gaussian-process terms: exact (a dense squared-exponential
  block over the unique positions, with a standard 1e-6 nugget) or the
  Hilbert-space approximation via `gp(x, k = 30)` (sine basis with
  spectral-density prior SDs), in any linear predictor. Exact fits match
  direct GP marginal ML; the approximation converges to the exact answer
  and predicts at arbitrary positions (the exact form predicts at
  observed positions only).
  [`confint_varcorr()`](../reference/confint_varcorr.md) reports
  `sd(gp)` and `range(gp)`.
- `mo()` and `mi()` interactions: `mo(x):z`, `mo(x)*z`, `mi(x):z`,
  `mi(x)*z` with numeric multipliers; `mo()` interactions share their
  variable’s simplex (brms convention). Both validated exact against
  direct ML (the `mi()` interaction stays linear in the latent value, so
  the closed-form marginal still applies).
- Group-level latent-class mixtures: `mixture(..., groups = ~g)` sums
  the class assignment per group (growth-mixture / latent-class
  regression; the tractable nested case of brms#1905). Exact against
  direct ML; [`mixture_probs()`](../reference/mixture_probs.md) returns
  posterior class probabilities per group, or per observation for
  ordinary mixtures. Restrictions: no random effects or smooths
  alongside, mixing predictors evaluated at each group’s first row, no
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) yet.
  Crossed-design group mixtures remain out of the Laplace class.

## frmtmb 0.17.0

Clearing the roadmap: mixtures, measurement error, and the remaining
grammar gaps.

- `mixture(fam1, fam2, ...)` finite-mixture families: each component
  keeps its own suffixed dpars (`mu1`, `sigma1`, …), mixing proportions
  are multinomial-logit dpars with their own linear predictors
  (mixture-of-experts comes free), and the likelihood is a branch-free
  logsumexp - so random effects and smooths work in any component
  formula. Exact against direct ML for gaussian and poisson mixtures.
  Component means initialize on spread response quantiles; the usual
  finite-mixture multimodality caveats apply
  ([`frm_allfit()`](../reference/frm_allfit.md), or order intercepts via
  bounds).
- `mi(sdx)` measurement error (brms `me()`): known per-observation
  measurement SDs make every true value latent, with the observed values
  entering through a measurement model. Exact against the closed-form
  bivariate-normal marginal; attenuation bias is corrected and NAs
  (missing + mismeasured) combine.
- `cs()` category-specific effects for `sratio`, `cratio`, and `acat`
  (refused for `cumulative`, as in brms). Exact against direct ML.
- `equalto(x + 0 | g, V)` covariance structure: fully fixed known V,
  zero parameters (meta-analytic sampling covariances).
- rr() fits now support `se.fit`: the delta method routes the factor
  columns through the loadings Jacobian and adds loading-parameter
  columns; verified parameterization-invariant against the
  full-rank-vs-`us()` equivalence.
- Registered insight methods (`find_formula`, `find_random`,
  `get_parameters`, `get_varcov`, `find_statistic`, link accessors):
  [`insight::is_mixed_model()`](https://easystats.github.io/insight/reference/is_mixed_model.html)
  and the easystats accessor layer now see the random-effect structure.

## frmtmb 0.16.0

The two deferred architecture features.

- `rr(x | g, d = k)` reduced-rank / factor-analytic covariance
  (glmmTMB’s GLVM structure): per-level coefficients are loadings times
  iid standard-normal factors, so `b` holds the factors while the design
  matrices span the coefficient space. Matches glmmTMB exactly and
  reduces to `us()` at full rank. `se.fit` on rr models is deferred (the
  loadings enter the coefficients nonlinearly); `predict`, `simulate`,
  `ranef`, and `VarCorr` (the rank-k covariance) all work.
- `mi()` one-step imputation for continuous predictors, reversing the
  original exclusion: `bf(y ~ mi(x) + z) + bf(x | mi() ~ z)` turns
  missing `x` values into latent parameters that Laplace integrates
  (exactly, in the linear-gaussian case - validated against the
  closed-form marginal likelihood). Imputation uncertainty propagates
  into the coefficient standard errors automatically. Rows with NAs in
  non-`mi()` variables still drop; imputation models must be gaussian or
  student; discrete predictors remain impossible (as in Stan). `me()`
  measurement error is the natural follow-on.

## frmtmb 0.15.0

Tier-2 sweep of dev/feature-gaps.md plus a method surface for sampled
fits.

- `mo()` monotonic effects (brms syntax): scale coefficient times an
  estimated simplex, exact against direct ML; works in prediction,
  [`conditional_effects()`](../reference/conditional_effects.md), and
  with ordered factors. No frequentist package offers this. Standalone
  additive terms only for now.
- Sequential ordinal families `sratio`, `cratio`, `acat`, validated
  against direct ML and collapsing to logistic regression at K = 2.
- Covariance structures: `hetar1`, `homcs`, `homtoep` (glmmTMB
  parameterizations, exact where glmmTMB converges) and spatial
  [`exp()`](https://rdrr.io/r/base/Log.html), `gau()`, `mat()` over
  `num_factor(x, y)` coordinates (Matern uses bounded internal
  transforms for stability where both we and glmmTMB otherwise diverge).
- Fixed a latent parser bug: `exp(x)` (or any covariance-structure name)
  used as a plain function in a fixed formula was silently stripped by
  the formula splitter - since v0.1. Such calls are now protected
  automatically.
- `influence(fit, groups = )` and
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html):
  case-deletion diagnostics over warm-started refits.
- [`frm_multiple()`](../reference/frm_multiple.md): fits across
  multiply-imputed datasets (list or
  [`mice::mids`](https://amices.org/mice/reference/mids.html)) pooled by
  Rubin’s rules with Barnard-Rubin df.
- [`conditional_effects()`](../reference/conditional_effects.md) gains
  `method = "predict"` (prediction intervals) and data-frame
  `conditions` (one condition set per row, brms style).
- `vint()`/`vreal()` addition terms pass arbitrary data vectors to
  custom families (brms custom-family convention). The test suite
  includes a full Wiener drift-diffusion model written as a
  [`custom_family()`](../reference/frmtmb_family.md) in plain R - the
  workflow brms needs raw Stan code for.
- Method surface for [`frm_sample()`](../reference/frm_sample.md) draws:
  [`summary()`](https://rdrr.io/r/base/summary.html) (with Rhat/ESS),
  [`fixef()`](../reference/fixef.md),
  [`VarCorr()`](../reference/VarCorr.md) (natural scale),
  [`hypothesis()`](../reference/hypothesis.md) (exact posterior
  version), [`posterior_epred()`](../reference/posterior_epred.md) /
  [`posterior_predict()`](../reference/posterior_epred.md) (each draw
  runs the full prediction machinery), `pp_check()`, and
  [`posterior::as_draws()`](https://mc-stan.org/posterior/reference/draws.html).
- Examples added across the reference (run by R CMD check).

## frmtmb 0.14.0

Gap-closing milestone from a sweep of the lme4 / glmmTMB / brms
vignettes (dev/feature-gaps.md).

- `se()` addition term: known per-observation sampling SDs
  (meta-analysis), gaussian and student. `se(x)` fixes the residual SD
  (sigma is mapped out); `se(x, sigma = TRUE)` adds an estimated sigma
  in quadrature. With `(1 | obs)` this is random-effects meta-analysis,
  and with `gr(g, cov = A)` the phylogenetic version.
- Binomial (and beta-binomial) responses may be proportions of
  `trials()`, the glm/glmer idiom, converted to counts internally.
- New families: `bernoulli`, `geometric`, `exponential`, `weibull`
  (mean-parameterized, matches survreg’s likelihood),
  `shifted_lognormal`, `hurdle_gamma`, `hurdle_lognormal`,
  `zero_inflated_binomial`, `zero_inflated_beta` (validated against
  glmmTMB), and `asym_laplace` for quantile regression (fixed `quantile`
  dpar reproduces
  [`quantreg::rq`](https://rdrr.io/pkg/quantreg/man/rq.html) estimates).
- `ranef(condVar = TRUE)`: conditional SDs of the modes (TMB/glmmTMB
  convention: fixed-effect uncertainty propagates), plus tidy
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) methods
  for `ranef` and `VarCorr` output.
- [`frm_allfit()`](../reference/frm_allfit.md): refit under every
  available optimizer (nlminb, optim, bobyqa, NLopt) and compare, the
  lme4 `allFit()` analog.
- [`frm_simulate()`](../reference/frm_simulate.md): de novo simulation
  from a [`bf()`](../reference/bf.md) formula and supplied parameters
  with no fitted model (power analysis; the glmmTMB `simulate_new()`
  analog).
- `frmtmb_control(profile = TRUE)` (experimental): profile the
  primary-coefficient vector into the inner Laplace problem, the
  `glmmTMBControl(profile = TRUE)` / `glmer(nAGQ = 0)` speed
  approximation for many-coefficient models; covariance comes from the
  joint precision.
- [`anova()`](https://rdrr.io/r/stats/anova.html) no longer fails on
  models sharing a primary formula (distributional vs plain fits); row
  labels now include dpar formulas.

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
