# frmtmb

frmtmb fits regression models that you specify with a brms-style
formula grammar. Estimation is maximum likelihood, with the Laplace
approximation for latent effects. The backend is
[RTMB](https://cran.r-project.org/package=RTMB). The package
generates each model's objective as an R closure and differentiates
it on the TMB AD tape. There is no MCMC, no Stan, and no compilation
at run time. A refit re-tapes in milliseconds. The bootstrap,
influence, and multi-start machinery build on that speed.

Documentation: <https://aforren1.github.io/frmtmb/>

## Example

```r
library(frmtmb)
data(sleepstudy, package = "lme4")

fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
           data = sleepstudy)
summary(fit)

# distributional regression: model the residual SD too
fit2 <- frm(bf(Reaction ~ Days + (Days | Subject),
               sigma ~ Days + (1 | Subject)) + gaussian(),
            data = sleepstudy)
anova(fit, fit2)

hypothesis(fit, "sd_Subject__Days^2 / (sd_Subject__Days^2 + sigma^2)",
           method = "boot")
```

## Status

Pre-release (v0.30). The goal is a CRAN release. Validation has
three layers:

- The suite contains about 2700 tests. Every model class is compared
  with an exact reference: glmmTMB, lme4, mgcv, MASS, survival,
  nnet, GLMMadaptive, quantreg, mice, closed-form marginals, or
  hand-written ML.
- The model-building layer is compared with brms itself. Design
  matrices, random-effect structures, and special-term data agree
  with `brms::make_standata()` to near machine precision. An opt-in
  tier verifies that our estimates equal the mode of the Stan
  programs that brms generates.
- A pairwise grammar fuzzer sweeps feature combinations against
  metamorphic invariants. The resulting compatibility map is
  queryable with `frm_compat()` and is published as the
  [feature compatibility](https://aforren1.github.io/frmtmb/articles/compatibility.html)
  article.

## Model grammar

- `bf()` uses the brms spelling. Families attach with `+`. To port
  brms code, remove the priors and change `brm()` to `frm()`.
- Every distributional parameter can take its own formula with the
  full predictor grammar
  (`bf(y ~ s(x) + (1 | g), sigma ~ s(z) + (1 | g))`), or a constant.
- Random effects use the lme4 syntax and add structured covariances:
  `us`, `diag`, `homdiag`, `cs`, `homcs`, `ar1`, `hetar1`, `toep`,
  `homtoep`, `ou`, spatial `exp`/`gau`/`mat` over `num_factor(x, y)`
  coordinates, reduced-rank `rr(d =)`, and known structure
  `gr(g, cov = A)` / `gr(g, prec = Q)` / `equalto()`. The `|ID|`
  syntax correlates effects across formulas.
- `car(M, gr = g, type =)` fits spatial Gaussian Markov random
  fields with brms's spelling and all four of its types (`escar`,
  `esicar`, `icar`, and `bym2` with brms's scaling convention).
  `spde(fem, gr = node)` fits a Matern field over a finite-element
  mesh, with `fmesher`/INLA matrices as fixed data. The package
  assembles precisions on the tape from fixed sparse matrices.
  Every normalizing constant is analytic.
- mgcv smooths `s()`/`t2()` work in any linear predictor.
  Matrix-covariate terms give scalar-on-function,
  function-on-scalar, and function-on-function regression.
- `gp(x)` fits an exact Gaussian process and predicts with kriging.
  `gp(x, k = 30)` uses the Hilbert-space approximation with brms's
  exact input convention, so the same call is the same
  approximation. Up to three dimensions are supported, with
  per-dimension lengthscales; `iso = TRUE` shares one.
- `mo()` fits monotonic effects. `mi()` imputes continuous
  predictors in one step, and `mi(sdx)` models measurement error.
  `cs()` fits category-specific ordinal effects. `mo()` and `mi()`
  support two-way interactions.
- Multivariate models use `mvbf`, per-response families, and
  `rescor`. Nonlinear formulas use `nl = TRUE`. `multinomial(K)`
  takes a matrix response.
- The families include the usual GLM(M) set plus student, tweedie,
  compois, beta-binomial, skew-normal, ex-gaussian, weibull, shifted
  lognormal, and quantile regression (`asym_laplace`), with
  zero-inflated and hurdle variants and four ordinal families.
  `mixture()` fits finite mixtures, including group-level latent
  classes (`mixture(..., groups = ~g)`) with class-specific random
  effects for growth-mixture models. `mixture_mvn(K, D, model =)`
  fits multivariate gaussian components over mclust's covariance
  taxonomy, with covariate-dependent means and gating.
  `custom_family()` takes a plain R log-density. The test suite fits
  a Wiener drift-diffusion model in about 15 lines.
- The addition terms are `weights()`, `trials()` (counts or
  proportions), `cens()`, `trunc()`, `se()` (meta-analysis), `mi()`,
  and `vint()`/`vreal()`.

## Estimation and inference

- Estimation is ML or REML. `quadrature = TRUE` uses adaptive
  quadrature for scalar random effects. `set_prior()` gives MAP
  estimation with the brms spelling. Hard bounds and pluggable
  optimizers are available. `frmtmb_control()` offers
  `profile = TRUE` for many-coefficient models, `sparse_x = TRUE`
  for many-level fixed factors, `autoscale = TRUE` for badly scaled
  predictors, and `verbose =` for timed progress on slow fits.
- `confint()` offers Wald, profile, likelihood-root, and bootstrap
  intervals. `hypothesis()` tests expressions of coefficients, with
  natural-scale `sd_`/`cor_` names for ICC-type quantities.
  `anova()` and `drop1()` give likelihood-ratio tests.
  `frm_bootstrap()`, `influence()` with `cooks.distance()` and
  `dfbetas()`, `frm_allfit()`, and `frm_multiple()` (Rubin pooling,
  including variance components and hypotheses) complete the set.
- Simulation goes both ways. `simulate()` draws from a fit.
  `frm_simulate()` draws from a bare design, with natural-scale
  parameters (`sd_g__Intercept = 0.7`) or with parameters drawn per
  replicate from `set_prior()` specifications. That is the
  `sample_prior = "only"` prior-predictive workflow, without MCMC.
- `frm_sample()` runs NUTS on the fitted objective and returns a
  full draws surface (`posterior_epred()`, `posterior_predict()`,
  `hypothesis()`, `pp_check()`). `check_laplace()` audits the
  approximation.
- Diagnostics include one-step-ahead residuals calibrated for
  censored, truncated, and ordinal responses, deviance residuals
  across the GLM families, and response-scale `se.fit` for every
  family through the joint delta method. DHARMa, `pp_check()`,
  `plot()`, and `conditional_effects()` are supported. Hooks for
  emmeans, marginaleffects, and insight/easystats are registered.

## Why

brms defined a formula grammar that covers a wide model space.
frmtmb brings that grammar to maximum likelihood. The package
generates each model's objective from the formula; it does not
select from a fixed set of likelihoods. Because of that, features
that would each need a bespoke implementation become ordinary code
paths, and they compose with each other: random effects in any
distributional parameter, nonlinear predictors, per-response
families, monotonic effects, in-model imputation, latent-class
mixtures, and custom families written as plain R log-densities.

Related work: [glmmTMB](https://glmmtmb.github.io/glmmTMB/) is the
closest relative, a mature and fast TMB-based mixed-model package.
frmtmb matches its fits where the models overlap and follows its
conventions in several places.
[BayesRTMB](https://github.com/norimune/BayesRTMB) is a
Bayesian-first, Stan-like modeling layer on the same RTMB backend.
[brms](https://paulbuerkner.com/brms/) itself is the right tool when
you want priors and full posterior inference; the shared grammar is
meant to make movement between the two easy.

## Life cycle

frmtmb is maturing. The package is not yet on CRAN.

- **The model grammar is stable.** It follows brms, so it changes
  only when brms changes. Formulas that fit today will fit in later
  versions.
- **The fitted-object API is stable.** `frm()`, the accessor methods
  (`coef()`, `confint()`, `vcov()`, `predict()`, and the rest), and
  the family constructors keep their current behavior.
- **Internal structure can change.** Fields of the fitted object
  that no exported method reaches are not part of the API. Use the
  accessors.
- **Some parts are still moving.** Multivariate coverage of the
  post-fit methods, the mixture families, and `frm_sample()` gain
  features between releases.

Version numbers stay below 1.0 until the CRAN release. Breaking
changes are listed in `NEWS.md`.

The package is actively developed and maintained. Planned work: full
multivariate support in every post-fit method, more reference
comparisons, and CRAN submission.

## Citation

Run `citation("frmtmb")` in R. Please also cite RTMB, TMB, and tmbstan.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Everyone taking part must
follow the [Code of Conduct](CODE_OF_CONDUCT.md).

See [SPEC.md](SPEC.md) for the design and `NEWS.md` for the
changelog.
