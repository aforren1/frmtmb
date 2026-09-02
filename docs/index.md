# frmtmb

[![R-CMD-check](https://github.com/aforren1/frmtmb/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/aforren1/frmtmb/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/test-coverage.yaml)
[![pkgcheck](https://github.com/aforren1/frmtmb/actions/workflows/pkgcheck.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/pkgcheck.yaml)

frmtmb fits regression models that you specify with a brms-style formula
grammar. Estimation is maximum likelihood, with the Laplace
approximation for latent effects. The backend is
[RTMB](https://cran.r-project.org/package=RTMB). The package generates
each model’s objective as an R closure and differentiates it on the TMB
AD tape. There is no MCMC, no Stan, and no compilation at run time. A
refit re-tapes in milliseconds. The bootstrap, influence, and
multi-start machinery build on that speed.

Documentation: <https://aforren1.github.io/frmtmb/>

## Example

``` r

library(frmtmb)
data(sleepstudy, package = "lme4")

fit <- frm(bf(Reaction ~ Days + (Days | Subject)), family = gaussian(),
           data = sleepstudy)
summary(fit)

# distributional regression: model the residual SD too
fit2 <- frm(bf(Reaction ~ Days + (Days | Subject),
               sigma ~ Days + (1 | Subject)), family = gaussian(),
            data = sleepstudy)
anova(fit, fit2)

hypothesis(fit, "sd_Subject__Days^2 / (sd_Subject__Days^2 + sigma^2)",
           method = "boot")
```

## Status

Pre-release (v0.35). The goal is a CRAN release. Validation has three
layers:

- The suite contains about 4000 tests. Every model class is compared
  with an exact reference: glmmTMB, lme4, mgcv, MASS, survival, nnet,
  GLMMadaptive, quantreg, mice, closed-form marginals, or hand-written
  ML.
- The model-building layer is compared with brms itself. Design
  matrices, random-effect structures, and special-term data agree with
  [`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html)
  to near machine precision. An opt-in tier verifies that our estimates
  equal the mode of the Stan programs that brms generates.
- A pairwise grammar fuzzer sweeps feature combinations against
  metamorphic invariants. The resulting compatibility map is queryable
  with
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md)
  and is published as the [feature
  compatibility](https://aforren1.github.io/frmtmb/articles/compatibility.html)
  article.

## Model grammar

- [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) uses the
  brms spelling, and the family is a separate argument:
  `frm(bf(y ~ x), family = poisson(), data = d)`. With no family,
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) fits
  gaussian, as brms and glmmTMB do. To port brms code, remove the priors
  and change `brm()` to
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md); a
  measured audit of the brms vignettes puts about 7 of 10 of their model
  calls through that mechanical transform unchanged.
- Every distributional parameter can take its own formula with the full
  predictor grammar (`bf(y ~ s(x) + (1 | g), sigma ~ s(z) + (1 | g))`),
  or a constant.
- Random effects use the lme4 syntax and add structured covariances:
  `us`, `diag`, `homdiag`, `cs`, `homcs`, `ar1`, `hetar1`, `toep`,
  `homtoep`, `ou`, spatial `exp`/`gau`/`mat` over `num_factor(x, y)`
  coordinates, reduced-rank `rr(d =)`, and known structure
  `gr(g, cov = A)` / `gr(g, prec = Q)` / `equalto()`. The `|ID|` syntax
  correlates effects across formulas, including known-matrix blocks (a
  multi-trait animal model is one Kronecker block). Multi-membership
  terms use `mm(g1, g2)` with weights and member-specific `mmc()`
  covariates.
- Within-group residual correlation uses brms’s R-side terms:
  [`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`, `cosy()`,
  and `unstr()`, for gaussian and student responses, validated against
  [`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html) under ML and
  REML.
- `car(M, gr = g, type =)` fits spatial Gaussian Markov random fields
  with brms’s spelling and all four of its types (`escar`, `esicar`,
  `icar`, and `bym2` with brms’s scaling convention).
  `spde(fem, gr = node)` fits a Matern field over a finite-element mesh,
  with `fmesher`/INLA matrices as fixed data. The package assembles
  precisions on the tape from fixed sparse matrices. Every normalizing
  constant is analytic.
- mgcv smooths `s()`/`t2()` work in any linear predictor.
  Matrix-covariate terms give scalar-on-function, function-on-scalar,
  and function-on-function regression.
- `gp(x)` fits an exact Gaussian process and predicts with kriging.
  `gp(x, k = 30)` uses the Hilbert-space approximation with brms’s exact
  input convention, so the same call is the same approximation. Up to
  three dimensions are supported, with per-dimension lengthscales;
  `iso = TRUE` shares one.
- `mo()` fits monotonic effects. `mi()` imputes continuous predictors in
  one step, and `mi(sdx)` models measurement error. `cs()` fits
  category-specific ordinal effects. `mo()` and `mi()` support two-way
  interactions.
- Multivariate models use `mvbf`, per-response families, and `rescor`.
  Nonlinear formulas use `nl = TRUE`. `multinomial(K)` takes a matrix
  response.
- The families include the usual GLM(M) set plus student, tweedie,
  compois, beta-binomial, skew-normal, ex-gaussian, weibull, shifted
  lognormal, and quantile regression (`asym_laplace`), with
  zero-inflated and hurdle variants and four ordinal families.
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  fits finite mixtures, including group-level latent classes
  (`mixture(..., groups = ~g)`) with class-specific random effects for
  growth-mixture models. `mixture_mvn(K, D, model =)` fits multivariate
  gaussian components over mclust’s covariance taxonomy, with
  covariate-dependent means and gating.
  [`categorical()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  fits nominal responses,
  [`von_mises()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  circular ones, and
  [`cox()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
  proportional hazards with a flexible baseline and Laplace frailties.
  `hmm(K)` fits hidden Markov models with covariate-dependent
  transitions and forward-backward decoding
  ([`hmm_probs()`](https://aforren1.github.io/frmtmb/reference/hmm_probs.md),
  [`hmm_viterbi()`](https://aforren1.github.io/frmtmb/reference/hmm_viterbi.md));
  `lca(K)` fits poLCA-style latent class analysis with class-membership
  regression.
  [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  takes a plain R log-density. The test suite fits a Wiener
  drift-diffusion model in about 15 lines.
- Ordinary differential equations:
  [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
  solves compartment models inside nonlinear formulas (population
  pharmacokinetics), with repeated dosing, infusions, and estimated
  bioavailability.
- The addition terms are
  [`weights()`](https://rdrr.io/r/stats/weights.html), `trials()`
  (counts or proportions), `cens()`,
  [`trunc()`](https://rdrr.io/r/base/Round.html), `se()`
  (meta-analysis), `mi()`, and `vint()`/`vreal()`.

## Estimation and inference

- Estimation is ML or REML. `quadrature = TRUE` uses adaptive quadrature
  for scalar random effects.
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  gives MAP estimation with the brms spelling. Hard bounds and pluggable
  optimizers are available.
  [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md)
  offers `profile = TRUE` for many-coefficient models, `sparse_x = TRUE`
  for many-level fixed factors, `autoscale = TRUE` for badly scaled
  predictors, and `verbose =` for timed progress on slow fits.
- [`confint()`](https://rdrr.io/r/stats/confint.html) offers Wald,
  profile, likelihood-root, and bootstrap intervals.
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  tests expressions of coefficients, with natural-scale `sd_`/`cor_`
  names for ICC-type quantities.
  [`anova()`](https://rdrr.io/r/stats/anova.html) and
  [`drop1()`](https://rdrr.io/r/stats/add1.html) give likelihood-ratio
  tests.
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md),
  [`influence()`](https://rdrr.io/r/stats/lm.influence.html) with
  [`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html)
  and [`dfbetas()`](https://rdrr.io/r/stats/influence.measures.html),
  [`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md),
  and
  [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  (Rubin pooling, including variance components and hypotheses) complete
  the set.
- Simulation goes both ways.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws from a
  fit.
  [`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
  draws from a bare design, with natural-scale parameters
  (`sd_g__Intercept = 0.7`) or with parameters drawn per replicate from
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  specifications. That is the `sample_prior = "only"` prior-predictive
  workflow, without MCMC.
- [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  runs NUTS on the fitted objective and returns a full draws surface
  ([`posterior_epred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md),
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md),
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md),
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)).
  [`check_laplace()`](https://aforren1.github.io/frmtmb/reference/check_laplace.md)
  audits the approximation.
- Diagnostics include one-step-ahead residuals calibrated for censored,
  truncated, and ordinal responses, deviance residuals across the GLM
  families, and response-scale `se.fit` for every family through the
  joint delta method. DHARMa,
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  are supported. Hooks for emmeans, marginaleffects, and
  insight/easystats are registered.

## Why

brms defined a formula grammar that covers a wide model space. frmtmb
brings that grammar to maximum likelihood. The package generates each
model’s objective from the formula; it does not select from a fixed set
of likelihoods. Because of that, features that would each need a bespoke
implementation become ordinary code paths, and they compose with each
other: random effects in any distributional parameter, nonlinear
predictors, per-response families, monotonic effects, in-model
imputation, latent-class mixtures, and custom families written as plain
R log-densities.

Related work: [glmmTMB](https://glmmtmb.github.io/glmmTMB/) is the
closest relative, a mature and fast TMB-based mixed-model package.
frmtmb matches its fits where the models overlap and follows its
conventions in several places.
[BayesRTMB](https://github.com/norimune/BayesRTMB) is a Bayesian-first,
Stan-like modeling layer on the same RTMB backend.
[brms](https://paulbuerkner.com/brms/) itself is the right tool when you
want priors and full posterior inference; the shared grammar is meant to
make movement between the two easy.

## Life cycle

frmtmb is maturing. The package is not yet on CRAN.

- **The model grammar is stable.** It follows brms, so it changes only
  when brms changes. Formulas that fit today will fit in later versions.
- **The fitted-object API is stable.**
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md), the
  accessor methods ([`coef()`](https://rdrr.io/r/stats/coef.html),
  [`confint()`](https://rdrr.io/r/stats/confint.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`predict()`](https://rdrr.io/r/stats/predict.html), and the rest),
  and the family constructors keep their current behavior.
- **Internal structure can change.** Fields of the fitted object that no
  exported method reaches are not part of the API. Use the accessors.
- **Some parts are still moving.** Multivariate coverage of the post-fit
  methods, the mixture families, and
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
  gain features between releases.

Version numbers stay below 1.0 until the CRAN release. Breaking changes
are listed in `NEWS.md`.

The package is actively developed and maintained. Planned work: full
multivariate support in every post-fit method, more reference
comparisons, and CRAN submission.

## Citation

Run `citation("frmtmb")` in R. Please also cite RTMB, TMB, and tmbstan.

## Contributing

See
[CONTRIBUTING.md](https://aforren1.github.io/frmtmb/CONTRIBUTING.md).
Everyone taking part must follow the [Code of
Conduct](https://aforren1.github.io/frmtmb/CODE_OF_CONDUCT.md).

See [SPEC.md](https://aforren1.github.io/frmtmb/SPEC.md) for the design
and `NEWS.md` for the changelog.
