# frmtmb

[![R-CMD-check](https://github.com/aforren1/frmtmb/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/aforren1/frmtmb/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/test-coverage.yaml)
[![pkgcheck](https://github.com/aforren1/frmtmb/actions/workflows/pkgcheck.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/pkgcheck.yaml)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

frmtmb fits regression models that you specify with a brms-style
formula grammar. Estimation is maximum likelihood, with the Laplace
approximation for latent effects. The backend is
[RTMB](https://cran.r-project.org/package=RTMB). The package
generates each model's objective as an R closure and differentiates
it on the TMB AD tape. There is no MCMC, no Stan, and no compilation
at run time. A refit re-tapes in milliseconds. The bootstrap,
influence, and multi-start machinery build on that speed.

Documentation: <https://aforren1.github.io/frmtmb/>

## Statement of need

brms gave R a formula grammar that covers a very wide model space.
It reaches that space only through Bayesian estimation in Stan. A
user who wants the same grammar under maximum likelihood must divide
one model across several packages, accept a smaller model, or wait
for a sampler. frmtmb closes that gap. It fits the brms grammar by
maximum likelihood with the Laplace approximation, and it compiles
nothing at run time.

The gap is real because the existing frequentist packages each stop
at a different place.

- [lme4](https://cran.r-project.org/package=lme4) fits mixed models
  in a small set of GLM families. Only the mean gets a formula, and
  each grouping factor gets one unstructured covariance.
- [glmmTMB](https://glmmtmb.github.io/glmmTMB/) adds many families,
  structured random-effect covariances, and separate formulas for
  dispersion and for zero inflation. The predictor grammar is not
  brms's: there are no nonlinear formulas, no multivariate responses
  with residual correlation, no monotonic or measurement-error
  terms, and no correlation of effects across formulas.
- [gamlss](https://cran.r-project.org/package=gamlss), and its
  successor gamlss2, give every distributional parameter its own
  additive predictor, which is the part lme4 and glmmTMB lack. Random
  effects are additive terms there rather than an lme4 grammar with
  structured covariances, and the spelling is not brms's.

frmtmb supplies the combination: distributional regression on every
parameter of the family, nonlinear formulas, multivariate responses
with `rescor`, effects correlated across formulas with `|ID|`,
custom families written as plain R log-densities, and the lme4
random-effect grammar with structured and spatial covariances, all
under maximum likelihood or REML. The package generates each model's
objective as an R closure and differentiates it on the RTMB tape, so
it does not select from a fixed set of likelihoods. Features that
would each need a bespoke implementation become ordinary code paths,
and they compose with each other.

The second need is migration. brms code ports by removing the priors
and changing `brm()` to `frm()`. A measured audit of the brms
vignettes puts about 7 of 10 of their model calls through that
transform unchanged. A user can therefore screen models at
millisecond speed, then return to brms for the final Bayesian fit
with the same formula. `vignette("brms-migration")` maps the
features and states, under "When you still want brms", the cases
that belong in brms.

**Audience.** Applied statisticians and quantitative researchers who
already write brms formulas, and who need a likelihood-based answer:
for a maximum-likelihood workflow, for a model that must fit in
seconds inside a simulation or a bootstrap, for a report that asks
for confidence intervals and likelihood-ratio tests, or for a
teaching setting where no toolchain can be installed. Fields where
this comes up include psychology and psychophysics, ecology,
pharmacometrics, and meta-analysis.

Related work: glmmTMB is the closest relative, a mature and fast
TMB-based mixed-model package. frmtmb matches its fits where the
models overlap and follows its conventions in several places.
[BayesRTMB](https://github.com/norimune/BayesRTMB) is a
Bayesian-first, Stan-like modeling layer on the same RTMB backend.
[brms](https://paulbuerkner.com/brms/) itself is the right tool when
you want priors and full posterior inference; the shared grammar is
meant to make movement between the two easy.

## Installation

frmtmb is not on CRAN yet. Install the development version from
GitHub:

```r
# install.packages("remotes")
remotes::install_github("aforren1/frmtmb")
```

One optional dependency, RTMBode, is not on CRAN either. It is
needed only by `frm_ode()`. Install it from r-universe:

```r
install.packages("RTMBode", repos = c(
  "https://kaskr.r-universe.dev",
  "https://cloud.r-project.org"))
```

## Example

```r
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

Pre-release (v0.39). The goal is a CRAN release. Validation has
three layers:

- The suite contains more than 5000 tests. Every model class is compared
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

- `bf()` uses the brms spelling, and the family is a separate
  argument: `frm(bf(y ~ x), family = poisson(), data = d)`. With no
  family, `frm()` fits gaussian, as brms and glmmTMB do. See
  "Statement of need" above for the port from `brm()`.
- Every distributional parameter can take its own formula with the
  full predictor grammar
  (`bf(y ~ s(x) + (1 | g), sigma ~ s(z) + (1 | g))`), or a constant.
- Random effects use the lme4 syntax and add structured covariances:
  `us`, `diag`, `homdiag`, `cs`, `homcs`, `ar1`, `hetar1`, `toep`,
  `homtoep`, `ou`, spatial `exp`/`gau`/`mat` over `num_factor(x, y)`
  coordinates, reduced-rank `rr(d =)`, and known structure
  `gr(g, cov = A)` / `gr(g, prec = Q)` / `equalto()`. The `|ID|`
  syntax correlates effects across formulas, including known-matrix
  blocks (a multi-trait animal model is one Kronecker block).
  Multi-membership terms use `mm(g1, g2)` with weights and
  member-specific `mmc()` covariates.
- Within-group residual correlation uses brms's R-side terms:
  `ar()`, `ma()`, `arma()`, `cosy()`, and `unstr()`, for gaussian
  and student responses, validated against `nlme::gls` under ML and
  REML.
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
  `categorical()` fits nominal responses, `von_mises()` circular
  ones, and `cox()` proportional hazards with a flexible baseline
  and Laplace frailties. `hmm(K)` fits hidden Markov models with
  covariate-dependent transitions and forward-backward decoding
  (`hmm_probs()`, `hmm_viterbi()`); `lca(K)` fits poLCA-style
  latent class analysis with class-membership regression.
  `custom_family()` takes a plain R log-density. The test suite fits
  a Wiener drift-diffusion model in about 15 lines.
- Ordinary differential equations: `frm_ode()` solves compartment
  models inside nonlinear formulas (population pharmacokinetics),
  with repeated dosing, infusions, and estimated bioavailability.
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
