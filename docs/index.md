# frmtmb

[![R-CMD-check](https://github.com/aforren1/frmtmb/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/aforren1/frmtmb/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/test-coverage.yaml)
[![pkgcheck](https://github.com/aforren1/frmtmb/actions/workflows/pkgcheck.yaml/badge.svg)](https://github.com/aforren1/frmtmb/actions/workflows/pkgcheck.yaml)
[![Project Status:
Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

frmtmb fits regression models that you specify with a brms-style formula
grammar. Estimation is maximum likelihood, with the Laplace
approximation for latent effects. The backend is
[RTMB](https://cran.r-project.org/package=RTMB). The package generates
each model’s objective as an R closure and differentiates it on the TMB
AD tape. There is no MCMC, no Stan, and no compilation at run time. A
refit re-tapes in milliseconds. The bootstrap, influence, and
multi-start machinery build on that speed.

Documentation: <https://aforren1.github.io/frmtmb/>

## Statement of need

brms gave R a formula grammar that covers a very wide model space. It
reaches that space only through Bayesian estimation in Stan. A user who
wants the same grammar under maximum likelihood must divide one model
across several packages, accept a smaller model, or wait for a sampler.
frmtmb closes that gap. It fits the brms grammar by maximum likelihood
with the Laplace approximation, and it compiles nothing at run time.

The gap is real because the existing frequentist packages each stop at a
different place.

- [lme4](https://cran.r-project.org/package=lme4) fits mixed models in a
  small set of GLM families. Only the mean gets a formula, and each
  grouping factor gets one unstructured covariance.
- [glmmTMB](https://glmmtmb.github.io/glmmTMB/) adds many families,
  structured random-effect covariances, and separate formulas for
  dispersion and for zero inflation. The predictor grammar is not
  brms’s: there are no nonlinear formulas, no multivariate responses
  with residual correlation, no monotonic or measurement-error terms,
  and no correlation of effects across formulas.
- [gamlss](https://cran.r-project.org/package=gamlss), and its successor
  gamlss2, give every distributional parameter its own additive
  predictor, which is the part lme4 and glmmTMB lack. Random effects are
  additive terms there rather than an lme4 grammar with structured
  covariances, and the spelling is not brms’s.

frmtmb supplies the combination: distributional regression on every
parameter of the family, nonlinear formulas, multivariate responses with
`rescor`, effects correlated across formulas with `|ID|`, custom
families written as plain R log-densities, and the lme4 random-effect
grammar with structured and spatial covariances, all under maximum
likelihood or REML. The package generates each model’s objective as an R
closure and differentiates it on the RTMB tape, so it does not select
from a fixed set of likelihoods. Features that would each need a bespoke
implementation become ordinary code paths, and they compose with each
other.

The second need is migration. brms code ports by changing `brm()` to
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md). The
priors can stay where they are: `frm(prior = )` takes brms’s own
spelling,
[`prior()`](https://aforren1.github.io/frmtmb/reference/prior.md) builds
the specification
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
does, and a prior object brms itself built is translated. The fit is
then MAP, so a prior is a penalty on the likelihood rather than a
posterior. A measured audit of the brms vignettes puts about 7 of 10 of
their model calls through that transform unchanged; over every call the
vignettes make, post-processing included, about 4 of 10 run unchanged,
and most of the rest need one spelling change or get a refusal that
names the replacement. A user can therefore screen models at millisecond
speed, then return to brms for the final Bayesian fit with the same
formula.
[`vignette("brms-migration")`](https://aforren1.github.io/frmtmb/articles/brms-migration.md)
maps the features and states, under “When you still want brms”, the
cases that belong in brms.

**A scoped core, with the rest in companion packages.** The core package
fits models and reports on them. Everything that needs another engine or
another literature ships beside it: NUTS sampling, ordinary differential
equations, discrete latent states, and drift-diffusion response times
each live in their own package in this repository. That keeps the core’s
dependencies on CRAN, keeps its promise of no compilation exact, and
makes the part that matters reviewable on its own. The seam is public:
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md),
[`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
and the registration functions are the same interface the companion
packages use, so a family written outside this repository reaches the
core the same way `hmm()` does.

**Audience.** Applied statisticians and quantitative researchers who
already write brms formulas, and who need a likelihood-based answer: for
a maximum-likelihood workflow, for a model that must fit in seconds
inside a simulation or a bootstrap, for a report that asks for
confidence intervals and likelihood-ratio tests, or for a teaching
setting where no toolchain can be installed. Fields where this comes up
include psychology and psychophysics, ecology, pharmacometrics, and
meta-analysis.

Related work: glmmTMB is the closest relative, a mature and fast
TMB-based mixed-model package. frmtmb matches its fits where the models
overlap and follows its conventions in several places.
[BayesRTMB](https://github.com/norimune/BayesRTMB) is a Bayesian-first,
Stan-like modeling layer on the same RTMB backend.
[qbrms](https://github.com/Tony-Myers/qbrms) also reads brms syntax and
sends it to a different engine, INLA, which fits approximate posteriors
for latent Gaussian models; frmtmb instead maximizes the likelihood and
takes families and nonlinear bodies that INLA’s model class does not
cover. [brms](https://paulbuerkner.com/brms/) itself is the right tool
when you want priors and full posterior inference; the shared grammar is
meant to make movement between the two easy.

## Installation

frmtmb is not on CRAN yet. Install the development version from GitHub:

``` r

# install.packages("remotes")
remotes::install_github("aforren1/frmtmb")
```

Five companion packages live in this repository, each installed the same
way and each documented on its own part of the site.

| package | adds |
|----|----|
| [frmtmb.sample](https://aforren1.github.io/frmtmb/frmtmb.sample/) | NUTS sampling through tmbstan, and the posterior method surface |
| [frmtmb.latent](https://aforren1.github.io/frmtmb/frmtmb.latent/) | discrete latent states: `hmm()` and `lca()` |
| [frmtmb.ode](https://aforren1.github.io/frmtmb/frmtmb.ode/) | ordinary differential equation dynamics: `frm_ode()` |
| [frmtmb.ddm](https://aforren1.github.io/frmtmb/frmtmb.ddm/) | response-time models: `wiener()`, `gddm()`, and the racing `lba()` |
| [frmtmb.spline](https://aforren1.github.io/frmtmb/frmtmb.spline/) | curve inference on any fitted smooth, and the flexible parametric survival family `royston_parmar()` |

``` r

remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.sample")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.latent")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.ode")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.ddm")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.spline")
```

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

[`vignette("habit")`](https://aforren1.github.io/frmtmb/articles/habit.md)
is a longer worked example: it replicates a published
response-preparation model, fits every participant in a few seconds, and
then fits the hierarchical version that the original per-participant
procedure could not express.

## Status

Pre-release. The goal is a CRAN release. Validation has three layers:

- Every model class is compared with an exact external reference:
  glmmTMB, lme4, mgcv, nlme, MASS, survival, nnet, GLMMadaptive,
  quantreg, mice, closed-form marginals, or hand-written maximum
  likelihood. The core suite and each companion package’s suite run on
  every change, and each companion package has its own check workflow.
- The model-building layer is compared with brms itself. Design
  matrices, random-effect structures, and special-term data agree with
  [`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html)
  to near machine precision. An opt-in tier verifies that our
  log-likelihood equals the Stan program’s log density at the estimate.
- A pairwise grammar fuzzer sweeps feature combinations against
  metamorphic invariants. The resulting compatibility map is queryable
  with
  [`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md),
  which the companion packages contribute their own rows to, and is
  published as the [feature
  compatibility](https://aforren1.github.io/frmtmb/articles/compatibility.html)
  article.

## Model grammar

- [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) uses the
  brms spelling, and the family is a separate argument:
  `frm(bf(y ~ x), family = poisson(), data = d)`. With no family,
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) fits
  gaussian, as brms and glmmTMB do. See “Statement of need” above for
  the port from `brm()`.
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
  Nonlinear formulas use `nl = TRUE`. Inside a nonlinear body, a bare
  [`pnorm()`](https://rdrr.io/r/stats/Normal.html) or
  [`qgamma()`](https://rdrr.io/r/stats/GammaDist.html) is the
  tape-capable version, so a process model reads as it would on paper.
  `multinomial(K)` takes a matrix response.
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
  [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  takes a plain R log-density; section 11 of
  [`vignette("case-studies")`](https://aforren1.github.io/frmtmb/articles/case-studies.md)
  writes one with a data-bounded link and checks it against the built-in
  family it shadows.
- The addition terms are
  [`weights()`](https://rdrr.io/r/stats/weights.html), `trials()`
  (counts or proportions), `cens()`,
  [`trunc()`](https://rdrr.io/r/base/Round.html), `se()`
  (meta-analysis), `mi()`, and `vint()`/`vreal()`. A companion package
  can register one of its own.
- The companion packages extend the same grammar. `frmtmb.latent` adds
  `hmm(K)` for hidden Markov models, with covariate-dependent
  transitions and forward-backward decoding, and `lca(K)` for latent
  class analysis with class-membership regression. `frmtmb.ode` adds
  `frm_ode()`, which solves compartment models inside nonlinear formulas
  (population pharmacokinetics), with repeated dosing, infusions, and
  estimated bioavailability. `frmtmb.ddm` adds `wiener()`, the
  drift-diffusion first-passage density, with a formula for each of the
  drift rate, boundary separation, non-decision time and starting bias,
  and Ratcliff’s across-trial variability as three more; `gddm()`, the
  generalized drift-diffusion model with time-varying drift and
  collapsing bounds, solved on a fixed grid; and `lba()`, the linear
  ballistic accumulator, which races any number of accumulators and so
  reaches choices with more than two alternatives.

## Estimation and inference

- Estimation is ML or REML. `quadrature = TRUE` uses adaptive quadrature
  for scalar random effects.
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  and [`prior()`](https://aforren1.github.io/frmtmb/reference/prior.md)
  give MAP estimation with the brms spelling, `nlpar =` and `resp =`
  included, and the same call carries hard bounds through `lb`/`ub`: the
  prior vocabulary is the only way to bound a parameter, and it reaches
  every one of them, including the residual-correlation classes `ar`,
  `ma`, `cosy`, `cortime` and `rescor`.
  [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md)
  offers pluggable optimizers, `profile = TRUE` for many-coefficient
  models, `sparse_x = TRUE` for many-level fixed factors,
  `autoscale = TRUE` for badly scaled predictors, and `verbose =` for
  timed progress on slow fits.
  [`par_template()`](https://aforren1.github.io/frmtmb/reference/par_template.md)
  reports the names and starting values a model takes, before there is a
  fit to ask.
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
- Posterior inference lives in `frmtmb.sample`. Its `frm_sample()` runs
  NUTS on the fitted objective, or on a formula with brms’s default
  priors, and returns a full draws surface (`posterior_epred()`,
  `posterior_predict()`,
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md),
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md),
  [`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md)). Its
  `check_laplace()` measures the Laplace and Wald approximations against
  the sampler on the same objective.
- Diagnostics include one-step-ahead residuals calibrated for censored,
  truncated, and ordinal responses, deviance residuals across the GLM
  families, and response-scale `se.fit` for every family through the
  joint delta method. DHARMa,
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  are supported. Hooks for emmeans, marginaleffects, and
  insight/easystats are registered.

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
- **The extension interface is public and versioned.**
  [`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md),
  [`frmtmb_structure()`](https://aforren1.github.io/frmtmb/reference/frmtmb_structure.md)
  and the registration functions are documented as one contract, and the
  companion packages in this repository are written against it with no
  reach into internals. `frmtmb.ddm` was built that way as a test of it.
  Each companion package pins the core version it needs, and the pin
  moves when the interface does.
- **Internal structure can change.** Fields of the fitted object that no
  exported method reaches are not part of the API. Use the accessors.
- **Some parts are still moving.** Multivariate coverage of the post-fit
  methods and the mixture families gain features between releases.

Version numbers stay below 1.0 until the CRAN release. Breaking changes
are listed in `NEWS.md`.

The package is actively developed and maintained. Planned work: full
multivariate support in every post-fit method, more reference
comparisons, and CRAN submission.

## Citation

Run `citation("frmtmb")` in R. Please also cite RTMB and TMB, and
tmbstan if you sampled.

## Contributing

See
[CONTRIBUTING.md](https://aforren1.github.io/frmtmb/CONTRIBUTING.md).
Everyone taking part must follow the [Code of
Conduct](https://aforren1.github.io/frmtmb/CODE_OF_CONDUCT.md).

See [SPEC.md](https://aforren1.github.io/frmtmb/SPEC.md) for the design
and `NEWS.md` for the changelog.
