# frmtmb

frmtmb fits regression models specified with a brms-style formula
grammar, by maximum likelihood with the Laplace approximation for
latent effects. The backend is
[RTMB](https://cran.r-project.org/package=RTMB): each model's
objective is generated as an R closure and differentiated on the TMB
AD tape. No MCMC, no Stan, and no compilation at run time; refits
re-tape in milliseconds, which the bootstrap, influence, and
multi-start machinery exploit.

Documentation: <https://aforren1.github.io/frmtmb/>

## Status

Pre-release (v0.28), working toward CRAN. Validation is layered:

- about 2500 tests; every model class is checked against an exact
  reference (glmmTMB, lme4, mgcv, MASS, survival, nnet, GLMMadaptive,
  quantreg, mice, closed-form marginals, or hand-written ML);
- the model-building layer is cross-validated against brms itself:
  design matrices, random-effect structures, and special-term data
  agree with `brms::make_standata()` to near machine precision, and
  an opt-in tier verifies our estimates equal the mode of brms's own
  generated Stan programs;
- a pairwise grammar fuzzer sweeps feature combinations against
  metamorphic invariants, and the resulting compatibility map is
  queryable (`frm_compat()`) and published as the
  [feature compatibility](https://aforren1.github.io/frmtmb/articles/compatibility.html)
  article.

## Model grammar

- `bf()` with the brms spelling; families attach with `+`. Ports of
  brms code drop the priors and change `brm()` to `frm()`.
- Distributional regression: every distributional parameter takes its
  own formula with the full predictor grammar
  (`bf(y ~ s(x) + (1 | g), sigma ~ s(z) + (1 | g))`), or a constant.
- Random effects: lme4 syntax plus structured covariances - `us`,
  `diag`, `homdiag`, `cs`, `homcs`, `ar1`, `hetar1`, `toep`,
  `homtoep`, `ou`, spatial `exp`/`gau`/`mat` over `num_factor(x, y)`
  coordinates, reduced-rank `rr(d =)`, known structure
  `gr(g, cov = A)` / `gr(g, prec = Q)` / `equalto()`, and `|ID|`
  correlation across formulas.
- Spatial Gaussian Markov random fields: `car(M, gr = g, type =)` with
  brms's spelling and all four of its types (`escar`, `esicar`,
  `icar`, `bym2`, the last with brms's scaling convention), and
  `spde(fem, gr = node)` for a Matern field over a finite-element mesh
  (`fmesher`/INLA matrices as fixed data). Precisions are assembled on
  the tape from fixed sparse matrices and every normalizing constant
  is analytic.
- mgcv smooths `s()`/`t2()` in any linear predictor, including
  matrix-covariate terms: scalar-on-function, function-on-scalar, and
  function-on-function regression.
- Gaussian processes: `gp(x)` exact with kriging prediction,
  `gp(x, k = 30)` Hilbert-space (brms's exact input convention, so
  the same call is the same approximation), up to three dimensions
  (`gp(x1, x2)`, per-dimension lengthscales, `iso = TRUE` to share).
- Special terms: `mo()` monotonic effects, `mi()` one-step imputation
  of continuous predictors, `mi(sdx)` measurement error, `cs()`
  category-specific ordinal effects; `mo()`/`mi()` support two-way
  interactions.
- Multivariate models (`mvbf`, per-response families, `rescor`),
  nonlinear formulas (`nl = TRUE`), matrix-response `multinomial(K)`.
- Families: the usual GLM(M) set plus student, tweedie, compois,
  beta-binomial, skew-normal, ex-gaussian, weibull, shifted
  lognormal, quantile regression (`asym_laplace`), zero-inflated and
  hurdle variants, four ordinal families, and finite `mixture()`
  families - including group-level latent classes
  (`mixture(..., groups = ~g)`) with class-specific random effects
  (growth-mixture models) and multivariate gaussian components
  (`mixture_mvn(K, D, model =)`, model-based clustering over mclust's
  covariance taxonomy, with covariate-dependent means and gating).
  `custom_family()` takes a
  plain R
  log-density (the test suite fits a Wiener drift-diffusion model in
  about 15 lines).
- Addition terms: `weights()`, `trials()` (counts or proportions),
  `cens()`, `trunc()`, `se()` (meta-analysis), `mi()`,
  `vint()`/`vreal()`.

## Estimation and inference

- ML and REML; adaptive quadrature (`quadrature = TRUE`) for scalar
  random effects; MAP via brms-style `set_prior()`; hard bounds;
  pluggable optimizers; `frmtmb_control(profile = TRUE)` for
  many-coefficient models, `sparse_x = TRUE` for many-level fixed
  factors, `autoscale = TRUE` for badly scaled predictors,
  `verbose =` for timed stage progress on slow fits.
- `confint()` (Wald/profile/likelihood-root/bootstrap),
  `hypothesis()` (Wald/profile/bootstrap, with natural-scale
  `sd_`/`cor_` names for ICC-type quantities), `anova()` LRTs and
  `drop1()`, `frm_bootstrap()`, `influence()` with
  `cooks.distance()`/`dfbetas()`, `frm_allfit()`, `frm_multiple()`
  (Rubin pooling, including variance components and hypotheses).
- Simulation both ways: `simulate()` from a fit, and
  `frm_simulate()` from a bare design with natural-scale parameters
  (`sd_g__Intercept = 0.7`) or with parameters drawn from
  `set_prior()` specifications per replicate - the
  `sample_prior = "only"` prior-predictive workflow, without MCMC.
- Post-hoc NUTS on the fitted objective: `frm_sample()` with a full
  draws method surface (`posterior_epred()`, `posterior_predict()`,
  `hypothesis()`, `pp_check()`) and `check_laplace()` to audit the
  approximation.
- Diagnostics: one-step-ahead residuals calibrated for censored,
  truncated, and ordinal responses; deviance residuals across the
  GLM families; response-scale `se.fit` for every family through the
  joint delta method; DHARMa, `pp_check()`, `plot()`,
  `conditional_effects()`; ecosystem hooks for emmeans,
  marginaleffects, and insight/easystats.

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

## Why

No frequentist equivalent of brms exists. glmmTMB is one fixed C++
likelihood behind a formula front end; frmtmb inverts that design and
compiles the formula into the objective. That turns features that are
structural dead ends elsewhere - random effects in any distributional
parameter, nonlinear predictors, per-response families, monotonic
effects, in-model imputation, latent-class mixtures, custom families
without C++ or Stan code - into ordinary code paths.

The architecture also sidesteps long-standing open issues in the
reference packages: zero prior weights are exactly equivalent to
subsetting, nonlinear-model standard errors match nlme where
`nlmer`'s are orders of magnitude off, crossed grouping `(1 | a*b)`
and call-valued grouping factors `(1 | factor(x))` work, tensor
smooths match mgcv where glmmTMB has none, REML predictions agree
with `fixef()`, and lme4's penalized-IRLS pathologies are
structurally absent (the audit lives in `dev/test-backlog.md`).

Related work: [glmmTMB](https://glmmtmb.github.io/glmmTMB/) (fixed
likelihood, C++ TMB; frmtmb matches its fits where the models
overlap) and [BayesRTMB](https://github.com/norimune/BayesRTMB) (a
Bayesian-first, Stan-like modeling layer on the same RTMB backend).

## Life cycle

frmtmb is maturing. The package is not yet on CRAN.

- **The model grammar is stable.** It follows brms, so it changes only
  when brms changes. Formulas that fit today will fit in later
  versions.
- **The fitted-object API is stable.** `frm()`, the accessor methods
  (`coef()`, `confint()`, `vcov()`, `predict()`, and the rest), and
  the family constructors keep their current behavior.
- **Internal structure can change.** Fields of the fitted object that
  no exported method reaches are not part of the API. Use the
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

Run `citation("frmtmb")` in R. Please also cite RTMB and TMB
(Kristensen et al. 2016, <doi:10.18637/jss.v070.i05>).

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Everyone taking part must
follow the [Code of Conduct](CODE_OF_CONDUCT.md).

See [SPEC.md](SPEC.md) for the design and `NEWS.md` for the
changelog.
