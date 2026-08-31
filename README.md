# frmtmb

frmtmb fits regression models specified with a brms-style formula
grammar, by maximum likelihood with the Laplace approximation for
latent effects. The backend is
[RTMB](https://cran.r-project.org/package=RTMB): each model's
objective is generated as an R closure and differentiated on the TMB
AD tape. No MCMC, no Stan, and no compilation at run time; refits
re-tape in milliseconds, which the bootstrap, influence, and
multi-start machinery exploit.

## Status

Pre-release (v0.19), working toward CRAN. The test suite holds about
850 tests; every model class is validated against an exact reference
(glmmTMB, lme4, mgcv, MASS, survival, nnet, GLMMadaptive, quantreg,
closed-form marginals, or hand-written ML) - see
`tests/testthat/` and `dev/feature-gaps.md`.

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
- mgcv smooths `s()`/`t2()` in any linear predictor, including
  matrix-covariate terms: scalar-on-function, function-on-scalar, and
  function-on-function regression.
- Gaussian processes: `gp(x)` exact, `gp(x, k = 30)` Hilbert-space.
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
  (growth-mixture models). `custom_family()` takes a plain R
  log-density (the test suite fits a Wiener drift-diffusion model in
  about 15 lines).
- Addition terms: `weights()`, `trials()` (counts or proportions),
  `cens()`, `trunc()`, `se()` (meta-analysis), `mi()`,
  `vint()`/`vreal()`.

## Estimation and inference

- ML and REML; adaptive quadrature (`quadrature = TRUE`) for scalar
  random effects; MAP via brms-style `set_prior()`; hard bounds;
  pluggable optimizers; `frmtmb_control(profile = TRUE)` for
  many-coefficient models.
- `confint()` (Wald/profile/likelihood-root), `hypothesis()`
  (Wald/profile/bootstrap, with natural-scale `sd_`/`cor_` names for
  ICC-type quantities), `anova()` LRTs, `frm_bootstrap()`,
  `influence()`, `frm_allfit()`, `frm_multiple()` (Rubin pooling),
  `frm_simulate()` (power analysis).
- Post-hoc NUTS on the fitted objective: `frm_sample()` with a full
  draws method surface (`posterior_epred()`, `posterior_predict()`,
  `hypothesis()`, `pp_check()`) and `check_laplace()` to audit the
  approximation.
- Diagnostics: OSA residuals, DHARMa, `pp_check()`, `plot()`,
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

Related work: [glmmTMB](https://glmmtmb.github.io/glmmTMB/) (fixed
likelihood, C++ TMB; frmtmb matches its fits where the models
overlap) and [BayesRTMB](https://github.com/norimune/BayesRTMB) (a
Bayesian-first, Stan-like modeling layer on the same RTMB backend).

See [SPEC.md](SPEC.md) for the design and `NEWS.md` for the
changelog.
