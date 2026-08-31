# frmtmb

frmtmb fits regression models specified with a brms-style formula grammar,
by maximum likelihood with the Laplace approximation for random effects.
The backend is [RTMB](https://cran.r-project.org/package=RTMB): model
objectives are generated as R closures and differentiated on the TMB AD
tape. No MCMC, no Stan, and no compilation at run time.

## Status

Early development (v0.6). Current scope:

- Diagnostics: `diagnose()` convergence forensics, simulation-based
  residuals via DHARMa (`dharma_residuals()`), and full-Bayes checks by
  handing the fitted objective to NUTS (`as_tmbstan()`; bayesplot works
  on the draws)
- Robustness hardening mined from lme4/glmmTMB/brms issue history:
  frozen data-dependent bases (poly/ns/scale), rank-deficient designs,
  non-default contrasts, `na.exclude`, nested `(1 | a/b)` grouping,
  fits that outlive their calling environment
- Deferred `sdreport`: standard errors are computed on first use, which
  cuts about a quarter off fit time in fit-and-predict or bootstrap
  loops (`se = TRUE` restores eager behavior)

- Nonlinear formulas (`nl = TRUE`): arbitrary R expressions over named
  parameters, each with the full predictor grammar
  (`bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE)`)
- Custom families as plain R functions (`custom_family()`), with an AD
  safety checker (`check_custom_family()`)
- Ordinal regression: `cumulative(logit/probit)` with ordered
  thresholds, ordered-factor responses, random effects
- Censoring and truncation: `y | cens(c) ~ ...`, `y | trunc(lb = 0)`
  (gaussian, lognormal)
- emmeans support and a tmbstan bridge (`as_tmbstan()`) for MCMC on the
  same objective
- Function-on-scalar (functional response) regression via the
  long-format + smooth representation, validated against mgcv

- Multivariate models: `mvbf()`, `bf() + bf()`, `mvbind()`; a family per
  response; residual correlation (`rescor = TRUE`, gaussian); and brms
  `|ID|` syntax for random-effect correlation across responses
  (`(1 | p | g)` in several formulas) - open glmmTMB issue #1267
- Matrix responses: `multinomial(K)` takes an n x K count matrix with
  per-category linear predictors
- Zero-inflation and hurdle: `zero_inflated_poisson()`,
  `zero_inflated_negbinomial()`, `hurdle_poisson()`, with full `zi ~` /
  `hu ~` formulas

- Distributional regression: every distributional parameter takes its own
  formula with the full predictor grammar, including random effects and
  smooths (`bf(y ~ s(x) + (1 | g), sigma ~ s(z) + (1 | g))`), or a fixed
  constant (`bf(y ~ x, sigma = 1)`)
- mgcv smooths `s()` and `t2()` in any linear predictor; smoothing
  parameters are estimated as variance components by marginal ML/REML
  (matches `mgcv::gam(method = "ML")`)
- Families: gaussian, poisson, binomial (with `trials()`), Gamma,
  lognormal, student, negbinomial, nbinom1, beta, tweedie, compois
- `lme4`-style random effects: `(1 | g)`, `(1 + x | g)`, `(x || g)`,
  with `us()`, `diag()`, `homdiag()`, `cs()`, and `ar1()` covariance
  structures
- `weights()` addition term and `offset()`
- ML and REML
- Methods: `summary`, `logLik`, `AIC`, `vcov`, `fixef`, `ranef`,
  `VarCorr`, `predict` (newdata, `se.fit`, `re.form`), `fitted`,
  `residuals`, `simulate`, `confint` (Wald, profile, uniroot),
  `anova` (LRT), `update`, `diagnose`

See [SPEC.md](SPEC.md) for the full design and roadmap: ordinal families
and `cens()`/`trunc()` (v0.4.x), nonlinear formulas and custom families
as plain R functions (v0.5).

## Example

```r
library(frmtmb)
data(sleepstudy, package = "lme4")

fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
           data = sleepstudy)
summary(fit)

# distributional regression: model the residual SD, with group effects
fit2 <- frm(bf(Reaction ~ Days + (Days | Subject),
               sigma ~ Days + (1 | Subject)) + gaussian(),
            data = sleepstudy)
```

## Why

No frequentist equivalent of brms exists. glmmTMB is one fixed C++
likelihood behind a formula front end; frmtmb inverts that design and
compiles the formula into the objective. That makes features that are
structural dead ends in glmmTMB (random effects in any distributional
parameter formula, nonlinear predictors, per-response families, custom
families without C++) into ordinary code paths.

Fits are validated against glmmTMB, lme4, mgcv, and gamlss in the test
suite.
