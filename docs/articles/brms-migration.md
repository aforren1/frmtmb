# frmtmb for brms users

frmtmb reimplements a documented subset of the brms grammar with
identical spelling, so brms model code ports mechanically: drop the
priors, change `brm()` to [`frm()`](../reference/frm.md), and get a
maximum-likelihood fit in milliseconds-to-seconds instead of an MCMC
run. Estimation is ML (or REML) with the Laplace approximation for
latent effects; inference is frequentist (Wald and profile intervals,
likelihood-ratio tests, AIC).

## What ports directly

| brms | frmtmb | notes |
|----|----|----|
| `bf(y ~ x + (x \| g))` | same | plus `us()`, [`diag()`](https://rdrr.io/r/base/diag.html), `homdiag()`, `cs()`, `ar1()` wrappers |
| `bf(y ~ x, sigma ~ z + (1 \| g))` | same | any dpar, full grammar |
| `bf(y ~ x, sigma = 1)` | same | fixed via TMB `map` |
| `s(x)`, `t2(x, z)` | same | mgcv smooths in any dpar; `te()`/`ti()` unsupported |
| `(1 \| p \| g)` | same | cross-formula RE correlation |
| `mvbind(y1, y2) ~ x`, [`mvbf()`](../reference/mvbf.md), [`set_rescor()`](../reference/mvbf.md) | same | per-response families; `rescor` gaussian-only |
| `y \| trials(n)`, `weights(w)`, `cens(c)`, `trunc(lb=, ub=)` | same | `cens`/`trunc` need a CDF-carrying family |
| `nl = TRUE` | same | provide `start`; se.fit on the nonlinear mu not yet |
| [`custom_family()`](../reference/frmtmb_family.md) | same idea | the lpdf is plain R over RTMB advectors, not Stan code |
| [`cumulative()`](../reference/frmtmb-families.md), [`multinomial()`](../reference/frmtmb-families.md) | same | multinomial takes `K` explicitly |

## What changes

- `brm(...)` becomes `frm(...)`; there are no priors, chains, or warmup.
  `REML = TRUE` replaces priors as the small-sample correction for
  variance components.
- Posterior summaries become ML estimates:
  [`fixef()`](../reference/fixef.md) returns point estimates,
  [`confint()`](https://rdrr.io/r/stats/confint.html) gives Wald or
  profile intervals, [`anova()`](https://rdrr.io/r/stats/anova.html)
  gives likelihood-ratio tests, and
  [`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
  replace `loo()`.
- [`posterior_predict()`](../reference/posterior_epred.md) becomes
  [`simulate()`](https://rdrr.io/r/stats/simulate.html);
  [`posterior_epred()`](../reference/posterior_epred.md) over `newdata`
  becomes `predict(newdata, type = "response")`.
- `mo()` monotonic effects, `mi()` one-step imputation of continuous
  predictors, `mi(sdx)` measurement error (the `me()` replacement, as in
  current brms), `cs()` category-specific ordinal effects, `gp(x)` /
  `gp(x, k =)` Gaussian processes, and
  [`mixture()`](../reference/mixture.md) families (including group-level
  latent classes via `groups = ~g`) all work. `mo()`/`mi()` support
  two-way `:`/`*` interactions with numeric terms; discrete predictors
  cannot be imputed in-model, the same restriction Stan has. Multiply
  imputed data can instead go through
  [`frm_multiple()`](../reference/frm_multiple.md) (Rubin’s rules), the
  `brm_multiple()` analog. Mixture fits are ML: expect multimodality,
  compare starts ([`frm_allfit()`](../reference/frm_allfit.md)).
- Not supported: `ar()/ma()` residual autocorrelation terms (use the
  `ar1()` random-effect structure).
- glmer’s proportion-response idiom (`weights = size`) becomes
  `y | trials(size)` with either proportions or counts;
  [`weights()`](https://rdrr.io/r/stats/weights.html) in the formula
  stays a frequency weight, as in brms.
- `se()` works as in brms (meta-analysis), including
  `se(x, sigma = TRUE)` and the phylogenetic version with
  `gr(g, cov = A)`.

## Method conventions

frmtmb keeps the brms spelling for functions that originate in brms and
the stats/lme4 spelling for standard generics:

- brms-origin, same names and argument spellings:
  [`conditional_effects()`](../reference/conditional_effects.md)
  (includes what `conditional_smooths()` covers),
  [`hypothesis()`](../reference/hypothesis.md) (a delta-method Wald test
  here), `pp_check()` (simulate-based, through bayesplot’s `ppc_*`
  functions), [`prior_summary()`](../reference/prior_summary.md),
  [`fixef()`](../reference/fixef.md),
  [`ranef()`](../reference/ranef.md),
  [`VarCorr()`](../reference/VarCorr.md),
  [`ngrps()`](../reference/ngrps.md),
  [`variables()`](../reference/variables.md) (the usable parameter
  names, e.g. `sd_Subject__Days`).
- stats-origin generics follow stats/lme4/glmmTMB, not brms:
  [`predict()`](https://rdrr.io/r/stats/predict.html) returns a vector
  (with `se.fit = TRUE`, a list), not a draws matrix with
  `Estimate`/`Q2.5` columns;
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) return
  vectors; [`confint()`](https://rdrr.io/r/stats/confint.html) takes
  `level =`, not `probs =`; `re.form` replaces `re_formula`.
- [`coef()`](https://rdrr.io/r/stats/coef.html) is the per-group
  convention shared by all three packages: fixed effects plus
  conditional modes per grouping level.
- `plot(fit)` draws residual diagnostics, not MCMC traces; for
  simulation-based residual checks use
  [`dharma_residuals()`](../reference/dharma_residuals.md).

## When you still want brms

Genuine prior information, full posterior uncertainty for derived
quantities, models with discrete latent structure beyond
observation-level mixtures, or `loo`-based model comparison. A practical
workflow is model screening in frmtmb and a final fit in brms;
`as_tmbstan(fit)` also runs NUTS on the frmtmb objective directly.
