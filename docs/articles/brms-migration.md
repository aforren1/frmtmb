# frmtmb for brms users

frmtmb reimplements a documented subset of the brms grammar with
identical spelling, so brms model code ports mechanically: drop the
priors, change `brm()` to
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md), and get a
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
| `mvbind(y1, y2) ~ x`, [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md), [`set_rescor()`](https://aforren1.github.io/frmtmb/reference/mvbf.md) | same | per-response families; `rescor` gaussian-only |
| `y \| trials(n)`, `weights(w)`, `cens(c)`, `trunc(lb=, ub=)` | same | `cens`/`trunc` need a CDF-carrying family |
| `nl = TRUE` | same | provide `start`; se.fit on the nonlinear mu not yet |
| [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md) | same idea | the lpdf is plain R over RTMB advectors, not Stan code |
| [`cumulative()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md), [`multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md) | same | multinomial takes `K` explicitly |
| `car(M, gr = g, type =)` | same | all four types; `M` is looked up in the data or the calling environment, not in `data2` |

## What changes

- `brm(...)` becomes `frm(...)`; there are no priors, chains, or warmup.
  `REML = TRUE` replaces priors as the small-sample correction for
  variance components.
- Posterior summaries become ML estimates:
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md)
  returns point estimates,
  [`confint()`](https://rdrr.io/r/stats/confint.html) gives Wald or
  profile intervals, [`anova()`](https://rdrr.io/r/stats/anova.html)
  gives likelihood-ratio tests, and
  [`AIC()`](https://rdrr.io/r/stats/AIC.html)/[`BIC()`](https://rdrr.io/r/stats/AIC.html)
  replace `loo()`.
- [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  becomes [`simulate()`](https://rdrr.io/r/stats/simulate.html);
  [`posterior_epred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  over `newdata` becomes `predict(newdata, type = "response")`. On a
  `cens()` response both draw the LATENT uncensored value, as brms does;
  `simulate(censored = TRUE)` applies the censoring mechanism instead.
- `sample_prior = "only"` plus
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
  becomes `frm_simulate(formula, data, priors = set_prior(...))`, which
  draws a parameter vector per simulation and returns it alongside the
  responses.
- `mo()` monotonic effects, `mi()` one-step imputation of continuous
  predictors, `mi(sdx)` measurement error (the `me()` replacement, as in
  current brms), `cs()` category-specific ordinal effects, `gp(x)` /
  `gp(x, k =)` Gaussian processes, and
  [`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
  families (including group-level latent classes via `groups = ~g`) all
  work. `mo()`/`mi()` support two-way `:`/`*` interactions with numeric
  terms; discrete predictors cannot be imputed in-model, the same
  restriction Stan has. Multiply imputed data can instead go through
  [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  (Rubin’s rules), the `brm_multiple()` analog. Being frequentist, that
  fit also needs pooled tests, which brms gets for free from the
  posterior:
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  pools any function of the parameters, and
  [`anova()`](https://rdrr.io/r/stats/anova.html) on the result compares
  two nested
  [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md)
  fits with the D1, D2 and D3 rules of `mice` (D3 by default). Mixture
  fits are ML: expect multimodality, compare starts
  ([`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)).
- [`binomial()`](https://rdrr.io/r/stats/family.html) without `trials()`
  is accepted and means Bernoulli, the
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html) convention; brms
  rejects it and asks for `trials()` or
  [`bernoulli()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md).
  The divergence is permissive, so brms code ports unchanged, but frmtmb
  code written this way does not port back.
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
  [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
  (includes what `conditional_smooths()` covers),
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  (a delta-method Wald test here),
  [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md)
  (simulate-based, through bayesplot’s `ppc_*` functions),
  [`prior_summary()`](https://aforren1.github.io/frmtmb/reference/prior_summary.md),
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
  [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md),
  [`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md),
  [`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.md),
  [`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
  (the usable parameter names, e.g. `sd_Subject__Days`).
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
  [`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md).

## Seeing the model code

brms users inspect the generated model with `stancode(fit)`. The frmtmb
analog is the objective function itself: an ordinary R closure that
`RTMB::MakeADFun()` tapes, rebuilt deterministically from the fitted
frame.

``` r

nll <- frmtmb:::build_objective(fit$frame)
nll(fit$estimates)   # joint negative log density at the estimates
```

Two properties matter when reading it:

- The closure takes the named parameter list (the shape of
  `fit$frame$par_template`) and returns the JOINT negative log density:
  random-effect prior terms are inside it. The marginal (Laplace)
  likelihood the optimizer minimized is what `MakeADFun(random = )`
  makes of it, exposed as `fit$obj$fn` over the flat outer-parameter
  vector; at the optimum `fit$obj$fn(fit$opt$par)` equals
  `-logLik(fit)`.
- All data is baked into the closure’s environment (the reason refits
  re-tape instead of recompiling), so `environment(nll)` holds the
  assembled frame, and the closure evaluates on plain numeric parameter
  lists without any tape.

`RTMB::MakeADFun(nll, fit$frame$par_template, random = "b", map = fit$frame$map)`
reproduces `fit$obj`. There is no generated code to print: the closure
IS the model, and [`body()`](https://rdrr.io/r/base/body.html),
[`environment()`](https://rdrr.io/r/base/environment.html), and the
debugger work on it like on any R function.

## When you still want brms

Genuine prior information, full posterior uncertainty for derived
quantities, models with discrete latent structure beyond
observation-level mixtures, or `loo`-based model comparison. A practical
workflow is model screening in frmtmb and a final fit in brms;
`as_tmbstan(fit)` also runs NUTS on the frmtmb objective directly.
