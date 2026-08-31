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
- `posterior_predict()` becomes
  [`simulate()`](https://rdrr.io/r/stats/simulate.html);
  `conditional_effects()` is not implemented yet (use
  `predict(newdata, se.fit = TRUE)` over a grid, or `emmeans`).
- Not supported: `mi()`, `me()`, `mo()`, `gp()`, `cs()`, `mixture()`,
  `ar()/ma()` residual autocorrelation terms (use the `ar1()`
  random-effect structure), category-specific effects.

## When you still want brms

Genuine prior information, full posterior uncertainty for derived
quantities, models outside the latent-Gaussian/Laplace class (mixtures,
discrete latents), or `loo`-based model comparison. A practical workflow
is model screening in frmtmb and a final fit in brms; `as_tmbstan(fit)`
also runs NUTS on the frmtmb objective directly.
