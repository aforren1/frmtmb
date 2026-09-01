# Residuals from a frmtmb fit

`"osa"` gives one-step-ahead (conditional quantile) residuals via
[`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html):
standard-normal under a correctly specified model, valid under
correlated observations where pearson residuals mislead.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
residuals(
  object,
  type = c("response", "pearson", "deviance", "osa"),
  osa_method = NULL,
  ...
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- type:

  `"response"`, `"pearson"`, `"deviance"`, or `"osa"`.

- osa_method:

  Method for
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html);
  defaults to `"fullGaussian"` for gaussian models and
  `"oneStepGeneric"` otherwise. A truncated, censored or ordinal
  response always uses `"oneStepGeneric"` (a truncated gaussian is not
  gaussian) with the integration domain and discrete support taken from
  the [`trunc()`](https://rdrr.io/r/base/Round.html) bounds or the
  censoring window, which must then be the same for every row.

- ...:

  For `type = "osa"`: passed to
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html).

## Value

A numeric vector, `NA` on censored rows.

## Details

On a [`trunc()`](https://rdrr.io/r/base/Round.html)ed response,
`"response"` residuals are taken against the truncated mean
`E[Y | lb <= Y <= ub]`. `"pearson"` divides by the untruncated family
variance, so it is conservative there. `"osa"` builds its conditional
CDF on `[lb, ub]` (see `osa_method`).

On a `cens()`ed response, `"osa"` returns `NA` for every censored row:
what is observed there is an event (`Y > c`), not a value, and an event
has no one-step CDF. The uncensored rows get residuals conditional on
the censoring events, which needs one censoring point per side (type-I
censoring); row-varying censoring times and interval censoring are
refused.
[`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
is not a substitute on a censored fit, because
[`simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.md)
draws the latent uncensored response by default (as brms's
[`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
does) and those draws are not comparable with the observed censored
values; `simulate(censored = TRUE)` makes them comparable, but the
resulting point mass at each censoring point is not a distribution
DHARMa's rank transform can use.

Ordinal responses use `"oneStepGeneric"` over the discrete support
`1..K`, which makes the residuals randomized quantile residuals.

## Deviance residuals

`"deviance"` returns `sign(y - E[Y]) * sqrt(w * d)`, where the unit
deviance
`d = 2 * (loglik of the saturated fit - loglik at the fitted value)` is
taken with the dispersion parameter held at its estimate, and `w` is the
[`weights()`](https://rdrr.io/r/stats/weights.html) addition term (1 by
default). For the exponential-dispersion families this is the glm unit
deviance, so a fixed-effect fit reproduces
`residuals(glm(...), type = "deviance")` exactly.

Supported families: `gaussian`, `poisson`, `binomial`, `bernoulli`,
`Gamma`, `exponential`, `inverse.gaussian`, `negbinomial` (`nbinom2`),
`nbinom1`, `geometric`, `beta`, and `tweedie`. Every other family is
refused: ordinal, mixture, multinomial, hurdle, zero-inflated and
location-shift families have no standard unit deviance. `nbinom1`
follows glmmTMB and evaluates the negative-binomial size at the fitted
row's `mu / phi`; letting the size follow the saturated mean is not a
deviance (the difference goes negative).
[`trunc()`](https://rdrr.io/r/base/Round.html)ed and `cens()`ed
responses are refused as well, because the fitted likelihood there is
not the family's own density.

In a mixed model the residuals are conditional on the random-effect
modes, the glmmTMB convention: `E[Y]` is
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), not the
population-level mean.

`deviance(fit)` is unrelated: it stays `-2 * logLik(fit)` (the lme4
convention), which for a mixed model is the Laplace-approximated
marginal deviance and does **not** equal
`sum(residuals(fit, type = "deviance")^2)`.
