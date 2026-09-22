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
  type = c("response", "ordinary", "pearson", "deviance", "osa"),
  osa_method = NULL,
  ...,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975)
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- type:

  `"response"` (brms spells the same thing `"ordinary"`, and both are
  accepted), `"pearson"`, `"deviance"`, or `"osa"`. brms's
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) takes
  `newdata` in this position; this one has always taken `type` there,
  and `newdata` is refused by name.

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
  [`TMB::oneStepPredict()`](https://rdrr.io/pkg/TMB/man/oneStepPredict.html),
  and checked against that function's own formals. For every other type:
  refused, naming the argument. The `residuals.brmsfit()` arguments this
  one does not have (`newdata`, `re_formula`, `method`, `resp`,
  `ndraws`, `draw_ids`, `sort`, `summary`, `robust`, `probs`) are
  refused with the reason rather than reported as unknown names.

- ndraws, draw_ids, sort, summary, robust:

  brms's arguments. Each needs posterior draws and a maximum-likelihood
  fit has none, so each is refused by name with the reason. The default
  of each is accepted and changes nothing.

- probs:

  Probabilities of the two quantile columns, the ends of the Wald
  interval at those probabilities.

## Value

brms's summary matrix: `n` rows with `NULL` dimnames, as brms's are, and
the columns `Estimate`, `Est.Error` and one per entry of `probs`, `NA`
on censored rows. The observed response is fixed, so `Est.Error` is the
standard error of the fitted value; a `"deviance"` or `"osa"` residual
has none and reports `NA` there.

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
refused, and so is a DISCRETE family. A discrete censoring bound is
inclusive (right censoring at `k` is `Y >= k`; see
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)),
so an uncensored count's support is `[lo + 1, hi - 1]` rather than the
`[lo, hi]` this window is built on, and no reference has measured the
shifted window. Every other residual type works there.
[`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
is not a substitute on a censored fit, because
[`simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.md)
draws the latent uncensored response by default (as brms's
`posterior_predict()` does) and those draws are not comparable with the
observed censored values; `simulate(censored = TRUE)` makes them
comparable, but the resulting point mass at each censoring point is not
a distribution DHARMa's rank transform can use.

## Ordinal responses

An ordinal response has no mean, so `"response"` and `"pearson"` score
the categories by the integer codes `1..K` the likelihood itself uses:
`"response"` is `y - E[Y]` with `E[Y] = sum_k k * P(y = k)` taken from
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)'s category
probabilities, and `"pearson"` divides by the standard deviation of that
same distribution. This is the frequentist point-estimate form of what
brms's [`residuals()`](https://rdrr.io/r/stats/residuals.html) reports
on an ordinal fit (there, the observed category minus a drawn one). It
is a residual on a SCORE, not on the ordinal scale, so read it for gross
lack of fit and pattern, not as a calibrated quantity: `"osa"` and
[`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
give residuals that use only the order. `"deviance"` is refused, as it
is for every family without a standard unit deviance.

`"osa"` uses `"oneStepGeneric"` over the discrete support `1..K`, which
makes the residuals randomized quantile residuals.

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

A gaussian response with `se()` has no common dispersion for a raw
squared residual to be measured against, so the known variance enters as
a glm prior weight `sigma^2 / s_i^2` on top of `w`, where `s_i` is the
row's residual sd (the quantity `"pearson"` divides by). Without `se()`
that weight is 1 and nothing changes; `se(x)` alone maps `sigma` out at
1, leaving the familiar `1 / se_i^2` of a known-variance weighted fit.

In a mixed model the residuals are conditional on the random-effect
modes, the glmmTMB convention: `E[Y]` is
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), not the
population-level mean.

## Residual correlation terms

Under an [`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`,
`cosy()` or `unstr()` term (see
[frmtmb-autocor](https://aforren1.github.io/frmtmb/reference/frmtmb-autocor.md))
the residual covariance of a group is `D R D` with `R` unit-diagonal, so
the marginal residual SD of a row is still `sigma` and `"pearson"` is
unchanged - it divides by exactly that. What `"response"` and
`"pearson"` do NOT do is decorrelate: plotted against time within a
group they still show the fitted autocorrelation, which is the intended
reading. `"osa"` is refused, because the taped likelihood is a joint
density per group rather than a product of per-observation terms;
[`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
works, since [`simulate()`](https://rdrr.io/r/stats/simulate.html) draws
one correlated residual per group.

`deviance(fit)` is unrelated: it stays `-2 * logLik(fit)` (the lme4
convention), which for a mixed model is the Laplace-approximated
marginal deviance and does **not** equal
`sum(residuals(fit, type = "deviance")^2)`.

## See also

[frmtmb-scales](https://aforren1.github.io/frmtmb/reference/frmtmb-scales.md)
for which scale each type is on.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

# raw and variance-standardized residuals
head(residuals(fit))
#>        Estimate Est.Error       Q2.5      Q97.5
#> [1,] -0.4940281 0.1731106 -0.8333187 -0.1547376
#> [2,] -0.9276322 0.2799946 -1.4764116 -0.3788528
#> [3,] -0.7935414 0.2259217 -1.2363398 -0.3507429
#> [4,]  0.3403558 0.7525356 -1.1345869  1.8152986
#> [5,]  2.0031268 0.2886996  1.4372859  2.5689676
#> [6,]  0.4382665 0.4892102 -0.5205679  1.3971010
head(residuals(fit, type = "pearson"))
#>        Estimate Est.Error       Q2.5      Q97.5
#> [1,] -0.7028714 0.2462906 -1.1855920 -0.2201507
#> [2,] -0.9631367 0.2907112 -1.5329202 -0.3933531
#> [3,] -0.8908094 0.2536140 -1.3878837 -0.3937351
#> [4,]  0.2086995 0.4614401 -0.6957064  1.1131054
#> [5,]  2.0062658 0.2891520  1.4395383  2.5729933
#> [6,]  0.2738239 0.3056529 -0.3252448  0.8728926
# the usual overdispersion check for a poisson fit
pr <- residuals(fit, type = "pearson")[, "Estimate"]
sum(pr^2) / df.residual(fit)
#> [1] 0.9176762

# one-step-ahead quantile residuals are standard normal under a
# correctly specified model, whatever the family
r <- residuals(fit, type = "osa")[, "Estimate"]
qqnorm(r); qqline(r)
```
