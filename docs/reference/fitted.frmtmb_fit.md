# Fitted values

brms's [`fitted()`](https://rdrr.io/r/stats/fitted.values.html): a
summary of the expected response, in the columns `Estimate`,
`Est.Error`, `Q2.5` and `Q97.5`. The estimate is the modelled response
at the estimates, conditional on the random-effect modes, which is
`frm_linpred(object, type = "response")` on the training data for every
family.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
fitted(
  object,
  newdata = NULL,
  re_formula = NULL,
  scale = c("response", "linear"),
  resp = NULL,
  dpar = NULL,
  nlpar = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...,
  allow_new_levels = FALSE
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- newdata:

  Optional data frame to evaluate on. Defaults to the training data.

- re_formula:

  `NULL` (default) keeps the random effects, so the answer is
  conditional on the modes; `NA` or `~0` gives the population-level
  answer. brms's spelling, and the only one: lme4's `re.form` is not
  accepted here.

- scale:

  `"response"` (default) for the modelled response, or `"linear"` for
  the linear predictor. brms's spelling of what
  [`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
  calls `type`.

- resp:

  For multivariate fits: which response (defaults to the first).

- dpar:

  Which distributional parameter to report instead of the mean.

- nlpar:

  brms's name for a non-linear parameter, which is a distributional
  parameter here: a synonym for `dpar`, and giving both is an error.

- ndraws, draw_ids, sort, summary, robust:

  brms's arguments, in brms's positions so that a positional brms call
  asks the same question. Each needs posterior draws, and a
  maximum-likelihood fit has none, so each is refused by name with the
  reason and with the place it does work: `frmtmb.sample`'s
  `posterior_epred()`. The default of each is accepted and changes
  nothing.

- probs:

  Probabilities of the two quantile columns. brms takes the posterior
  quantiles there; here they are the ends of the Wald interval at those
  probabilities.

- ...:

  Refused. An argument this method does not have is an error naming it,
  because a swallowed `re_formula` returned the conditional fit and said
  nothing.

- allow_new_levels:

  Predict unseen grouping-factor levels at the population level instead
  of erroring, as in
  [`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md).

## Value

An `n x 4` matrix with the columns `Estimate`, `Est.Error` and one per
entry of `probs`. For an ordinal or categorical family an `n x 4 x K`
array, the third dimension named `P(Y = k)`, which is brms's shape. The
ROW dimnames are `NULL`, as brms's are; the data's row names are on
[`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md)
and [`model.frame()`](https://rdrr.io/r/stats/model.frame.html).

## Details

A maximum-likelihood fit has no draws to summarize, so `Est.Error` is
the delta-method standard error
[`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md)
reports for the same quantity, and the `Q` columns are the Wald interval
at those probabilities. The interval is around the EXPECTED response and
carries no observation noise;
[`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
is the predictive interval that does.

A predictor whose standard error this package cannot produce, such as a
nonlinear body or a structured likelihood, reports `NA` in `Est.Error`
and in the `Q` columns rather than a number it does not have. The
estimate is unaffected.

## Ordinal responses

An ordinal response has no mean, so
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) summarizes the
`K` category probabilities, with rows of `Estimate` summing to one,
which is the brms
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) convention.
`cs()` terms are honored. The standard error of a category probability
is the finite-difference delta method over the whole outer parameter
vector, because the probability depends on the thresholds and the `cs()`
coefficients as well as on the linear predictor. That route covers the
OUTER parameters only, so on a mixed ordinal fit it does not carry the
conditional variance of the random-effect modes, which the scalar route
does: `se.fit` reads the joint precision and this reads
`vcov(full = TRUE)`. The estimates are unaffected. The latent linear
predictor, which is where the coefficients live, is
`frm_linpred(object, type = "link")`.

## See also

[`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
for the predictive interval,
[`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md)
for the linear predictor,
[`residuals.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/residuals.frmtmb_fit.md),
and
[frmtmb-scales](https://aforren1.github.io/frmtmb/reference/frmtmb-scales.md)
for which scale each method reports

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x))
fit <- frm(bf(y ~ x) + poisson(), data = dd)
head(fitted(fit))
#>       Estimate Est.Error      Q2.5     Q97.5
#> [1,] 0.7751024 0.1130250 0.5535775 0.9966273
#> [2,] 1.1849807 0.1126043 0.9642804 1.4056811
#> [3,] 0.6946365 0.1132732 0.4726251 0.9166479
#> [4,] 2.4828747 0.3477984 1.8012024 3.1645469
#> [5,] 1.2791025 0.1161520 1.0514488 1.5067563
#> [6,] 0.7001766 0.1132798 0.4781522 0.9222009
max(abs(fitted(fit)[, "Estimate"] -
          frm_linpred(fit, type = "response")))
#> [1] 0
```
