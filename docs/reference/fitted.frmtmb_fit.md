# Fitted values

The modelled response at the estimates, conditional on the random-effect
modes. Equal to `predict(object, type = "response")` on the training
data, for every family.

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
  ...
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

- ...:

  Refused. An argument this method does not have is an error naming it,
  because a swallowed `re_formula` returned the conditional fit and said
  nothing.

## Value

A numeric vector of expected responses; for an ordinal family
([`cumulative()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`sratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`cratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md))
an `n x K` matrix of category probabilities.

## Arguments brms has and this does not

`fitted.brmsfit()` summarizes posterior draws, so it also takes
`ndraws`, `draw_ids`, `sort`, `summary`, `robust` and `probs`. A maximum
likelihood fit has no draws to thin or summarize, so each of those is
refused by name and the message says where the argument does work:
`frmtmb.sample`'s `posterior_epred()`. `nlpar` is refused too; a
non-linear parameter is reached through `dpar`.

## Ordinal responses

An ordinal response has no mean, so
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
`n x K` matrix of category probabilities, with the response's own factor
levels as column names and rows summing to one - the brms
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) convention.
`cs()` terms are honored. The latent linear predictor, which is where
the coefficients live and where `se.fit` is available, is
`predict(object, type = "link")`.

## See also

[`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md),
[`residuals.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/residuals.frmtmb_fit.md),
[frmtmb-scales](https://aforren1.github.io/frmtmb/reference/frmtmb-scales.md)
for which scale each method reports

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x))
fit <- frm(bf(y ~ x) + poisson(), data = dd)
max(abs(fitted(fit) - predict(fit, type = "response")))
#> [1] 0
```
