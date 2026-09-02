# Fitted values

The modelled response at the estimates, conditional on the random-effect
modes. Equal to `predict(object, type = "response")` on the training
data, for every family.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
fitted(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

## Value

A numeric vector of expected responses; for an ordinal family
([`cumulative()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`sratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`cratio()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md))
an `n x K` matrix of category probabilities.

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
[`residuals.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/residuals.frmtmb_fit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x))
fit <- frm(bf(y ~ x) + poisson(), data = dd)
max(abs(fitted(fit) - predict(fit, type = "response")))
#> [1] 0
```
