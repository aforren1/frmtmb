# Diagnostic plots for a fit

Panel 1: Pearson residuals against fitted values with a lowess trend.
Panel 2: normal QQ plot of the Pearson residuals. For simulation-based
residuals that are exact for discrete families, use
[`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
or `residuals(type = "osa")`.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
plot(x, which = 1:2, ask = NULL, ...)
```

## Arguments

- x:

  A `frmtmb_fit`.

- which:

  Subset of `1:2`.

- ask:

  Whether to prompt between plots; defaults to the usual
  interactive-device rule.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the plots it draws.

## Details

On an ordinal fit
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is a matrix of
category probabilities, so panel 1 uses the expected category index
`sum_k k * P(y = k)` - the same scalar the Pearson residual is taken
against - and labels the axis accordingly.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# both panels, side by side
op <- par(mfrow = c(1, 2))
plot(fit, ask = FALSE)

par(op)

# just the QQ panel
plot(fit, which = 2)
```
