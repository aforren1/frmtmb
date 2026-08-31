# Diagnostic plots for a fit

Panel 1: Pearson residuals against fitted values with a lowess trend.
Panel 2: normal QQ plot of the Pearson residuals. For simulation-based
residuals that are exact for discrete families, use
[`dharma_residuals()`](dharma_residuals.md) or
`residuals(type = "osa")`.

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
