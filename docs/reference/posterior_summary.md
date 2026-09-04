# Summaries and intervals of draws

`posterior_summary()` reduces a matrix of draws to estimate, error and
quantiles in brms's column layout (`Estimate`, `Est.Error`, `Q2.5`,
`Q97.5`). Variables are columns and draws are rows, which is the layout
every draws object in this ecosystem converts to.

## Usage

``` r
posterior_summary(object, ...)

# Default S3 method
posterior_summary(object, probs = c(0.025, 0.975), robust = FALSE, ...)
```

## Arguments

- object:

  A matrix of draws, variables in columns.

- ...:

  Passed to methods.

- probs:

  Quantiles to report.

- robust:

  If `TRUE`, median and MAD instead of mean and SD.

## Value

A matrix with one row per variable.

## Details

The method for posterior draws of a fitted model is in the
`frmtmb.sample` package, along with the sampler that produces them.

## Examples

``` r
# any matrix of draws: rows are draws, columns are variables
m <- cbind(a = rnorm(500), b = rnorm(500, 2))
posterior_summary(m)
#>     Estimate Est.Error       Q2.5    Q97.5
#> a 0.01077349  1.000574 -1.9062756 2.093029
#> b 1.94079266  0.979340  0.0894709 3.910332
posterior_summary(m, robust = TRUE)
#>      Estimate Est.Error       Q2.5    Q97.5
#> a -0.01590824 0.9905923 -1.9062756 2.093029
#> b  1.91763189 1.0087451  0.0894709 3.910332
```
