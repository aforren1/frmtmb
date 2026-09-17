# Summaries and intervals of draws

`posterior_summary()` reduces draws to estimate, error and quantiles in
brms's column layout (`Estimate`, `Est.Error`, `Q2.5`, `Q97.5`). A
matrix has variables in columns and draws in rows and gives one row per
variable. A three-dimensional array has draws in its first margin, as
`ranef(summary = FALSE)` returns them, and gives a
`variables x statistics x third-margin` array, as in brms.

## Usage

``` r
posterior_summary(x, ...)

# Default S3 method
posterior_summary(x, probs = c(0.025, 0.975), robust = FALSE, ...)

# S3 method for class 'frmtmb_fit'
posterior_summary(x, ...)

# S3 method for class 'frmtmb_multiple'
posterior_summary(x, ...)
```

## Arguments

- x:

  A matrix or three-dimensional array of draws.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

- probs:

  Quantiles to report.

- robust:

  If `TRUE`, median and MAD instead of mean and SD.

## Value

A matrix with one row per variable, or an array, as above.

## Details

The computation is brms's `posterior_summary.default()`, statistic for
statistic: `mean` and `sd`, or `median` and `mad` when `robust`, then
`quantile`, each with `na.rm = TRUE`.

The method for posterior draws of a fitted model is in the
`frmtmb.sample` package, along with the sampler that produces them. A
maximum-likelihood fit has no draws, and its method says so.

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

# a maximum-likelihood fit has no draws to summarize
dd <- data.frame(y = rnorm(40), x = rnorm(40))
try(posterior_summary(frm(bf(y ~ x) + gaussian(), data = dd)))
#> Error : posterior_summary() needs posterior draws and a frmtmb_fit has none: frm() is maximum likelihood, so it carries one parameter vector rather than chains. Install frmtmb.sample and sample first, with posterior_summary(frmtmb.sample::frm_sample(fit)); for the point estimates and their covariance use fixef(), vcov() or confint() on the fit itself
```
