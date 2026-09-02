# Summaries and intervals of draws

`posterior_summary()` reduces draws to estimate, error and quantiles in
brms's column layout (`Estimate`, `Est.Error`, `Q2.5`, `Q97.5`);
`posterior_interval()` gives the central interval alone, in rstantools'
layout. Both work on a `frmtmb_draws` object and on any matrix of draws,
which is what makes `posterior_summary(bayes_R2(ds, summary = FALSE))`
work.

## Usage

``` r
posterior_summary(object, ...)

# Default S3 method
posterior_summary(object, probs = c(0.025, 0.975), robust = FALSE, ...)

# S3 method for class 'frmtmb_draws'
posterior_summary(
  object,
  probs = c(0.025, 0.975),
  robust = FALSE,
  variable = NULL,
  ...
)

posterior_interval(object, ...)

# S3 method for class 'frmtmb_draws'
posterior_interval(object, prob = 0.95, variable = NULL, ...)

predictive_interval(object, ...)

# S3 method for class 'frmtmb_draws'
predictive_interval(
  object,
  prob = 0.9,
  newdata = NULL,
  resp = NULL,
  re.form = NULL,
  ndraws = NULL,
  ...
)

predictive_error(object, ...)

# S3 method for class 'frmtmb_draws'
predictive_error(object, resp = NULL, re.form = NULL, ndraws = NULL, ...)
```

## Arguments

- object:

  A `frmtmb_draws`, or a matrix of draws (variables in columns).

- ...:

  Unused.

- probs:

  Quantiles for `posterior_summary()`.

- robust:

  If `TRUE`, median and MAD instead of mean and SD.

- variable:

  Optional subset of variables, by name.

- prob:

  Central interval width for `posterior_interval()` and
  `predictive_interval()`.

- ndraws, newdata, re.form, resp:

  Passed to
  [`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md).

## Value

A matrix with one row per variable (or per observation, for the
predictive functions), except `predictive_error()`, which returns a
draws-by-observations matrix.

## Details

`predictive_interval()` is the same central interval of
[`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
draws, and `predictive_error()` is the matrix of predictive residuals
`y - yrep`, one row per draw.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                   data = dd, chains = 1, iter = 500, refresh = 0)
  posterior_summary(ds, variable = c("Intercept", "x"))
  posterior_interval(ds, prob = 0.9, variable = "x")
  head(predictive_interval(ds))
}
#> frm_sample(): default priors (brms 2.23 defaults; priors = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: The largest R-hat is 1.07, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>              5%      95%
#> [1,] -1.0497754 2.071685
#> [2,] -1.3898044 1.640983
#> [3,] -0.5838004 2.500210
#> [4,] -0.9800667 2.618998
#> [5,] -0.3606254 2.627856
#> [6,] -1.7300234 1.353915
# }
```
