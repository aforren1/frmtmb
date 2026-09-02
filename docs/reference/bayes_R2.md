# Bayesian R-squared

The proportion of the outcome's variance the model explains, computed
per draw and returned as a posterior distribution. This is the
residual-based estimator of Gelman, Goodrich, Gabry and Vehtari (2019),
*R-squared for Bayesian regression models*, The American Statistician
73(3), which brms implements as
`var(ypred_s) / (var(ypred_s) + var(y - ypred_s))` for each draw `s`,
with both variances taken over the observations. frmtmb computes the
same expression on
[`posterior_epred()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
draws, so the two agree up to Monte Carlo error on the same model.

## Usage

``` r
bayes_R2(object, ...)

# S3 method for class 'frmtmb_draws'
bayes_R2(
  object,
  resp = NULL,
  summary = TRUE,
  probs = c(0.025, 0.975),
  ndraws = NULL,
  ...
)
```

## Arguments

- object:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

- ...:

  Unused.

- resp:

  For a multivariate model, which response.

- summary:

  If `TRUE` (the default, as in brms), summarize the draws into
  estimate, error and quantiles; if `FALSE`, return the `ndraws x 1`
  matrix of R-squared draws.

- probs:

  Quantiles for the summary.

- ndraws:

  Number of draws to use (default: all).

## Value

A one-row summary matrix, or the matrix of draws when `summary = FALSE`.

## Details

Because the denominator is a per-draw variance rather than a fixed total
sum of squares, the value is bounded in `(0, 1)` by construction and
does not have the classical R-squared's habit of exceeding 1 on
posterior draws.

## See also

[`log_lik()`](https://aforren1.github.io/frmtmb/reference/log_lik.md),
[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.md)

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
  bayes_R2(ds)
  quantile(bayes_R2(ds, summary = FALSE), c(0.1, 0.9))
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
#>       10%       90% 
#> 0.1905620 0.4227541 
# }
```
