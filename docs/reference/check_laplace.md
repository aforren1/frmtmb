# Check the Laplace/Wald approximation against NUTS

Samples the fitted objective (see
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md))
and compares the ML estimates and sdreport standard errors against
posterior means and SDs. Close agreement supports the Laplace
approximation and Wald intervals; a posterior SD much larger than the
Wald SE, or a shifted mean, flags parameters where they are unreliable
(typically variance components with few groups).

## Usage

``` r
check_laplace(fit, chains = 2, iter = 1000, ...)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- chains, iter:

  Passed to
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

- ...:

  Passed to
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

## Value

A data frame (one row per outer parameter) with columns `ml`,
`post_mean`, `wald_se`, `post_sd`, `z_shift` ((post_mean - ml)/post_sd)
and `sd_ratio` (post_sd/wald_se).

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
# a binary GLMM with small clusters: the regime where the Laplace
# approximation and Wald intervals are least reliable
set.seed(4)
dd <- data.frame(x = rnorm(120), g = factor(rep(1:30, 4)))
dd$y <- rbinom(120, 1,
               plogis(0.3 + 0.5 * dd$x + rnorm(30, 0, 1)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + bernoulli(), data = dd)

cl <- check_laplace(fit, chains = 1, iter = 500, refresh = 0)
cl
# |z_shift| well above 0 or sd_ratio far from 1 marks the parameters
# whose Wald interval to replace with a profile or bootstrap one
cl[abs(cl$z_shift) > 0.3 | cl$sd_ratio > 1.3, ]
}
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.42, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Laplace/Wald approximation questionable for: theta_1
#>   parameter         ml post_mean   wald_se   post_sd    z_shift sd_ratio
#> 3   theta_1 0.09859697  0.105007 0.3522135 0.5948881 0.01077513 1.688998
# }
```
