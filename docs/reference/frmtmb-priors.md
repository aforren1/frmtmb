# Prior specifications for frm_sample

Priors apply on the INTERNAL parameter scale: coefficients are on their
link scale, and covariance parameters (`theta_*`) are the unconstrained
parameterization (log-SDs, scaled-Cholesky terms), so
`prior_normal(0, 1)` on `theta_1` is a lognormal prior on that SD.

## Usage

``` r
prior_normal(location = 0, scale = 1)

prior_t(df = 3, location = 0, scale = 1)

prior_lkj(eta = 1)
```

## Arguments

- location, scale, df:

  Prior parameters.

- eta:

  LKJ shape. `1` is uniform over correlation matrices, larger values
  concentrate toward the identity, and `0 < eta < 1` pushes toward the
  boundary.

## Value

A `frmtmb_prior` object.

## Examples

``` r
# the objects themselves are cheap descriptions
prior_normal(0, 2)
#> $kind
#> [1] "normal"
#> 
#> $location
#> [1] 0
#> 
#> $scale
#> [1] 2
#> 
#> attr(,"class")
#> [1] "frmtmb_prior"
prior_t(df = 3, location = 0, scale = 1)
#> $kind
#> [1] "t"
#> 
#> $df
#> [1] 3
#> 
#> $location
#> [1] 0
#> 
#> $scale
#> [1] 1
#> 
#> attr(,"class")
#> [1] "frmtmb_prior"

# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# names are parameter names as in the draws, or whole components.
# theta_1 is a log-SD, so a normal there is a lognormal on the SD.
ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0,
                 prior = list(beta = prior_normal(0, 5),
                               theta_1 = prior_t(3, 0, 1)))
summary(ds)
}
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>                        mean         sd       2.5%      97.5%    n_eff      Rhat
#> Intercept        0.85276138 0.13582704  0.5625223  1.1209474 111.4628 0.9961619
#> x                0.61282556 0.10893028  0.4246522  0.8233411 219.5038 1.0043144
#> sigma_Intercept -0.04262745 0.08395983 -0.2119155  0.1364832 200.6217 1.0152033
#> theta_1         -1.50527252 0.96086982 -3.9788596 -0.3289263 141.6968 1.0097685
# }
```
