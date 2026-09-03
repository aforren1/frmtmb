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
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 1, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>                        mean         sd       2.5%       97.5%     n_eff
#> Intercept        0.84891010 0.16309215  0.5687758  1.25173865 184.25129
#> x                0.60945161 0.11938398  0.4017843  0.83215316 371.49544
#> sigma_Intercept -0.03979135 0.08232127 -0.1897263  0.09976085 202.39185
#> theta_1         -1.60326585 0.96638523 -4.0529050 -0.27698417  81.16065
#>                      Rhat
#> Intercept       0.9991088
#> x               1.0009587
#> sigma_Intercept 1.0026390
#> theta_1         1.0276044
# }
```
