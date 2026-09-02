# Size of a draws object

`ndraws()` counts the post-warmup draws (all chains pooled),
`niterations()` the draws per chain, `nchains()` the chains and
`nvariables()` the sampled parameters. The names and meanings are
posterior's; frmtmb registers methods with posterior so that the
generics work whether or not that package is attached.

## Usage

``` r
ndraws(x)

# S3 method for class 'frmtmb_draws'
ndraws(x)

nchains(x)

# S3 method for class 'frmtmb_draws'
nchains(x)

niterations(x)

# S3 method for class 'frmtmb_draws'
niterations(x)

nvariables(x)

# S3 method for class 'frmtmb_draws'
nvariables(x)
```

## Arguments

- x:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

## Value

A single integer.

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
  c(ndraws(ds), niterations(ds), nchains(ds), nvariables(ds))
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
#> [1] 250 250   1  11
# }
```
