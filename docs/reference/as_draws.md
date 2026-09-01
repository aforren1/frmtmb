# Convert draws to a posterior draws object

Convert draws to a posterior draws object

## Usage

``` r
as_draws(x, ...)

# S3 method for class 'frmtmb_draws'
as_draws(x, ...)
```

## Arguments

- x:

  A `frmtmb_draws` object.

- ...:

  Unused.

## Value

A
[`posterior::draws_matrix`](https://mc-stan.org/posterior/reference/draws_matrix.html):
one column per sampled variable and one row per draw.

## Examples

``` r
# \donttest{
if (requireNamespace("posterior", quietly = TRUE)) {
  set.seed(9)
  dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
  dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
  fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)

  # hands the draws to the posterior package, keeping the frmtmb
  # parameter names
  dm <- as_draws(ds)
  posterior::summarise_draws(dm)
  # which is what variables() lists
  head(variables(ds))
}
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> [1] "(Intercept)"       "x"                 "sigma_(Intercept)"
#> [4] "b[1]"              "b[2]"              "b[3]"             
# }
```
