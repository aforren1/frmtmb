# Convert draws to a posterior draws object

Convert draws to a posterior draws object

## Usage

``` r
as_draws(x, ...)

# S3 method for class 'frmtmb_draws'
as_draws(x, ...)

# S3 method for class 'frmtmb_draws'
as.array(x, ...)

as_draws_matrix(x, ...)

# S3 method for class 'frmtmb_draws'
as_draws_matrix(x, ...)

as_draws_array(x, ...)

# S3 method for class 'frmtmb_draws'
as_draws_array(x, ...)

as_draws_df(x, ...)

# S3 method for class 'frmtmb_draws'
as_draws_df(x, ...)

as_draws_list(x, ...)

# S3 method for class 'frmtmb_draws'
as_draws_list(x, ...)

as_draws_rvars(x, ...)

# S3 method for class 'frmtmb_draws'
as_draws_rvars(x, ...)

as.mcmc(x, ...)

# S3 method for class 'frmtmb_draws'
as.mcmc(x, combine_chains = FALSE, ...)

# S3 method for class 'frmtmb_multiple'
as_draws(x, ...)
```

## Arguments

- x:

  A `frmtmb_draws` object.

- ...:

  Unused.

- combine_chains:

  If `TRUE`, one `mcmc` object over the pooled draws; otherwise an
  `mcmc.list` with one component per chain, which is what coda's
  diagnostics (`gelman.diag()`) need.

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
#> frm_sample(): sampling stays centered: no random-effect block of this model has a non-centered form:
#>   1 | g [us]: its variance parameter has a flat prior here, and a non-centered chain walks the flat tail that opens at sd = 0. Give it a prior, set_prior(class = "sd"), which the formula interface supplies for you
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
#> [1] "Intercept"       "x"               "sigma_Intercept" "b[1]"           
#> [5] "b[2]"            "b[3]"           
# }
```
