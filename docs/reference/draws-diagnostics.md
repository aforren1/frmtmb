# Sampler diagnostics and MCMC plots

`nuts_params()`, `log_posterior()`, `rhat()` and `neff_ratio()` delegate
to bayesplot's `stanfit` methods on the `stanfit` inside the draws
object, so every `bayesplot::mcmc_nuts_*()` display works. `mcmc_plot()`
is brms's spelling for "call a bayesplot `mcmc_*` function on these
draws"; [`pairs()`](https://rdrr.io/r/graphics/pairs.html) is
[`bayesplot::mcmc_pairs()`](https://mc-stan.org/bayesplot/reference/MCMC-scatterplots.html).

## Usage

``` r
mcmc_plot(object, ...)

# S3 method for class 'frmtmb_draws'
mcmc_plot(object, type = "intervals", variable = NULL, ...)

# S3 method for class 'frmtmb_draws'
pairs(x, variable = NULL, ...)

nuts_params(object, ...)

# S3 method for class 'frmtmb_draws'
nuts_params(object, ...)

log_posterior(object, ...)

# S3 method for class 'frmtmb_draws'
log_posterior(object, ...)

rhat(object, ...)

# S3 method for class 'frmtmb_draws'
rhat(object, ...)

neff_ratio(object, ...)

# S3 method for class 'frmtmb_draws'
neff_ratio(object, ...)
```

## Arguments

- object, x:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

- ...:

  Passed to the bayesplot function.

- type:

  The bayesplot function to call, without the `mcmc_` prefix (default
  `"intervals"`).

- variable:

  Variables to plot; defaults to everything except the group-level modes
  and `lp__`.

## Value

A ggplot object, or the diagnostic data frame / vector bayesplot
returns.

## Details

The parameter names bayesplot sees are the frmtmb draws-side names (no
parentheses), not Stan's `par[1]`, because
[`as.array()`](https://rdrr.io/r/base/array.html) relabels them, except
in `nuts_params()`, `rhat()` and `neff_ratio()`, which read the
`stanfit` directly and therefore show Stan's own names.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    requireNamespace("bayesplot", quietly = TRUE)) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                   data = dd, chains = 1, iter = 500, refresh = 0)
  mcmc_plot(ds)
  mcmc_plot(ds, type = "trace", variable = "x")
  head(rhat(ds))
}
#> frm_sample(): default priors (brms 2.23 defaults; priors = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#>   beta[1]   beta[2]     betad      b[1]      b[2]      b[3] 
#> 1.0496060 0.9964106 0.9961408 0.9991621 0.9973939 1.0111695 
# }
```
