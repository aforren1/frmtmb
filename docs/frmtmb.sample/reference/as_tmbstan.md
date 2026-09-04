# Sample from a frmtmb fit with tmbstan (NUTS)

Hands the fitted RTMB object to
[`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html);
from Stan's point of view the model is the (Laplace-free) joint density,
so all parameters including random effects are sampled unless
`laplace = TRUE` is passed through.

## Usage

``` r
as_tmbstan(fit, ...)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- ...:

  Passed to
  [`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html)
  (chains, iter, laplace, ...).

## Value

A `stanfit` object.

## Details

It hands over the fit's OWN objective and nothing else: no default
priors, no non-centering, no named draws.
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
is the route that applies brms's default priors and returns the draws
surface; this one is the escape hatch to tmbstan's own arguments.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# The raw stanfit, for rstan and bayesplot code that wants one.
# Use frm_sample() instead when you want frmtmb parameter names.
# This run is deliberately short, so expect sampler warnings.
sf <- as_tmbstan(fit, chains = 1, iter = 400, refresh = 0)
class(sf)
# every parameter is sampled, random effects included, because Stan
# sees the joint density. Pass laplace = TRUE to integrate them out.
dim(as.matrix(sf))
}
#> Warning: There were 5 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: There were 1 chains where the estimated Bayesian Fraction of Missing Information was low. See
#> https://mc-stan.org/misc/warnings.html#bfmi-low
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: The largest R-hat is 1.1, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> [1] 200  13
# }
```
