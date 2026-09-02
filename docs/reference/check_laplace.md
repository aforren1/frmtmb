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

## Details

This is a diagnostic tool: it explores the LIKELIHOOD, with flat priors,
which is what makes the comparison against the ML mode and its Wald
standard errors meaningful.
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
on a formula is the sampling tool instead: it samples a POSTERIOR, under
brms's default priors. A default prior here would change the very thing
being measured, so `check_laplace()` never sets one.

That is also why it samples CENTERED.
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)'s
non-centered parameterization (see Reparameterization there) is offered
only for blocks whose variance parameters carry a prior, and here none
do: the flat prior that makes the comparison meaningful is exactly the
one that leaves a flat tail at `sd = 0` for a non-centered chain to walk
into. So the default costs this function nothing and changes nothing
about it. Give the variance parameters a prior through `priors =` and
the run non-centers; but then it is measuring the Laplace approximation
of a different posterior, which is usually not the question.

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
#> frm_sample(): sampling stays centered: no random-effect block of this model has a non-centered form:
#>   1 | g [us]: its variance parameter has a flat prior here, and a non-centered chain walks the flat tail that opens at sd = 0. Give it a prior, set_prior(class = "sd"), which the formula interface supplies for you
#> Warning: The largest R-hat is 1.08, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> check_laplace(): the chain mixed too poorly to judge the approximation (bulk ESS under 100 for Intercept, theta_1). Rerun with more iterations before reading z_shift or sd_ratio
#> [1] parameter ml        post_mean wald_se   post_sd   z_shift   sd_ratio 
#> [8] ess_bulk 
#> <0 rows> (or 0-length row.names)
# }
```
