# Posterior mixture-component probabilities

For a
[`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html),
[`frmtmb::mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.html)
or
[`frmtmb.latent::lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.html)
fit, the posterior probability that each observation came from each
component, propagating the uncertainty in the parameters: the fit-side
[`frmtmb::mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.html)
computation is run at every draw. brms calls this `pp_mixture()`.

## Usage

``` r
pp_mixture(x, ...)

# S3 method for class 'frmtmb_draws'
pp_mixture(x, summary = TRUE, ndraws = NULL, ...)
```

## Arguments

- x:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).

- ...:

  Unused.

- summary:

  If `TRUE` (the default, as in brms), an
  `observations x statistics x components` array of summaries; otherwise
  the raw `draws x observations x components` array.

- ndraws:

  Number of draws to use (default: all).

## Value

An array; see `summary`. For a group-level mixture
(`mixture(groups = )`,
[`frmtmb.latent::lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.html))
the rows are groups, as in
[`frmtmb::mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.html).

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
  set.seed(4)
  dd <- data.frame(y = c(rnorm(60, -2), rnorm(60, 3)))
  fit <- frm(bf(y ~ 1), family = frmtmb::mixture(gaussian(), gaussian()),
             data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
  head(pp_mixture(ds)[, "Estimate", ])
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0.6, 3.5)
#>   Intercept          student_t(3, 0.6, 3.5)
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#>         class1       class2
#> [1,] 0.9999716 2.842349e-05
#> [2,] 0.9999974 2.586659e-06
#> [3,] 0.9996354 3.646268e-04
#> [4,] 0.9998869 1.131103e-04
#> [5,] 0.9899028 1.009723e-02
#> [6,] 0.9998378 1.622369e-04
# }
```
