# Posterior mixture-component probabilities

For a
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md),
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
or [`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md) fit,
the posterior probability that each observation came from each
component, propagating the uncertainty in the parameters: the fit-side
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
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
  [`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md).

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
[`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md)) the rows
are groups, as in
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md).

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
  set.seed(4)
  dd <- data.frame(y = c(rnorm(60, -2), rnorm(60, 3)))
  fit <- frm(bf(y ~ 1), family = mixture(gaussian(), gaussian()),
             data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
  head(pp_mixture(ds)[, "Estimate", ])
}
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>         class1       class2
#> [1,] 0.9999825 1.752471e-05
#> [2,] 0.9999991 9.463510e-07
#> [3,] 0.9996757 3.242515e-04
#> [4,] 0.9999128 8.722546e-05
#> [5,] 0.9889037 1.109632e-02
#> [6,] 0.9998687 1.313094e-04
# }
```
