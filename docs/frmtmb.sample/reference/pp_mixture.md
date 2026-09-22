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
pp_mixture(
  x,
  newdata = NULL,
  re_formula = arg_unset(),
  resp = NULL,
  ndraws = NULL,
  draw_ids = NULL,
  log = FALSE,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- x:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).

- ...:

  Unused.

- newdata, re_formula:

  Accepted for brms's signature and refused:
  [`frmtmb::mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.html)
  is a statement about the rows the model was fitted on.

- resp:

  The response whose mixture to report, for a multivariate model.

- ndraws:

  Number of draws to use (default: all).

- draw_ids:

  The draws to use, by row index, instead of the evenly spaced subsample
  `ndraws` takes.

- log:

  If `TRUE`, log probabilities.

- summary:

  If `TRUE` (the default, as in brms), an
  `observations x statistics x components` array of summaries; otherwise
  the raw `draws x observations x components` array.

- robust:

  If `TRUE`, median and MAD instead of mean and SD.

- probs:

  The two quantiles the summary reports.

## Value

An array; see `summary`. For a group-level mixture
(`mixture(groups = )`,
[`frmtmb.latent::lca()`](https://aforren1.github.io/frmtmb/frmtmb.latent/reference/lca.html))
the rows are groups, as in
[`frmtmb::mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.html).

## Details

The argument order is brms's, so `summary` sits in brms's own eighth
position and not in the second: the second is `newdata`.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
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
#> [1,] 0.9999820 1.798661e-05
#> [2,] 0.9999990 9.966469e-07
#> [3,] 0.9996895 3.104748e-04
#> [4,] 0.9999137 8.634107e-05
#> [5,] 0.9901016 9.898355e-03
#> [6,] 0.9998713 1.286582e-04
# }
```
