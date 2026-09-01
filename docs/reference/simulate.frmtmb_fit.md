# Simulate responses from a frmtmb fit

A [`trunc()`](https://rdrr.io/r/base/Round.html)ed response simulates by
rejection within its bounds, so every draw lies in `[lb, ub]` and
posterior-predictive checks
([`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md),
`pp_check()`) see the same support the likelihood was normalized on.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
simulate(object, nsim = 1, seed = NULL, re.form = NULL, censored = FALSE, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- nsim:

  Number of simulated response vectors.

- seed:

  Optional RNG seed. Follows the
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) contract:
  the global RNG state is restored afterwards, and the seed used is
  attached as the `"seed"` attribute.

- re.form:

  `NULL` (default) conditions on the estimated random effects; `NA`
  redraws them from their estimated distribution (marginal simulation).

- censored:

  Apply the fitted `cens()` mechanism to the draws (see Censored
  responses). Ignored without `cens()`.

- ...:

  Unused.

## Value

A data frame with `nsim` columns and a `"seed"` attribute.

## Censored responses

On a `cens()` fit the default draws the LATENT, uncensored response: the
model describes the latent distribution, and censoring is a property of
the observation process, not of the response. This matches brms, whose
[`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md)
also ignores `cens()` (and whose `pp_check()` therefore drops the
censored rows). The draws are then not comparable with the observed
values on censored rows, which is why
[`dharma_residuals()`](https://aforren1.github.io/frmtmb/reference/dharma_residuals.md)
and `residuals(type = "osa")` refuse or skip them.

`censored = TRUE` applies the fitted censoring mechanism to each draw
instead, so the draws are directly comparable with the observed data:
every draw is recorded at the edge of the observation window it falls
outside, capped above by the right-censoring point and below by the
left-censoring point. Those points are the response values of the
censored rows, and they must be the same for every censored row on a
side (type-I censoring): with row-varying censoring times an uncensored
row's censoring point is unknown, so the mechanism cannot be applied to
its draws and the call is refused. Interval censoring has no
single-value representation and is refused too.
