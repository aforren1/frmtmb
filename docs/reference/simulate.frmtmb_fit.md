# Simulate responses from a frmtmb fit

A [`trunc()`](https://rdrr.io/r/base/Round.html)ed response simulates by
rejection within its bounds, so every draw lies in `[lb, ub]` and
posterior-predictive checks
([`dharma_residuals()`](dharma_residuals.md), `pp_check()`) see the same
support the likelihood was normalized on.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
simulate(object, nsim = 1, seed = NULL, re.form = NULL, ...)
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

- ...:

  Unused.

## Value

A data frame with `nsim` columns and a `"seed"` attribute.
