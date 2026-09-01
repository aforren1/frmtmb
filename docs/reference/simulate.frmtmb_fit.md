# Simulate responses from a frmtmb fit

Simulate responses from a frmtmb fit

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
