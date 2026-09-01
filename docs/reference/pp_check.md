# Predictive check against simulated responses

The frequentist analog of brms's `pp_check()`: responses are simulated
from the fitted model (marginally over the random effects) and handed to
the corresponding bayesplot `ppc_*` function (bayesplot must be
installed, but not necessarily attached).

## Usage

``` r
pp_check(object, ...)

# S3 method for class 'frmtmb_fit'
pp_check(object, type = "dens_overlay", ndraws = 10, re.form = NA, ...)

# S3 method for class 'frmtmb_draws'
pp_check(object, type = "dens_overlay", ndraws = 50, ...)
```

## Arguments

- object:

  A `frmtmb_fit` for a univariate model.

- ...:

  Passed to the `ppc_*` function.

- type:

  The bayesplot check, i.e. the part after `ppc_` (`"dens_overlay"`,
  `"hist"`, `"stat"`, `"scatter_avg"`, ...).

- ndraws:

  Number of simulated response vectors.

- re.form:

  Passed to [`simulate()`](https://rdrr.io/r/stats/simulate.html); the
  default `NA` simulates new random effects.

## Value

A ggplot object, as returned by the bayesplot `ppc_*` function that
`type` selects.
