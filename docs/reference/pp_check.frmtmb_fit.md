# Predictive check against simulated responses

The frequentist analog of brms's `pp_check()`: responses are simulated
from the fitted model (marginally over the random effects) and handed to
the corresponding bayesplot `ppc_*` function. Requires bayesplot; call
as `bayesplot::pp_check(fit)` or load bayesplot first.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
pp_check(object, type = "dens_overlay", ndraws = 10, re.form = NA, ...)
```

## Arguments

- object:

  A `frmtmb_fit` for a univariate model.

- type:

  The bayesplot check, i.e. the part after `ppc_` (`"dens_overlay"`,
  `"hist"`, `"stat"`, `"scatter_avg"`, ...).

- ndraws:

  Number of simulated response vectors.

- re.form:

  Passed to [`simulate()`](https://rdrr.io/r/stats/simulate.html); the
  default `NA` simulates new random effects.

- ...:

  Passed to the `ppc_*` function.
