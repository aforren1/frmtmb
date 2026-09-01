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

## Examples

``` r
if (requireNamespace("bayesplot", quietly = TRUE)) {
  set.seed(1)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

  # the observed density against draws from the fit
  pp_check(fit, ndraws = 20)

  # any bayesplot ppc_* check, named by its suffix. A statistic the
  # model was not fitted to is the informative one: here, the share
  # of zeros, which is how zero inflation shows up.
  pp_check(fit, type = "stat", stat = function(y) mean(y == 0),
           ndraws = 50)
}
#> `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```
