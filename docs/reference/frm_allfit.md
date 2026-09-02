# Refit a model with every available optimizer

The lme4 `allFit()` idea: reuse the assembled design and refit under
each optimizer, to separate optimizer trouble from model
misspecification. Uses nlminb and optim (L-BFGS-B) always, plus bobyqa
(minqa) and NLopt L-BFGS (nloptr) when those packages are installed.

## Usage

``` r
frm_allfit(fit, optimizers = NULL, ...)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- optimizers:

  Named list of optimizers (names or functions, as in
  [`frmtmb_control()`](https://aforren1.github.io/frmtmb/reference/frmtmb_control.md)).
  Default: everything available.

- ...:

  Unused.

## Value

A `frmtmb_allfit` object: `$fits` (the refits, `NULL` where one errored)
and a printed comparison of logLik, convergence, and fixed-effect
spread.

## Examples

``` r
set.seed(6)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rpois(60, exp(0.3 + 0.4 * dd$x + rnorm(6, 0, 0.4)[dd$g]))
fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)
frm_allfit(fit)
#>     optimizer    logLik convergence seconds
#>        nlminb -100.9269           0    0.02
#>         optim -100.9269           0    0.02
#>        bobyqa -100.9269           0    0.03
#>  nloptr_lbfgs -100.9269           0    0.03
#> 
#> logLik spread: 3.82e-11 
#> max fixed-effect spread: 4.27e-07 
```
