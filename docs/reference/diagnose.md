# Convergence diagnostics for a frmtmb fit

Reports the optimizer's own verdict plus four checks that a converged
fit can still fail: non-finite standard errors, complete separation in a
binomial-type fit, predictor columns scaled far from one, and variance
components on the boundary of their parameter space (lme4's
`isSingular()`, read off the estimates rather than the Hessian).

## Usage

``` r
diagnose(fit, quiet = FALSE)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- quiet:

  If `TRUE`, return the diagnostics without printing.

## Value

Invisibly, a list of diagnostics.
