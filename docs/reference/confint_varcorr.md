# Natural-scale confidence intervals for covariance parameters

Wald intervals for random-effect standard deviations (on the log scale,
back-transformed) and correlations (on the Fisher-z scale,
back-transformed), delta-method-propagated from the internal `theta`
covariance. One row per SD and per correlation of every block.

## Usage

``` r
confint_varcorr(fit, level = 0.95)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- level:

  Confidence level.

## Value

A data frame with columns `block`, `term`, `type`, `estimate`, `lwr`,
`upr`.
