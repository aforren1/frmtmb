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

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)

# one row per SD and per correlation, on the scale they are read on
confint_varcorr(fit)
#>   block               term type   estimate        lwr       upr
#> 1 x | g        (Intercept)   sd 0.88616477  0.6246517 1.2571614
#> 2 x | g                  x   sd 0.45121347  0.2633493 0.7730933
#> 3 x | g cor((Intercept),x)  cor 0.05273605 -0.4893613 0.5654374

# confint() reports the same parameters on their internal scale, so
# the bounds there are log-SDs and Fisher-z correlations
confint(fit)[grep("^theta", rownames(confint(fit))), ]
#>                lwr        upr         est
#> theta_1 -0.4705610  0.2288563 -0.12085237
#> theta_2 -1.3342739 -0.2573555 -0.79581472
#> theta_3 -0.5360151  0.6416342  0.05280954

# a fit with no random effects has no covariance parameters
confint_varcorr(frm(bf(y ~ x) + gaussian(), data = dd))
#> NULL
```
