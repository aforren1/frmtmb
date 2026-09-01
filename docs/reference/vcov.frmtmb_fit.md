# Covariance matrix of the fixed-effect estimates

Covers the estimated coefficients of every linear predictor; dpars fixed
to constants are excluded.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
vcov(object, full = FALSE, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- full:

  If `TRUE`, include covariance parameters (`theta`), named as in
  [`confint()`](https://rdrr.io/r/stats/confint.html) (the glmmTMB
  `vcov(full = TRUE)` convention).

- ...:

  Unused.

## Value

A covariance matrix.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# standard errors of the fixed effects
sqrt(diag(vcov(fit)))
#>       (Intercept)                 x sigma_(Intercept) 
#>        0.27456628        0.11328517        0.07454062 
# the covariance parameters join the block on their internal scale
rownames(vcov(fit, full = TRUE))
#> [1] "(Intercept)"       "x"                 "sigma_(Intercept)"
#> [4] "theta_1"          

# the matrix is what a delta-method calculation needs
V <- vcov(fit)
a <- c(1, 2)                       # prediction at x = 2, no group
sqrt(drop(t(a) %*% V[1:2, 1:2] %*% a))
#> [1] 0.348039
```
