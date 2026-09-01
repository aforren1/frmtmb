# Likelihood profiles

Wraps [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html)
per parameter. The returned objects have
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) methods (from TMB).

## Usage

``` r
# S3 method for class 'frmtmb_fit'
profile(fitted, parm, ...)
```

## Arguments

- fitted:

  A `frmtmb_fit`.

- parm:

  Parameter names (as in
  [`confint()`](https://rdrr.io/r/stats/confint.html) rownames) or
  indices. Required; profiling is not free, so there is no
  all-parameters default.

- ...:

  Passed to
  [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html).

## Value

A `tmbprofile` data frame, or a named list of them when `parm` has
length above one.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# parameter names are the confint() row names
rownames(confint(fit))
#> [1] "(Intercept)"       "x"                 "sigma_(Intercept)"
#> [4] "theta_1"          
pr <- profile(fit, "theta_1")
plot(pr)

# TMB's confint() reads the interval off the profile
confint(pr)
#>            lower     upper
#> theta -0.6835298 0.3479822

# several parameters at once return a named list
prs <- profile(fit, c("x", "theta_1"))
names(prs)
#> [1] "x"       "theta_1"
```
