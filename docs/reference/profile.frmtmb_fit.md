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

  Parameter names or indices. Required; profiling is not free, so there
  is no all-parameters default. The names are the ones
  [`confint()`](https://rdrr.io/r/stats/confint.html) takes, in any of
  its three spellings: internal (`theta_1`, `tarsus_(Intercept)`),
  parenthesis-free (`tarsus_Intercept`), or a one-to-one natural-scale
  alias (`sd_dam__Intercept`). The profile is of the internal parameter
  either way - a log standard deviation, not a standard deviation - and
  the returned element keeps the internal name. For a profile of a
  natural-scale quantity itself, including one that mixes several
  parameters, use `hypothesis(method = "profile")`.

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

# parameter names are the confint() row names, with or without the
# parentheses, or a one-to-one natural-scale alias
rownames(confint(fit))
#> [1] "(Intercept)"       "x"                 "sigma_(Intercept)"
#> [4] "theta_1"          
pr <- profile(fit, "theta_1")
identical(profile(fit, "(Intercept)"), profile(fit, "Intercept"))
#> [1] TRUE
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
