# Confidence intervals for frmtmb fits

Covariance parameters (`theta_*`) are reported on their internal
(unconstrained) scale.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
confint(
  object,
  parm = NULL,
  level = 0.95,
  method = c("wald", "Wald", "profile", "uniroot", "boot"),
  nsim = 500,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- parm:

  Parameter names (see `rownames` of the Wald result) or indices.
  Required for `"profile"` and `"uniroot"`; defaults to all parameters
  for `"wald"` and `"boot"`.

- level:

  Confidence level.

- method:

  `"wald"` (fast, from the sdreport covariance), `"profile"` (likelihood
  profile via
  [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html)),
  `"uniroot"` (likelihood-root search via
  [`TMB::tmbroot()`](https://rdrr.io/pkg/TMB/man/tmbroot.html)), or
  `"boot"` (parametric-bootstrap percentile intervals through
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md),
  the `lme4::confint(method = "boot")` analog; like the other methods it
  works on the internal parameter scale). `"Wald"` is accepted as an
  alias for `"wald"`.

- nsim, seed:

  Bootstrap draws and seed for `method = "boot"`.

- ...:

  Passed to the TMB profiling functions, or to
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  for `method = "boot"` (e.g. `re.form`).

## Value

A matrix with columns `lwr`, `upr`, `est`.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# Wald intervals for every parameter, covariance ones included
confint(fit)
#>                          lwr       upr         est
#> (Intercept)        0.7014216 1.7777016  1.23956159
#> x                  0.4486521 0.8927218  0.67068694
#> sigma_(Intercept) -0.1808711 0.1113227 -0.03477419
#> theta_1           -0.7094981 0.2925625 -0.20846781

# the likelihood profile does not assume a quadratic log-likelihood,
# so it is the one to trust for a variance component
confint(fit, parm = "theta_1", method = "profile")
#>                lwr       upr        est
#> theta_1 -0.6835298 0.3479822 -0.2084678

# confint_varcorr() puts the same information on the SD scale
confint_varcorr(fit)
#>   block        term type  estimate      lwr      upr
#> 1 1 | g (Intercept)   sd 0.8118272 0.491891 1.339856
# \donttest{
# a parametric bootstrap, the lme4 confint(method = "boot") analog
confint(fit, parm = "x", method = "boot", nsim = 50, seed = 1)
#>         lwr       upr       est
#> x 0.4314735 0.8223053 0.6706869
# }
```
