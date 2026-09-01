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
