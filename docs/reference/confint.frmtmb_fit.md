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
  method = c("wald", "profile", "uniroot"),
  ...
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- parm:

  Parameter names (see `rownames` of the Wald result) or indices.
  Required for `"profile"` and `"uniroot"`; defaults to all parameters
  for `"wald"`.

- level:

  Confidence level.

- method:

  `"wald"` (fast, from the sdreport covariance), `"profile"` (likelihood
  profile via
  [`TMB::tmbprofile()`](https://rdrr.io/pkg/TMB/man/tmbprofile.html)),
  or `"uniroot"` (likelihood-root search via
  [`TMB::tmbroot()`](https://rdrr.io/pkg/TMB/man/tmbroot.html)).

- ...:

  Passed to the TMB profiling functions.

## Value

A matrix with columns `lwr`, `upr`, `est`.
