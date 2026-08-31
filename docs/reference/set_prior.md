# Set up priors brms-style

Builds prior specifications with brms spelling:
`set_prior("normal(0, 5)", class = "b")`. Combine several with `+` or
[`c()`](https://rdrr.io/r/base/c.html). Distributions: `normal(mu, sd)`,
`student_t(df, mu, sd)`, `cauchy(mu, sd)`, `exponential(rate)`; an empty
string sets bounds only.

## Usage

``` r
set_prior(
  prior = "",
  class = "b",
  coef = "",
  dpar = "",
  group = "",
  lb = NA,
  ub = NA
)
```

## Arguments

- prior:

  Distribution string, e.g. `"normal(0, 5)"`, or a
  [`prior_normal()`](frmtmb-priors.md)/[`prior_t()`](frmtmb-priors.md)
  object, or `""` for bounds only.

- class:

  `"b"`, `"Intercept"`, `"sd"`, or `"theta"`.

- coef:

  Restrict to one coefficient (classes `"b"`/`"Intercept"`).

- dpar:

  Distributional parameter (default: the location parameters).

- group:

  Restrict class `"sd"` to one grouping factor.

- lb, ub:

  Optional hard bounds.

## Value

A `frmtmb_priorlist`.

## Details

Classes and their scales:

- `"b"`: population-level coefficients of `dpar` (default: the location
  parameters), excluding the intercept; narrow to one coefficient with
  `coef`. Link scale.

- `"Intercept"`: the intercept of `dpar`. Link scale.

- `"sd"`: random-effect standard deviations (and smoothing SDs), on the
  NATURAL sd scale with the log-Jacobian applied, so
  `set_prior("exponential(1)", class = "sd")` means what it says; narrow
  with `group`. Correlation parameters are not covered (use class
  `"theta"` on the internal scale if you must).

- `"theta"`: raw internal covariance parameters (escape hatch).

When priors overlap, later specifications override earlier ones, so put
class-wide priors first and coefficient-specific ones after. `lb`/`ub`
become hard bounds (via Stan's constrained transforms in
[`frm_sample()`](frm_sample.md)); for class `"sd"` they apply on the sd
scale.
