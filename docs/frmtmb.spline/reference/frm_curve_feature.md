# Features of a fitted curve: a peak, a trough, or a level crossing

Locates a stationary point or a level crossing of a fitted curve and
gives its position a standard error, by the implicit-function delta
method.

## Usage

``` r
frm_curve_feature(
  object,
  var,
  type = c("maximum", "minimum", "extremum", "crossing"),
  at = 0,
  newdata = NULL,
  dpar = NULL,
  resp = NULL,
  re.form = NA,
  level = 0.95,
  eps = NULL,
  maxit = 50L,
  tol = 1e-06
)
```

## Arguments

- object:

  A `frmtmb_fit`, or a `frmtmb_curve` from
  [`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md).

- var:

  Name of the covariate the feature is located along.

- type:

  `"maximum"`, `"minimum"` or `"extremum"` for a stationary point of the
  given kind; `"crossing"` for the points where the curve passes `at`.

- at:

  The level to cross, for `type = "crossing"`.

- newdata:

  The grid the search scans, and the values every other covariate is
  held at. Required when `object` is a fit.

- dpar:

  Distributional parameter to read the curve off. `NULL`, the default,
  is the location parameter `mu`.

- resp:

  Response name, for a multivariate fit.

- re.form:

  `NA` (the default) evaluates the population curve, the convention
  `mgcv` and `gratia` plot. `NULL` keeps every random effect, so the
  grid must carry the grouping columns and the curve is that group's
  own.

- level:

  Coverage of both intervals.

- eps:

  Step size for the differences. `NULL` is the measured default of
  [`frm_curve_deriv()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_deriv.md).

- maxit:

  Newton iterations allowed per root.

- tol:

  Largest relative disagreement with `predict(se.fit = TRUE)` the
  assembled covariance may show before the call refuses.

## Value

A data frame with one row per root: `.feature`, `.var`, `.estimate` (the
located position), `.se`, `.lower_ci`, `.upper_ci`, `.value` (the curve
there) and `.value_se`. The `"check"` attribute carries the covariance
agreement, as
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)'s
does.

## Details

The search is a scan of the grid for a sign change, then Newton
refinement on the fitted curve, then one design pass at the located root
for the variance. Every root the grid brackets is returned, so a curve
with two peaks gives two rows; a curve with none gives a zero-row result
rather than an error, because "this curve does not peak in this window"
is an answer.

At a stationary point the curve's own value has a delta-method
simplification worth knowing: the derivative of `f(t*)` with respect to
the coefficients is `df/dc + f'(t*) dt*/dc`, and `f'(t*)` is zero there,
so the standard error of the PEAK HEIGHT is just the pointwise standard
error of the curve at `t*`. It is reported as `.value_se`, and it is not
inflated by the uncertainty in the peak's location.

## See also

[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md),
[`frm_curve_deriv()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_deriv.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(300)))
dd$y <- 2 * sin(pi * dd$x) + rnorm(300, 0, 0.3)
fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
                   family = stats::gaussian(), data = dd)
# 2 sin(pi x) peaks at x = 0.5
frm_curve_feature(fit, var = "x", type = "maximum",
                  newdata = data.frame(x = seq(0.05, 0.95, length.out = 41)))
#> <frmtmb curve feature> maximum, 1 found, level 0.95
#>   covariance checked against predict(se.fit = TRUE) to 5.11e-15 relative
#>   .feature .var .estimate       .se .lower_ci .upper_ci   .value .value_se
#> 1  maximum    x 0.4911416 0.0163341 0.4591274 0.5231559 2.053619 0.0409223
```
