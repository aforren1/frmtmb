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
  contrast = NULL,
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

- contrast:

  A second grid with the same number of rows, or `NULL` for an ordinary
  curve. With it the curve is `newdata` minus `contrast`, row by row.

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

## On a difference curve

A `frmtmb_curve` from `frm_curve(contrast = )` carries its second grid
here, so the feature located is a feature OF THE DIFFERENCE:
`type = "crossing"` with `at = 0` is the question a difference curve is
usually drawn to answer, the place where two curves meet. The delta
method is the same one, on the same covariance: for a crossing the
variance of the located position is the variance of the difference at
that position over the squared slope of the difference.

`var` moves in BOTH grids together, so `contrast` must hold the same
values of it as `newdata` does, and a contrast that does not is refused
rather than overwritten.

A new `newdata` on a difference curve needs a new `contrast` with it,
and a call that gives one grid and not the other is refused. Without the
refusal the search would run on the FIRST curve alone, for an object
whose every row is a difference.

## Past a `ps()` knot span

This function REFUSES rather than warns. A
[`frmtmb::ps()`](https://aforren1.github.io/frmtmb/reference/ps.html)
basis decays to zero past its frozen knot span, so a curve drawn there
still has peaks and still crosses levels, and a root found among them is
a root of the decaying partial sum. Unlike a band, which shows the
reader what it is doing, that root leaves the function as a number with
a standard error beside it and nothing to say which curve it came off.
The bracket is checked at the grid scan, before any root is refined.

What is checked is the grid you passed, not the difference stencil the
scan widens it into: a grid whose endpoint sits exactly on a knot is
inside the span, and was refused for being a millionth of its range
outside the stencil's.

The search holds every column but `var` at row 1's value, so a grid
whose other columns change from row to row is checked a second time,
against the whole grid. That is where a second `ps()` term can leave its
span in a row the search itself never predicts at.

A difference curve has two grids and the search pins row 1 of each, so
the second check asks the question of each of them. It is not asked
whether the two grids differ from ONE ANOTHER: they always do, because
that is what a difference is.

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
#>   covariance checked against predict(se.fit = TRUE) to 2.22e-16 relative
#>   .feature .var .estimate       .se .lower_ci .upper_ci   .value .value_se
#> 1  maximum    x 0.4911416 0.0163341 0.4591274 0.5231559 2.053619 0.0409223
```
