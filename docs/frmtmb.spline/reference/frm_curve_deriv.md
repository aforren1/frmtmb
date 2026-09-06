# Derivatives of a fitted curve

The first or second derivative of a fitted linear predictor with respect
to one covariate, on a grid, with delta-method standard errors and the
same pair of intervals
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
gives: pointwise, and simultaneous over the whole grid.

## Usage

``` r
frm_curve_deriv(
  object,
  var,
  order = 1L,
  newdata = NULL,
  dpar = NULL,
  resp = NULL,
  re.form = NA,
  level = 0.95,
  simultaneous = TRUE,
  nsim = 10000L,
  eps = NULL,
  seed = NULL,
  tol = 1e-06
)
```

## Arguments

- object:

  A `frmtmb_fit`, or a `frmtmb_curve` from
  [`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md),
  in which case its grid and its `dpar`, `resp` and `re.form` are
  reused.

- var:

  Name of the covariate to differentiate with respect to. It must be a
  numeric column of the grid.

- order:

  1 or 2.

- newdata:

  The grid. Required when `object` is a fit; taken from the curve
  otherwise.

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

- simultaneous:

  Compute the simultaneous band. `FALSE` returns the pointwise interval
  alone and skips the simulation.

- nsim:

  Draws in the max-deviation simulation. The default 10000 puts the
  Monte Carlo error of the critical value near 0.013; 200000 puts it
  near 0.003.

- eps:

  Step size. `NULL`, the default, is `1e-6` of the grid's range at
  `order = 1` and `1e-4` of it at `order = 2`.

- seed:

  Seed for the simulation, for a reproducible band.

- tol:

  Largest relative disagreement with `predict(se.fit = TRUE)` the
  assembled covariance may show before the call refuses.

## Value

A `frmtmb_curve` data frame, as
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)
returns, whose `.estimate` is the derivative.

## Details

The derivative is taken by central differences of the DESIGN, not of the
fitted values, so the estimate and its standard error describe one
function rather than two. Writing `D` for the differenced design and `V`
for the joint coefficient covariance, the reported curve is `D c` and
its covariance is `D V D'`, which is exact for the differenced basis;
the only approximation is the difference itself, and `eps` controls
that.

## Past a [`ps()`](https://aforren1.github.io/frmtmb/reference/ps.html) knot span

Warned once per
[`frmtmb::ps()`](https://aforren1.github.io/frmtmb/reference/ps.html)
term per call, as
[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md)'s
section describes, and counted over the grid you passed rather than over
the three-point stencil the design is built on.

Past the outer knot the basis is exactly zero, so the DERIVATIVE design
is exactly zero and those rows carry a standard error of exactly zero.
They are covered with probability one, so they leave the max-deviation
simulation instead of standardizing a zero deviation by a zero divisor;
`simultaneous = TRUE` therefore still returns a band, computed over the
rows that carry uncertainty. A grid on which EVERY row is past the outer
knot has no such row, and the call refuses rather than returning a band
with no content.

## See also

[`frm_curve()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve.md),
[`frm_curve_feature()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_feature.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200)))
dd$y <- 2 * sin(pi * dd$x) + rnorm(200, 0, 0.4)
fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
                   family = stats::gaussian(), data = dd)
g <- data.frame(x = seq(0.05, 0.95, length.out = 19))
d1 <- frm_curve_deriv(fit, var = "x", order = 1, newdata = g,
                      simultaneous = FALSE)
# the derivative of 2 sin(pi x) is 2 pi cos(pi x)
head(cbind(g, fitted = d1$.estimate, truth = 2 * pi * cos(pi * g$x)))
#>      x   fitted    truth
#> 1 0.05 4.992535 6.205829
#> 2 0.10 4.808133 5.975664
#> 3 0.15 4.574524 5.598359
#> 4 0.20 4.494088 5.083204
#> 5 0.25 4.589734 4.442883
#> 6 0.30 4.605194 3.693164
```
