# A fitted curve on a grid, with pointwise and simultaneous bands

Evaluates a fitted linear predictor on a grid of covariate values and
returns it with two intervals: the usual pointwise interval, and a
SIMULTANEOUS band that covers the whole curve at once.

## Usage

``` r
frm_curve(
  object,
  newdata,
  contrast = NULL,
  dpar = NULL,
  resp = NULL,
  re.form = NA,
  level = 0.95,
  simultaneous = TRUE,
  nsim = 10000L,
  transform = FALSE,
  seed = NULL,
  tol = 1e-06
)
```

## Arguments

- object:

  A `frmtmb_fit` from
  [`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html).

- newdata:

  The grid, as a data frame. Every variable the linear predictor reads
  must be a column, held at the value the curve is wanted at.

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

- simultaneous:

  Compute the simultaneous band. `FALSE` returns the pointwise interval
  alone and skips the simulation.

- nsim:

  Draws in the max-deviation simulation. The default 10000 puts the
  Monte Carlo error of the critical value near 0.013; 200000 puts it
  near 0.003.

- transform:

  Return the curve and both bands through the link inverse. The bands
  are transformed end to end rather than rebuilt, which keeps their
  coverage under any monotone link.

- seed:

  Seed for the simulation, for a reproducible band.

- tol:

  Largest relative disagreement with `predict(se.fit = TRUE)` the
  assembled covariance may show before the call refuses.

## Value

A data frame of class `frmtmb_curve`: the columns of `newdata`, then
`.estimate`, `.se`, `.crit`, `.lower_ci`, `.upper_ci`, and when
`simultaneous = TRUE` also `.crit_sim`, `.lower_sim` and `.upper_sim`.
The grid covariance is the `"Sigma"` attribute, the fit is the `"fit"`
attribute, and `"check"` carries the covariance agreement and the
[`predict()`](https://rdrr.io/r/stats/predict.html) call count.

## Details

A pointwise interval is the wrong tool for the question a curve usually
raises. "Is the velocity above zero at 300 ms" is pointwise; "does this
curve have the shape I claim" is a statement about every point at once,
and a 95 percent pointwise band covers the whole curve far less than 95
percent of the time. The simultaneous band is the max-deviation
simulation of Ruppert, Wand and Carroll (2003, ch. 6): draw the curve's
own deviation process from its joint covariance, standardize each draw
by the pointwise standard error, take the largest absolute value over
the grid, and use the `level` quantile of those maxima in place of
`qnorm(0.975)`. It is the construction
`gratia::confint(type = "simultaneous")` uses on an mgcv fit, and this
package's simulation reproduces gratia's critical value inside its Monte
Carlo error.

## What the covariance is, and how it is checked

A penalized smooth's wiggly part is a random-effect block in the fitted
objective even when the smooth is a population term, so the covariance
of a curve needs the joint covariance of the fixed AND random
coefficients. frmtmb exports no route to it: `vcov(full = TRUE)` returns
the outer parameter vector, which excludes `b` under both of its
branches, and `predict(se.fit = TRUE)` forms the grid covariance
internally and returns only its diagonal.

So this function rebuilds it. The linear predictor is LINEAR in the
coefficients, so the difference between a prediction and the same
prediction with one coefficient raised by one is that coefficient's
design column, exactly. The joint covariance comes from the fit's own
joint precision matrix.

Neither piece was handed over by an exported function, so neither is
trusted. Every call recomputes `sqrt(diag(Sigma))` and compares it with
`predict(se.fit = TRUE)`, and refuses when the two disagree by more than
`tol`. The measured agreement is in the `"check"` attribute and is
reported by [`print()`](https://rdrr.io/r/base/print.html). On the
package's own test models it is at the tenth significant figure or
better.

## The route to the covariance

Both halves come from frmtmb's own exported seam,
[`frmtmb::frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.html),
which returns the design `A` over the coefficient vector, the joint
covariance `V` at exactly the rows `A`'s columns sit at, and the
variance that is NOT coefficient uncertainty (a new grouping level's
marginal variance, an exact `gp()`'s kriging variance) as a separate
element. `Sigma` is `A V A'`.

Up to frmtmb 0.51.0 there was no such seam. This package rebuilt `A` by
unit perturbation, one
[`predict()`](https://rdrr.io/r/stats/predict.html) call per
contributing coefficient, and read `V` out of `fit$cache$Vjoint`, which
was an internal with no precedent to point at. Both are gone: the
reconstruction and the reach were replaced by one call, this package now
requires frmtmb (\>= 0.52.0), and what the covariance check verifies has
changed from "the reconstruction reproduced core's number" to "the seam
is being read correctly".

The check itself stays. Every call recomputes `sqrt(diag(Sigma))` and
compares it with `predict(se.fit = TRUE)`, and refuses when the two
disagree by more than `tol`. The measured agreement is in the `"check"`
attribute and is reported by
[`print()`](https://rdrr.io/r/base/print.html).

The one case with nothing to check against is a nonlinear (`nl = TRUE`)
body: `predict(se.fit = TRUE)` is refused there, so
[`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.html)
is the only route to the number and `cov_rel_error` is `NA`.
[`print()`](https://rdrr.io/r/base/print.html) says so rather than
reporting a check that never ran.

## Cost

What this call costs is dominated by ONE thing: the single
`predict(se.fit = TRUE)` check call, inside which core inverts the fit's
joint precision matrix over EVERY coefficient, including the ones this
curve does not touch. Measured at `re.form = NA` on a 20-point grid, one
process each:

- `s(x, k = 10)`, 8 random coefficients: 0.29 s.

- `s(t, k = 8) + (1 + t | subject)`, 1000 subjects and 2006 random
  coefficients: 0.98 s.

- the same over 4000 subjects, 8006 random coefficients: 6.87 s.

The design rebuild that used to sit beside those figures, and that was a
tenth of them at every size, is gone:
[`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.html)
returns the design core already had, so the
[`predict()`](https://rdrr.io/r/stats/predict.html) call count no longer
depends on the number of coefficients at all. The joint-precision solve
is now the whole cost, it is paid once because core memoizes it, and it
grows with the total number of coefficients in the fit rather than with
the grid.

## A difference curve

`contrast` is a second grid of the same height. The curve returned is
then the DIFFERENCE of the two linear predictors, row by row, and its
covariance is `(A1 - A2) V (A1 - A2)'`, so both bands describe the
difference and the simultaneous one answers "is this difference anywhere
other than zero" over the whole grid at once. It is the quantity
`gratia::difference_smooths(group_means = TRUE)` reports for a factor-by
smooth, and on the same mgcv fit the two agree. Name that argument when
you compare: gratia's default, `group_means = FALSE`, zeroes the
intercept and the parametric group columns and reports the smooth-only
difference, which is a different quantity. On the fixture
`test-difference.R` uses it is 0.52 away.

A difference is NOT the difference of two calls to this function. The
two curves share coefficients, so their covariance is what the
difference is made of, and adding two standard errors in quadrature
would ignore it.

What the difference path cannot do, and refuses by name:

- `transform = TRUE`. A difference of linear predictors is not a
  difference of responses under any link but the identity, so there is
  nothing to map it through.

- Two grids that load on different coefficients.

- Two grids that load DIFFERENT draws of a latent field whose variance
  is not coefficient uncertainty. See the next section.

The covariance check also means less here, and
[`print()`](https://rdrr.io/r/base/print.html) says so.
`predict(se.fit = TRUE)` returns a marginal standard error per row and
never the covariance between the grids, so the check runs on each half
and `cov_rel_error` is the worse of the two: what it licenses is that
both designs were read correctly.

## An exact `gp()` under a difference

An exact `gp()` evaluated off the observed positions carries a kriging
residual that is not coefficient uncertainty.
[`frmtmb::frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.html)
returns its variance one number per ROW and returns no covariance
BETWEEN two grids, so `var(g1 - g2)` has no public route.

It needs none in the ordinary case. A contrast taken across a factor at
ONE `gp()` position leaves both grids loading the same residual, so it
cancels exactly and the difference is `(A1 - A2) V (A1 - A2)'` with
nothing left over. That case is computed rather than refused.

Sameness is decided on the design and not on the numbers: the two grids
must agree bit for bit on EVERY column of `A` outside the fixed effects.
Equality of the variances would not be enough, because two levels of one
grouping block have identical marginal variances by construction and are
different draws.

That test is stricter than the mathematics needs, and it is worth
knowing where the extra strictness bites. Only the block carrying the
kriging residual has to match for the residual to cancel, but the test
asks it of every latent column, so a contrast across `fac` on
`y ~ fac + s(x, by = fac) + gp(x)` is REFUSED even though the `gp()`
columns are identical: the by-factor smooth's own columns differ, which
is what a by-factor smooth is for. It fails closed, so the cost is an
answer you do not get rather than one you should not trust.

The rest is refused because the SEAM cannot supply it, not because the
mathematics is missing. The conditional cross-covariance is
`k(x1, x2) - Xr1 K Xr2'`, and core forms every piece of it while
predicting, but reduces the result to one variance per row before the
seam returns. Hold every latent term equal between the grids and
contrast a fixed effect, or read the two curves separately.

## Past a [`ps()`](https://aforren1.github.io/frmtmb/reference/ps.html) knot span

A [`frmtmb::ps()`](https://aforren1.github.io/frmtmb/reference/ps.html)
basis is a partition of unity only between its frozen outer knots. Past
them it is a partial sum that decays to zero, so a curve drawn there
bends smoothly to whatever the rest of the body gives, which is exactly
the shape a reader does not question. `predict(newdata = )` says so, and
so does
[`frmtmb::frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.html),
the seam this function reads. It is surfaced again here, ONCE per
[`ps()`](https://aforren1.github.io/frmtmb/reference/ps.html) term per
call and carrying the span, so that the sentence names the function you
called: this one reads the seam on the grid, but
[`frm_curve_deriv()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_deriv.md)
reads it on a three-point difference stencil and
[`frm_curve_feature()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_feature.md)
on a five-point one, and core counts the rows it was handed.

[`frm_curve_feature()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_feature.md)
REFUSES instead of warning. A band past the span is visibly wrong on the
page; a peak located past it is a number with a standard error beside it
and nothing to give it away.

## References

Ruppert, D., Wand, M. P. and Carroll, R. J. (2003) *Semiparametric
Regression*. Cambridge University Press, ch. 6.

## See also

[`frm_curve_deriv()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_deriv.md)
for the derivative of the same curve,
[`frm_curve_feature()`](https://aforren1.github.io/frmtmb/frmtmb.spline/reference/frm_curve_feature.md)
for the location of a peak or a crossing.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200)))
dd$y <- 2 * sin(pi * dd$x) + rnorm(200, 0, 0.4)
fit <- frmtmb::frm(frmtmb::bf(y ~ s(x, k = 8)),
                   family = stats::gaussian(), data = dd)
cv <- frm_curve(fit, newdata = data.frame(x = seq(0, 1, length.out = 25)),
                nsim = 2000)
head(cv[, c("x", ".estimate", ".se", ".lower_ci", ".lower_sim")])
#> <frmtmb curve> , 6 grid points, level 
#>   critical value: pointwise NULL
#>   covariance NOT checked: predict(se.fit = TRUE) is refused for a nonlinear predictor, so there is no second route to compare against
#>            x .estimate        .se  .lower_ci  .lower_sim
#> 1 0.00000000 0.1320147 0.15875291 -0.1791353 -0.33430711
#> 2 0.04166667 0.3413067 0.11074391  0.1242526  0.01600684
#> 3 0.08333333 0.5478148 0.07848158  0.3939938  0.31728253
#> 4 0.12500000 0.7472177 0.06874628  0.6124775  0.54528198
#> 5 0.16666667 0.9386514 0.06827917  0.8048267  0.73808771
#> 6 0.20833333 1.1261965 0.06628101  0.9962881  0.93150224
```
