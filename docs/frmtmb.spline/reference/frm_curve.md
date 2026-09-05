# A fitted curve on a grid, with pointwise and simultaneous bands

Evaluates a fitted linear predictor on a grid of covariate values and
returns it with two intervals: the usual pointwise interval, and a
SIMULTANEOUS band that covers the whole curve at once.

## Usage

``` r
frm_curve(
  object,
  newdata,
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

## The one internal this reaches into

The covariance is read from `fit$cache$Vjoint`, and that is a read into
frmtmb's INTERNALS rather than a sanctioned seam. The function that
writes it, `get_joint_cov()`, is `@noRd`; neither `fit$cache` nor the
`list(V =, names =)` shape of the memo appears in
`?frmtmb::`frmtmb-extension-api“. Unlike `fit$obj` and
\`fit\$estimates\`, which other extensions already read, this one has no
precedent to point at. It is stated here rather than buried in a
development note because a user is entitled to know which of a package's
dependencies are contractual and which are not.

It is made anyway, because every alternative was worse. Computing the
covariance here instead means a second `sdreport()` per call, a dense
Schur complement over coefficients the curve never touches (114 s and
2.1 GB at 8000 random coefficients, against 8.2 s and 1.2 GB for the
cache read), and a covariance that goes round `autoscale_sdreport()` and
is therefore WRONG on an autoscaled fit. A public-but-wrong route was
traded for a private-but-correct one.

What happens if core changes it, in full:

- **The name or the shape changes.** The read returns `NULL`, the sparse
  fallback runs, and the answer is the same one about 16 times slower.
  No wrong number.

- **The meaning of `V` changes without the name changing.** The
  covariance check catches it: every call compares `sqrt(diag(Sigma))`
  against `predict(se.fit = TRUE)` and refuses above `tol`. No wrong
  number.

- **It is absent.** It is absent on every FRESH fit, because it is a
  memo rather than a slot. This is not a dependency on the cache being
  warm: the `predict(se.fit = TRUE)` check runs first and warms it,
  which is why that call is ordered ahead of the covariance and why the
  ordering is enforced by an argument rather than by a comment.

So the reach cannot produce a wrong answer; it can only become slow, or
refuse. The standing ask is an exported accessor,
`dev/spline-seam-proposal.md` Part 1a, which would remove it.

## Cost

What this call costs is dominated by ONE thing, and it is not the
[`predict()`](https://rdrr.io/r/stats/predict.html) call count that the
`"check"` attribute reports. It is the single `predict(se.fit = TRUE)`
call, inside which core inverts the fit's joint precision matrix over
EVERY coefficient, including the ones this curve does not touch.
Measured at `re.form = NA` on a 20-point grid, one process each:

- `s(x, k = 10)`, 8 random coefficients: design rebuild 0.01 s,
  `predict(se.fit = TRUE)` 0.29 s.

- `s(t, k = 8) + (1 + t | subject)`, 1000 subjects and 2006 random
  coefficients: 0.07 s against 0.98 s.

- the same over 4000 subjects, 8006 random coefficients: 0.28 s against
  6.87 s.

The design rebuild is a tenth of the cost at every size, and the
covariance itself is FREE: it is read from the object
`predict(se.fit = TRUE)` has already cached, not recomputed. Size a job
from the joint-precision solve, which grows with the total number of
coefficients in the fit, and not from the grid or the call count.

The design is rebuilt with one
[`predict()`](https://rdrr.io/r/stats/predict.html) call per
contributing coefficient, plus one probe per block of `24` that
contributes nothing, so the count does not grow with the number of
LEVELS of a grouping factor. On the model
[`vignette("curve-inference")`](https://aforren1.github.io/frmtmb/frmtmb.spline/articles/curve-inference.md)
fits, `v ~ s(t, k = 12) + s(t, subject, bs = "fs", k = 5)` over 20
subjects, that is 32 calls against 110 random coefficients; a
factor-smooth model with NO population smooth needs far fewer, because
at `re.form = NA` the `fs` term contributes nothing and the population
curve is a constant. `tests/testthat/test-curve.R` pins both. The call
count is in the `"check"` attribute.

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
#>   covariance checked against predict(se.fit = TRUE) to NULL relative, in  predict() calls
#>            x .estimate        .se  .lower_ci  .lower_sim
#> 1 0.00000000 0.1320147 0.15875291 -0.1791353 -0.33430711
#> 2 0.04166667 0.3413067 0.11074391  0.1242526  0.01600685
#> 3 0.08333333 0.5478148 0.07848158  0.3939938  0.31728253
#> 4 0.12500000 0.7472177 0.06874628  0.6124775  0.54528198
#> 5 0.16666667 0.9386514 0.06827917  0.8048267  0.73808771
#> 6 0.20833333 1.1261965 0.06628101  0.9962881  0.93150224
```
