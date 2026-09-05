# A penalized spline whose value a nonlinear body consumes

`ps()` declares a penalized coefficient block inside the body of a
nonlinear formula (`bf(..., nl = TRUE)` or
[`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md)) and
evaluates to the spline's VALUE at `expr`. It is not added to a linear
predictor, which is what separates it from
[`mgcv::s()`](https://rdrr.io/pkg/mgcv/man/s.html): the body may
multiply it, exponentiate it, or pass it on.

## Usage

``` r
ps(expr, k = 10, degree = 3, pad = 0.1, center = TRUE)
```

## Arguments

- expr:

  The expression the spline is evaluated at. It may name nonlinear
  parameters, so the evaluation point moves with the fit -
  `ps(age + shift)` is a curve each subject sees shifted along its own
  time axis.

- k:

  Number of basis functions (default 10). At least `max(degree + 1, 4)`
  and at most 50; see Accuracy.

- degree:

  Degree of the B-spline pieces (default 3, cubic).

- pad:

  Fraction of the data-time range of `expr` added at each end before the
  knots are placed (default 0.1). The basis is exactly zero outside its
  knot span, so an evaluation point that the parameters push past the
  padded range reads as a curve value of zero rather than as an
  extrapolation. Widen `pad` when the transformation is large; the fit
  reports the coverage it achieved.

- center:

  If `TRUE` (default), the coefficients are constrained to sum to zero.
  B-splines are a partition of unity, so the constraint removes the
  curve's overall level, which is otherwise confounded with an intercept
  in the body.

## Value

`ps()` is a term, not a function to call. Evaluating one outside a
nonlinear body is an error.

## Where the pieces go

The second-difference penalty `S` is eigendecomposed and the basis is
split the way
[`mgcv::smooth2random()`](https://rdrr.io/pkg/mgcv/man/smooth2random.html)
splits `s()`: the null space of `S` joins the fixed coefficients
(`beta`, or `betad` for a distributional parameter), and the range space
becomes one random-effect block with a single variance in `theta`. The
smoothing parameter is that variance's inverse, so smoothness is
estimated jointly with every other variance component rather than
chosen.

[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
and [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md)
report the block under the term's own label. Its standard deviation is a
smoothing parameter, not a subject effect.

## The knots are frozen

The knot vector is built once, from the range of `expr` evaluated on the
model frame with every nonlinear parameter set to zero, padded by `pad`.
It is then frozen on the frame and reused on `newdata`, which is
[`poly()`](https://rdrr.io/r/stats/poly.html)'s rule and for the same
reason: a basis rebuilt from new data is a different basis, and the
coefficients would not mean what they were fitted to mean.

Boundary knots are DISTINCT and lie outside the padded range. The
divided-difference form degenerates on the repeated boundary knots
[`splines::splineDesign()`](https://rdrr.io/r/splines/splineDesign.html)
normally uses, so `ps()` places them apart rather than exposing the
choice.

This is NOT the rule of D'Alessandro, Thoresen and Sorensen (2026),
section 2.3.1, and `ps()` is not a reimplementation of it. That paper
rescales the spline argument to `[0, 1]` between bounds
`min(t) - 3 sigma_b` and `max(t) + 3 sigma_b`, which are smooth
functions of the CURRENT variance estimates, so its knots move with the
fit. `ps()` freezes instead, for two reasons. Moving knots make the
penalty matrix `S` and its eigendecomposition functions of the
parameters, so the split into fixed and random coefficients would have
to be redone inside the objective; and a basis that is rebuilt between
the fit and the prediction is a different basis, which is the rule
[`poly()`](https://rdrr.io/r/stats/poly.html) follows and the reason
`pad =` exists. On the paper's own application the two agree closely
enough that six of its seven reported parameters come back at the
published precision, but they are related constructions rather than one.

## Accuracy

The basis agrees with
[`splines::splineDesign()`](https://rdrr.io/r/splines/splineDesign.html)
to 3.1e-13 absolute at `k = 12` and 1.8e-11 at `k = 40`, with an AD
input giving the same values bit for bit and the taped derivative
matching `splineDesign(derivs = 1)` to 4.5e-13 and 4.5e-11 respectively.
The error grows like `k^(degree - 1)` because a divided difference of
truncated powers cancels terms of order `(range / spacing)^degree`. `k`
above 50 is refused for that reason; a curve that needs more than 50
basis functions needs a different basis, not a looser tolerance.

## What the rest of the package does with it

[`predict()`](https://rdrr.io/r/stats/predict.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) work: the body is
re-evaluated at the new data or at the drawn coefficients, and the
frozen basis is evaluated at whatever argument comes out.

`predict(se.fit = TRUE)` stays refused for a nonlinear predictor.
[`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.md)
is the route: it tapes the body and returns `d eta / d coef` as a
Jacobian, which is what a delta method over a warped curve needs.

`REML = TRUE`, `quadrature = TRUE`, `frmtmb_control(profile = TRUE)` and
multivariate
([`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md)) models
are refused by name. The first three integrate out something a `ps()`
block has already put in the Laplace approximation; the fourth is
untested.

The importance correction
([`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)'s
`importance =`) is refused, and not by anything `ps()` declares: core
refuses the correction for EVERY nonlinear predictor, because a body
mixes parameter values with raw data columns and the corrected objective
evaluates the predictor once per draw.

## References

D'Alessandro, M., Thoresen, M. and Sorensen, O. (2026). A Semiparametric
Nonlinear Mixed Effects Model with Penalized Splines Using Automatic
Differentiation. arXiv:2603.11728.

Wood, S. N. (2004). Stable and efficient multiple smoothing parameter
estimation for generalized additive models. Journal of the American
Statistical Association 99, 673-686.

## See also

[`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.md)
for standard errors on a curve a `ps()` block builds,
[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)

## Examples

``` r
set.seed(1)
n_id <- 40
d <- expand.grid(t = seq(0, 1, length.out = 8), id = factor(1:n_id))
sh <- stats::rnorm(n_id, 0, 0.05)[d$id]
d$y <- 2 + sin(2 * pi * (d$t + sh)) + stats::rnorm(nrow(d), 0, 0.1)
fit <- frm(bf(y ~ lev + ps(t + shift, k = 8),
              lev ~ 1, shift ~ 0 + (1 | id), nl = TRUE),
           data = d)
#> Warning: 2 of 320 fitted values of t + shift fall outside the knot span of ps(t + shift, k = 8) [-0.1, 1.1], where the basis is exactly zero and its gradient with it. Refit with a larger pad =
fixef(fit)
#> $lev
#> (Intercept) 
#>    2.002918 
#> 
#> $shift
#> numeric(0)
#> 
#> $mu
#> numeric(0)
#> 
#> $sigma
#> (Intercept) 
#>   -2.310741 
#> 
```
