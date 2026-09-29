# emmeans support

frmtmb registers
[`emmeans::emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.html)
support for `frmtmb_fit` objects. The arguments that select what is
averaged are brms's: `dpar`, `nlpar`, `resp`, `epred` and `re_formula`,
with brms's defaults. Give them to `emmeans()` or `ref_grid()` directly,
and emmeans passes them on.

## Value

This page documents the emmeans methods and returns nothing. `emmeans()`
on a `frmtmb_fit` returns an `emmGrid` object.

## What is averaged

- default:

  The linear predictor of `mu`, on the link scale, at the population
  level (`re_formula = NA`).

- `dpar = "sigma"`:

  The linear predictor of that distributional parameter, on its own link
  scale.

- `nlpar = "a"`:

  The linear predictor of that nonlinear parameter. It is linear in its
  own coefficients, so the basis is the usual design, coefficients and
  covariance.

- `resp = "y1"`:

  One response of a multivariate fit.

- no `resp` on a multivariate fit:

  Every response, stacked as the multivariate factor `rep.meas` whose
  levels are the response names. This is brms's layout. Average over it,
  condition on it (`by = "rep.meas"`) or compare its levels with
  `contrast()`.

- `epred = TRUE`:

  The expected value of the response, as
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) reports it,
  on the response scale.

- `re_formula = NULL`:

  Keeps the group-level effects. The grouping factors join the reference
  grid, so the marginal mean averages over their observed levels, as in
  brms.

## How the uncertainty is computed

brms computes each marginal mean from the posterior draws. frmtmb uses
the Wald covariance of the estimates instead. For a linear predictor
with only parametric terms, the basis is emmeans's usual one: the design
at the reference grid, the coefficients and their covariance. Otherwise
the basis has brms's shape: one column per grid point, the predictions
at the grid as the estimates, and the covariance `J V J'`, where `V` is
the joint covariance of the estimates and `J` is the Jacobian of the
grid predictions that
[`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.md)
computes. That route is used for a nonlinear predictor, for
`epred = TRUE`, for `re_formula` other than `NA`, and for a predictor
with a `s()`, `t2()`, `gp()`, `mo()` or `mi()` term. For a linear
predictor it is exact. For a nonlinear predictor or the expected
response it is the delta method.

## Scales and averaging

emmeans averages the grid estimates with linear weights. For the default
route that is the average of linear predictors, as for any generalized
linear model. With `epred = TRUE` it is the average of the predicted
means, which is how brms marginalizes the expected response.
`type = "response"` applies the inverse link of the selected predictor
once, after the average. With `epred = TRUE` there is no link to apply,
because the estimates are already on the response scale, so
`type = "response"` changes nothing. A nonlinear `mu` is reported on its
link scale, where the body's value lives; `type = "response"` applies
the family's inverse link to it. A response transformation such as
`log(y)` is detected for `mu` and for `epred = TRUE` only, since it does
not apply to a distributional or nonlinear parameter.

An ordinal fit is averaged on its latent scale, which is emmeans's
`mode = "latent"` convention for `clm`-like models. The thresholds are
not in the basis, so the means carry no threshold offset and
`type = "response"` does not transform them. Contrasts are not affected
by the offset.

## Refusals

A refusal names its reason. emmeans hides an error that `recover_data()`
raises behind "Perhaps a 'data' or 'params' argument is needed", so
frmtmb returns the reason as text instead, which emmeans reports as the
error. Refused:

- `dpar` and `nlpar` together, and a name the model does not have (brms
  refuses both).

- `epred = TRUE` for an ordinal or categorical family, whose expected
  response is a distribution over categories. Use the latent predictor,
  or
  [`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md)
  with `type = "response"`.

- No `resp` on a multivariate fit whose responses have different links,
  since the stacked predictors would be on different scales. brms
  refuses this too. Name one response, or use `epred = TRUE`.

- An exact `gp()` term predicted at a position the fit did not see. Its
  kriging variance has no covariance between two grid points here. Use
  an approximate `gp(..., k = )`, or hold the covariate at an observed
  value with the `at` argument.

## Divergence from brms

For a nonlinear `mu` without `nlpar`, the brms reference grid holds only
the covariates of the nonlinear body. frmtmb also puts in the covariates
of the nonlinear parameters' own formulas, so that emmeans averages over
them rather than holding them at an unstated value.

## See also

[`frm_lp_basis()`](https://aforren1.github.io/frmtmb/reference/frm_lp_basis.md)
for the Jacobian behind the delta method, and
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) for the
expected response.

## Examples

``` r
# \donttest{
if (requireNamespace("emmeans", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(f = factor(rep(c("a", "b", "c"), 40)),
                  x = runif(120))
  d$y <- 1 + as.numeric(d$f) + 2 * d$x + rnorm(120, 0, 0.3)
  d$z <- as.numeric(d$f) + rnorm(120)

  # a nonlinear model: one of its parameters, then the whole mean
  fnl <- frm(bf(y ~ a + b * x, a ~ f, b ~ 1, nl = TRUE), data = d)
  emmeans::emmeans(fnl, "f", nlpar = "a")
  emmeans::emmeans(fnl, "f")
  pairs(emmeans::emmeans(fnl, "f", epred = TRUE))

  # a multivariate model: one response, then both as rep.meas
  fmv <- frm(mvbf(bf(y ~ f), bf(z ~ f)), data = d)
  emmeans::emmeans(fmv, "f", resp = "y")
  emmeans::emmeans(fmv, ~ f | rep.meas)
}
#> rep.meas = y:
#>  f emmean     SE  df asymp.LCL asymp.UCL
#>  a  3.105 0.0938 Inf     2.921      3.29
#>  b  4.000 0.0938 Inf     3.816      4.18
#>  c  5.004 0.0938 Inf     4.820      5.19
#> 
#> rep.meas = z:
#>  f emmean     SE  df asymp.LCL asymp.UCL
#>  a  0.964 0.1590 Inf     0.653      1.27
#>  b  2.130 0.1590 Inf     1.819      2.44
#>  c  2.899 0.1590 Inf     2.588      3.21
#> 
#> Confidence level used: 0.95 
# }
```
