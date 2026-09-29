# Usable parameter names

brms's names for the parameters of a fit, which are the names
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
expressions accept: fixed-effect coefficients with brms's `b_` prefix
(`b_Intercept`, `b_x`, and `b_sigma_Intercept` for a coefficient of a
`sigma` formula), natural-scale random-effect summaries
(`sd_<group>__<coef>`, `cor_<group>__<coef1>__<coef2>`, with a
distributional or nonlinear parameter's name in the coefficient part,
`sd_g__sigma_Intercept`), and a distributional parameter nobody wrote a
formula for, on its natural scale (`sigma`, `shape`, `sigma_ya`). Every
name is spelled through brms's renaming: `b_IxE2` for `I(x^2)`,
`sd_g:h__Intercept` for `(1 | g:h)`. For sampled fits, `variables()` on
the
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
result lists the draw columns, which follow the same convention EXCEPT
for an ordinal fit, where they are the internal names: `tau_raw_1`,
`tau_raw_2` for the thresholds this page calls `b_Intercept[1]`,
`b_Intercept[2]`, and `bcs2_1` for a `cs()` coefficient this page calls
`bcs_<column>[1]`.
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) on the
draws object does report brms's rows.

## Usage

``` r
variables(x, ...)

# S3 method for class 'frmtmb_fit'
variables(x, ...)
```

## Arguments

- x:

  A `frmtmb_fit` or `frmtmb_draws`.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

## Value

A character vector.

## Details

brms's `variables()` also lists what a fit has no counterpart of:
group-level coefficients `r_<group>[<level>,<coef>]`, the centered
`Intercept`, `lprior` and `lp__`. A maximum-likelihood fit has no draws
of those, and
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md)
reports the conditional modes.

An ordinal fit's thresholds are `b_Intercept[k]`, and a `cs()`
category-specific term contributes `bcs_<column>[k]`, one name per
DESIGN COLUMN per category boundary. The column is the one
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
builds, as it is in brms, so `cs(x)` on a numeric predictor gives
`bcs_x[1]`, `bcs_x[2]`, while `cs(f)` on a factor with levels `a`, `b`
and `c` gives the treatment-contrast pairs `bcs_fb[k]` and `bcs_fc[k]`.
See the "Category-specific effects, cs()" section of
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md).

A residual correlation term
([frmtmb-autocor](https://aforren1.github.io/frmtmb/reference/frmtmb-autocor.md))
contributes its natural-scale parameters under brms's names: `ar[1]`,
`ar[2]`, `ma[1]`, `cosy`, `cortime__<t1>__<t2>`.

`gr(cov = )`, `gr(prec = )` and `equalto()` blocks contribute
`sd_`/`cor_` names for their within-level covariance. Smooths,
`gp()`/`hsgp()`, `car()` and `spde()` blocks contribute none: their
parameters are not standard deviations. See the "Which random-effect
blocks contribute names" section of
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md).

## Examples

``` r
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
variables(fit)
#> [1] "b_Intercept"     "b_x"             "sd_g__Intercept" "sigma"          
```
