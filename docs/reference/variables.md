# Usable parameter names

The names that
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
expressions (and
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
targeting) accept: fixed-effect coefficients under their
[`vcov()`](https://rdrr.io/r/stats/vcov.html) names with parentheses
stripped, natural-scale random-effect summaries (`sd_<group>__<term>`,
`cor_<group>__<t1>__<t2>`), and `sigma` when the residual SD is a
scalar. The brms spelling; for sampled fits, `variables()` on the
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
result lists the draw columns instead.

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

  Unused.

## Value

A character vector.

## Details

A residual correlation term
([frmtmb-autocor](https://aforren1.github.io/frmtmb/reference/frmtmb-autocor.md))
contributes its natural-scale parameters under brms's names, sanitized
the same way: `ar1`, `ar2`, `ma1`, `cosy`, `cortime__<t1>__<t2>`.

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
#> [1] "Intercept"       "x"               "sigma_Intercept" "sd_g__Intercept"
#> [5] "sigma"          
```
