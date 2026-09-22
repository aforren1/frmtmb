# Fixed effects per linear predictor

The coefficients of each linear predictor, as a named list keyed by
distributional parameter (and by response on a multivariate fit), each
entry named by DESIGN COLUMN: `fixef_by_dpar(fit)$sigma`,
`fixef_by_dpar(fit)$mu[["(Intercept)"]]`.

## Usage

``` r
fixef_by_dpar(object)
```

## Arguments

- object:

  A `frmtmb_fit`.

## Value

A named list of coefficient vectors, one per linear predictor.

## Details

This is the shape
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md)
returned before frmtmb took brms's summary matrix. It is kept because
nothing else has it:
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) names
its rows as brms does and flattens every predictor into one table, and
`fixef(flatten = TRUE)` is one vector in the
[`confint()`](https://rdrr.io/r/stats/confint.html) spelling. A model
with several linear predictors (a mixture, a multivariate response, a
nonlinear body) is the case this reads cleanly.

It carries ESTIMATES only. For a standard error or an interval use
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md), whose
rows are the same coefficients under brms's names.

## See also

[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
[`coef.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/coef.frmtmb_fit.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100))
dd$y <- rnorm(100, 1 + 0.5 * dd$x, exp(0.2 + 0.1 * dd$x))
fit <- frm(bf(y ~ x, sigma ~ x) + gaussian(), data = dd)

# one entry per distributional parameter, named by design column
fixef_by_dpar(fit)
#> $mu
#> (Intercept)           x 
#>   0.9540016   0.5184273 
#> 
#> $sigma
#> (Intercept)           x 
#>   0.1515397   0.1026294 
#> 
exp(fixef_by_dpar(fit)$sigma[["(Intercept)"]])
#> [1] 1.163624

# the same coefficients under brms's names, with their errors
fixef(fit)
#>                  Estimate  Est.Error        Q2.5     Q97.5
#> Intercept       0.9540016 0.11693512  0.72481294 1.1831902
#> sigma_Intercept 0.1515397 0.07154025  0.01132334 0.2917560
#> x               0.5184273 0.13753420  0.24886519 0.7879893
#> sigma_x         0.1026294 0.09977117 -0.09291848 0.2981773
```
