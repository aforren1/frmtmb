# Noise-free predictors with known measurement error

`me(x, sdx)` in a formula says that the column `x` is a noisy
measurement of a latent, noise-free value, with known measurement
standard deviation `sdx`, and that the model uses the latent value. It
is brms's `me()` term with brms's meaning.

## Details

The model is \$\$x_i \sim N(\tilde{x}\_i, sdx_i), \qquad \tilde{x}\_i
\sim N(meanme, sdme)\$\$

where `meanme` and `sdme` are estimated. The coefficient of the term is
`bsp_me<x><sdx>` in brms's names, for example `bsp_mexsx` for
`me(x, sx)`.

`me()` is not a function. It is a formula special that the parser reads,
like `mo()` and `mi()`.

## Arguments of `me()`

- `x`: the noisy variable, a numeric column or an expression of columns,
  such as `log(x)`.

- `sdx`: the known measurement SD, a column, an expression or a positive
  constant.

- `gr`: optional. The name of a grouping column. There is then one
  latent value per level of `gr`, and `x` and `sdx` must be constant
  within each level.

## Estimation

The latent values are integrated out by the Laplace approximation,
together with the random effects. For a gaussian response with a linear
`me()` term the marginal likelihood is multivariate normal, and the
approximation is exact. An interaction of two `me()` terms, or a
non-gaussian response, makes it approximate.

Several `me()` terms that share a grouping have correlated latent values
by default, as in brms. Use
[`set_mecor()`](https://aforren1.github.io/frmtmb/reference/set_mecor.md)
to make them independent.

`me()` can appear on its own, in a `:` or `*` interaction with numeric
variables or with other `me()` terms, and in the formula of any
distributional or nonlinear parameter. It cannot appear in a group-level
term, in a nonlinear formula body, inside another function call, or in
the same interaction as `mo()` or `mi()`.

## Parameters and methods

The hyperparameters are the template components `meanme`, `logsdme` (the
log of `sdme`) and `thetame` (the correlation parameters). The
Noise-free Terms section of
[`summary()`](https://rdrr.io/r/base/summary.html) reports them on
brms's scale as `meanme_<coef>`, `sdme_<coef>` and
`corme__<coef1>__<coef2>`. Priors on brms's classes `meanme`, `sdme` and
`corme` are refused; class `"b"` addresses the `me()` coefficients by
their names, for example `coef = "mexsx"`.

[`anova()`](https://rdrr.io/r/stats/anova.html) compares only fits with
the same `me()` calls, because each call puts its noisy variable's
measurements into the likelihood.
[`AIC()`](https://rdrr.io/r/stats/AIC.html) does not check this, so
compare by AIC only across fits that share their `me()` calls.

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) on the fitted data
use the estimated latent values (the conditional modes). On new data the
observed value of the noisy variable stands in for the latent value.
brms draws that value from `N(x, sdx)`, so the two agree on the mean of
a linear term, and frmtmb does not add the measurement noise to a
prediction interval.

## See also

[`set_mecor()`](https://aforren1.github.io/frmtmb/reference/set_mecor.md);
`mi()` for a measurement model that has covariates of its own,
`bf(x | mi(sdx) ~ z)`.

## Examples

``` r
set.seed(1)
n <- 200
tx <- rnorm(n, 1, 0.8)
d <- data.frame(x = tx + rnorm(n, 0, 0.4), sx = 0.4)
d$y <- 2 + 0.7 * tx + rnorm(n, 0, 0.5)

fit <- frm(bf(y ~ me(x, sx)) + gaussian(), data = d)
fixef(fit)
#>            Estimate  Est.Error     Q2.5     Q97.5
#> Intercept 1.9213924 0.08471428 1.755355 2.0874294
#> mexsx     0.7443341 0.06947988 0.608156 0.8805121

# the naive slope is attenuated by the measurement error
coef(lm(y ~ x, data = d))[["x"]]
#> [1] 0.5741128
```
