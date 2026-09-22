# Which scale each method reports

Every quantity a fitted model reports lives on one of three scales, and
nothing in the output says which. This page says it once, per method, so
that a ported script can be read without guessing.

## Value

This page documents a convention. It is not a function, so it returns no
value.

## Details

The three scales are:

- **Link.** The linear predictor of one distributional parameter, before
  its inverse link is applied. This is where the coefficients live,
  where the standard errors are symmetric, and where a random effect has
  its variance.

- **Response.** The units of the outcome itself.

- **Unitless.** A standardized quantity, such as a Pearson or quantile
  residual, which has no units by construction.

## What each method returns

|  |  |  |
|----|----|----|
| method | returns | scale |
| [`predict()`](https://rdrr.io/r/stats/predict.html) | a summary of the PREDICTIVE distribution | response |
| [`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md), default | the `mu` linear predictor | link |
| `frm_linpred(type = "response")` | the conditional mean | response |
| `frm_linpred(dpar = "sigma")` | that predictor | link; `"response"` gives response |
| [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) | a summary of the conditional mean | response |
| [`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md) | `estimate__` and its band | response |
| [`residuals()`](https://rdrr.io/r/stats/residuals.html), default | observed minus [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) | response |
| `residuals(type = "pearson")` | that, over the conditional SD | unitless |
| `residuals(type = "deviance")` | signed root unit deviance | unitless |
| `residuals(type = "osa")` | one-step-ahead quantile residual | unitless |
| [`simulate()`](https://rdrr.io/r/stats/simulate.html) | draws of the outcome | response |
| [`coef()`](https://rdrr.io/r/stats/coef.html), [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md), [`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md) | coefficients, modes | link, per dpar |
| [`confint()`](https://rdrr.io/r/stats/confint.html), [`vcov()`](https://rdrr.io/r/stats/vcov.html), [`summary()`](https://rdrr.io/r/base/summary.html) | the same | link, per dpar |
| [`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md) | SDs and correlations of a predictor | link, per dpar |
| [`sigma()`](https://rdrr.io/r/stats/sigma.html) | the residual SD, `NA` if it varies by row | response |
| [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md) | `Estimate` and its interval | link, per dpar |
| [`posterior_summary()`](https://aforren1.github.io/frmtmb/reference/posterior_summary.md) | a summary of a draws MATRIX | that matrix's own |
| [`pp_check()`](https://aforren1.github.io/frmtmb/reference/pp_check.md) | the outcome against replicates | response |
| [`bayes_R2()`](https://aforren1.github.io/frmtmb/reference/bayes_R2.md) | refuses on a `frmtmb_fit` | neither |
| [`logLik()`](https://rdrr.io/r/stats/logLik.html), [`AIC()`](https://rdrr.io/r/stats/AIC.html), [`BIC()`](https://rdrr.io/r/stats/AIC.html) | the fitted likelihood | the data's own |

## Where this differs from brms

[`predict()`](https://rdrr.io/r/stats/predict.html) used to be the
first, and is no longer: since item 2.6d it is brms's, a summary of the
predictive distribution on the response scale. The LINEAR PREDICTOR,
which is what it returned before and what the glmmTMB `type` vocabulary
reaches, moved to
[`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md)
unchanged. On a lognormal fit the two differ by about 8 against about
6,700, which is why they could not share a name. Measured with
`dev/generics-scale.R` and `dev/generics-scale-brms.R`.

`sigma` is the one that remains.
[`summary()`](https://rdrr.io/r/base/summary.html)'s
`Regression Coefficients` block and
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) report
every coefficient on its own link, so a `sigma` WITH A FORMULA and the
default log link is reported as `log(sigma)` and can be negative. A
`sigma` nobody wrote a formula for is not a coefficient at all: it is in
[`summary()`](https://rdrr.io/r/base/summary.html)'s
`Further Distributional Parameters` block on its own response scale, as
in brms, and [`sigma()`](https://rdrr.io/r/stats/sigma.html) returns
that number.

## What agrees with brms without any conversion

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.md),
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md),
[`conditional_effects()`](https://aforren1.github.io/frmtmb/reference/conditional_effects.md)
and the coefficients of a linear predictor whose link is the identity.

For a lognormal fit with a CONSTANT sigma and no truncation,
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is
`exp(mu + sigma^2 / 2)`, the mean of the distribution rather than its
median `exp(mu)`; the ratio between the two is `exp(sigma^2 / 2)`, which
is the whole of the discrepancy a ported script sees. **Both conditions
are load-bearing**, and neither is obvious from the formula:

- With a distributional sigma (`bf(y ~ x, sigma ~ x)`),
  [`sigma()`](https://rdrr.io/r/stats/sigma.html) returns `NA` with a
  warning, because there is no one number to return, and the formula
  gives `NA` with it. Use
  `frm_linpred(dpar = "sigma", type = "response")`, which reproduces
  `fitted(dpar = "sigma")[, "Estimate"]` exactly
  ([`identical()`](https://rdrr.io/r/base/identical.html) is `TRUE`).

- Under truncation,
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) is the
  TRUNCATED mean and the formula answers a different question. With a
  lower bound `lb`, the truncated mean is
  `exp(mu + sigma^2 / 2) * pnorm(sigma - a) / pnorm(-a)` with
  `a = (log(lb) - mu) / sigma`, so the naive formula falls short by a
  relative `1 - pnorm(-a) / pnorm(sigma - a)` in each row. That is an
  identity, not an estimate, and there is no single number for it: it
  depends on where the bound sits in each row's distribution, near zero
  for a row far above the bound and approaching one for a row whose mean
  is below it (`dev/generics-trunc.R` sweeps three bounds on one
  design).

Measured on one 400-row lognormal fit against a two-chain brms fit of
the same model (`dev/generics-scale-brms.R`, seed 2026): the largest
relative disagreement in
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) over the 400
rows is 0.0106, which is 0.113 of brms's own posterior standard
deviation for that row and 0.058 of it at the median.
[`predict()`](https://rdrr.io/r/stats/predict.html) disagrees by roughly
440 at the median row and 765 at row 1, because it is a different scale
rather than a different answer. Those two are given loosely on purpose:
brms's [`predict()`](https://rdrr.io/r/stats/predict.html) summarizes
fresh predictive draws, so it moves between calls on one fit, where the
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) figures above
are deterministic and repeat to the last digit.

## See also

[`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md),
[`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md),
[`fitted.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/fitted.frmtmb_fit.md),
[`residuals.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/residuals.frmtmb_fit.md),
[`sigma.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/sigma.frmtmb_fit.md),
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md),
[`confint.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/confint.frmtmb_fit.md)

## Examples

``` r
set.seed(2026)
dd <- data.frame(x = rnorm(200))
dd$y <- exp(rnorm(200, 8 + 0.4 * dd$x, 0.4))
fit <- frm(bf(y ~ x) + lognormal(), data = dd)
#> Warning: Large maximum absolute gradient at the optimum (0.00143); the fit may not have converged. diagnose() names the offending parameter; see the 'Convergence problems' section of vignette('diagnostics') for the remedies

# predict() summarizes the predictive distribution; frm_linpred()
# is the linear predictor and fitted() the expected response
head(predict(fit, ndraws = 200))
#>      Estimate Est.Error      Q2.5    Q97.5
#> [1,] 3832.295 1535.2672 1715.7373 7631.886
#> [2,] 2162.523  931.8451 1051.5997 4591.016
#> [3,] 3580.354 1627.6969 1459.8562 7605.899
#> [4,] 3114.909 1282.5714 1335.0573 6272.218
#> [5,] 2405.242  946.9341  979.7062 4669.090
#> [6,] 1219.641  509.4481  512.4532 2496.949
head(frm_linpred(fit))
#>        1        2        3        4        5        6 
#> 8.214703 7.569772 8.061015 7.970745 7.736237 6.990888 
head(fitted(fit))
#>      Estimate Est.Error     Q2.5    Q97.5
#> [1,] 3984.808 126.97321 3735.945 4233.671
#> [2,] 2090.824  87.35777 1919.606 2262.042
#> [3,] 3417.128  98.16207 3224.734 3609.522
#> [4,] 3122.179  89.39280 2946.972 3297.386
#> [5,] 2469.517  84.57853 2303.746 2635.288
#> [6,] 1171.955  89.39567  996.743 1347.168

# the last two are related by the family's own mean
head(exp(frm_linpred(fit) + sigma(fit)^2 / 2) -
       fitted(fit)[, "Estimate"])
#> 1 2 3 4 5 6 
#> 0 0 0 0 0 0 

# a sigma nobody wrote a formula for is a distributional parameter
# on its own scale, not a coefficient on its link
summary(fit)$spec_pars
#>        Estimate  Est.Error  l-95% CI  u-95% CI
#> sigma 0.3886928 0.01943464 0.3524085 0.4287129
sigma(fit)
#> [1] 0.3886928
```
