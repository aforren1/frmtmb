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
| [`predict()`](https://rdrr.io/r/stats/predict.html), default | the `mu` linear predictor | link |
| `predict(type = "link")` | the same | link |
| `predict(type = "response")` | the conditional mean | response |
| `predict(dpar = "sigma")` | that predictor | link; `type = "response"` gives response |
| [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) | the conditional mean | response |
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

## The two places this differs from brms

[`predict()`](https://rdrr.io/r/stats/predict.html) is the first. brms
returns posterior predictive draws summarized on the RESPONSE scale;
frmtmb returns the LINEAR PREDICTOR, because it has one parameter vector
rather than a posterior and the linear predictor is the quantity an ML
fit actually estimates. On a lognormal fit that is the difference
between about 8 and about 6,700. Ask for `type = "response"`, or call
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), to get the
quantity brms's [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
reports. Measured with `dev/generics-scale.R` and
`dev/generics-scale-brms.R`.

`sigma` is the second.
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) report
every coefficient on its own link, so a `sigma` with the default log
link is printed as `log(sigma)` and can be negative. brms samples
`sigma` itself and prints it on the response scale.
[`sigma()`](https://rdrr.io/r/stats/sigma.html) here back-transforms and
gives brms's number.

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
  gives `NA` with it. Use `predict(dpar = "sigma", type = "response")`,
  which reproduces
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) exactly
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

# the default is the linear predictor, not the outcome
head(predict(fit))
#>        1        2        3        4        5        6 
#> 8.214703 7.569772 8.061015 7.970745 7.736237 6.990888 
head(fitted(fit))
#>        1        2        3        4        5        6 
#> 3984.808 2090.824 3417.128 3122.179 2469.517 1171.955 

# and the two are related by the family's own mean
head(exp(predict(fit) + sigma(fit)^2 / 2) - fitted(fit))
#> 1 2 3 4 5 6 
#> 0 0 0 0 0 0 

# sigma is printed on its log link and back-transformed by sigma()
summary(fit)$coefficients$sigma[1, 1]
#> [1] -0.944966
sigma(fit)
#> [1] 0.3886928
```
