# Getting started with frmtmb

frmtmb fits regression models written in the brms formula grammar by
maximum likelihood, with the Laplace approximation for random effects.
The backend is RTMB: the model formula compiles to an R objective
function that is differentiated on the TMB AD tape. There is no MCMC and
no compilation at run time.

## A mixed model

``` r

library(frmtmb)
data(sleepstudy, package = "lme4")

fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
           data = sleepstudy)
summary(fit)
#> Family: gaussian 
#> Formula: Reaction ~ Days + (Days | Subject) 
#> Method: ML   nobs: 180 
#> logLik: -875.97  AIC: 1763.94  BIC: 1783.1 
#> 
#> Random effects:
#>   Days | Subject 
#>         Name Std.Dev. (Intercept)
#>  (Intercept)  23.7800            
#>         Days   5.7168      0.0813
#> 
#> Coefficients (mu):
#>             Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 251.4052     6.6323 37.9064 < 2.2e-16
#> Days         10.4673     1.5022  6.9678  3.22e-12
#> 
#> Coefficients (sigma):
#>             Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) 3.242273   0.058926  55.023 < 2.2e-16
```

The usual methods work: [`fixef()`](../reference/fixef.md),
[`ranef()`](../reference/ranef.md),
[`VarCorr()`](../reference/VarCorr.md),
[`predict()`](https://rdrr.io/r/stats/predict.html) (with `newdata`,
`se.fit`, and `re.form`),
[`confint()`](https://rdrr.io/r/stats/confint.html) (Wald, profile, or
likelihood-root), [`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`anova()`](https://rdrr.io/r/stats/anova.html) for likelihood-ratio
tests, and `emmeans` for marginal means.

## Distributional regression

Every distributional parameter can get its own formula with the full
predictor grammar, including random effects and smooths:

``` r

fit2 <- frm(bf(Reaction ~ Days + (Days | Subject),
               sigma ~ Days) + gaussian(),
            data = sleepstudy)
fixef(fit2)$sigma
#> (Intercept)        Days 
#>  2.81184862  0.08463401
```

## Nonlinear formulas

With `nl = TRUE` the model formula is an arbitrary R expression over
named parameters, each of which has its own linear predictor:

``` r

set.seed(1)
d <- data.frame(x = runif(200, 0, 5))
d$y <- 2.5 * exp(-0.7 * d$x) + rnorm(200, 0, 0.15)
nlfit <- frm(bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE) +
               gaussian(),
             data = d, start = list(beta = c(1, 0.3)))
unlist(fixef(nlfit)[c("a", "b")])
#> a.(Intercept) b.(Intercept) 
#>     2.5067769     0.7022003
```

## Custom families

A family is a plain R function returning an AD-safe log-density; no Stan
or C++ involved:

``` r

my_poisson <- custom_family(
  "my_poisson", dpars = "mu", links = list(mu = "log"),
  lpdf = function(y, dpars, aterms) {
    y * log(dpars$mu) - dpars$mu - lgamma(y + 1)
  }
)
check_custom_family(my_poisson, y = rpois(50, 3),
                    dpars = list(mu = rep(2.5, 50)))
```

## Diagnostics

`diagnose(fit)` reports convergence forensics. Simulation-based
residuals come through DHARMa:

``` r

plot(dharma_residuals(fit))
```

And the fitted objective can be handed to NUTS for full Bayes on the
same model, where bayesplot works on the draws:

``` r

sf <- as_tmbstan(fit, chains = 4)
bayesplot::mcmc_trace(rstan::extract(sf, permuted = FALSE))
```
