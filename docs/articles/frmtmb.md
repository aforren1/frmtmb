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

## Smooths and Gaussian processes

mgcv smooths and `gp()` terms work in any linear predictor; their
smoothing and kernel parameters are estimated as variance components:

``` r

set.seed(2)
ds <- data.frame(x = runif(200, 0, 10))
ds$y <- sin(ds$x) + rnorm(200, 0, 0.4)
fs <- frm(bf(y ~ s(x)) + gaussian(), data = ds)
fg <- frm(bf(y ~ gp(x, k = 25)) + gaussian(), data = ds)
c(smooth = AIC(fs), gp = AIC(fg))
#>   smooth       gp 
#> 250.6414 249.8730
```

## Special terms

`mo()` fits monotonic effects of ordinal predictors, `mi()` imputes
missing continuous predictors inside the model (with `mi(sdx)` for known
measurement error), and `se()` gives meta-analysis:

``` r

set.seed(3)
dm <- data.frame(inc = sample(0:3, 300, TRUE), z = rnorm(300))
dm$y <- 1 + c(0, 1, 1.6, 2)[dm$inc + 1] + 0.3 * dm$z + rnorm(300)
fmo <- frm(bf(y ~ mo(inc) + z) + gaussian(), data = dm)
predict(fmo, newdata = data.frame(inc = 0:3, z = 0))
#>         1         2         3         4 
#> 0.9227861 1.8599425 2.5575740 2.9094146
```

## Mixtures

Finite mixtures keep per-component parameters and covariate-dependent
mixing weights; `groups = ~g` moves the class draw to the group level
(latent classes, combinable with random effects):

``` r

set.seed(4)
g <- rep(1:40, each = 5)
cl <- rbinom(40, 1, 0.4)
dmix <- data.frame(y = rnorm(200, c(-1, 2)[cl + 1][g], 0.8),
                   g = factor(g))
fmix <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
            data = dmix)
head(mixture_probs(fmix), 3)
#>   class1       class2
#> 1      1 3.936344e-08
#> 2      1 1.174035e-12
#> 3      1 1.276275e-16
```

## Inference beyond Wald

[`hypothesis()`](../reference/hypothesis.md) tests arbitrary parameter
expressions (Wald, profile, or parametric bootstrap), with natural-scale
random-effect names:

``` r

hypothesis(fit, "sd_Subject__Days^2 / (sd_Subject__Days^2 + sigma^2)")
#> Hypothesis tests (method = wald)
#>                                           hypothesis estimate      se      lwr
#>  sd_Subject__Days^2 / (sd_Subject__Days^2 + sigma^2)  0.04753 0.01989 0.008539
#>      upr     z       p
#>  0.08652 2.389 0.01688
```

[`frm_bootstrap()`](../reference/frm_bootstrap.md) runs a warm-started
parametric bootstrap, and `confint(fit, method = "profile")` profiles
single parameters.

## Diagnostics

`diagnose(fit)` reports convergence forensics. Simulation-based
residuals come through DHARMa, predictive checks through `pp_check()`,
and effect displays through
[`conditional_effects()`](../reference/conditional_effects.md):

``` r

plot(dharma_residuals(fit))
plot(conditional_effects(fit))
```

And the fitted objective can be handed to NUTS for full Bayes on the
same model - [`frm_sample()`](../reference/frm_sample.md) returns draws
with the full method surface
([`posterior_epred()`](../reference/posterior_epred.md),
[`hypothesis()`](../reference/hypothesis.md), `pp_check()`), and
[`check_laplace()`](../reference/check_laplace.md) audits the Laplace
approximation against them:

``` r

ds <- frm_sample(fit, chains = 4)
check_laplace(fit)
```

See [`vignette("brms-migration")`](../articles/brms-migration.md) for
the brms feature map and
[`vignette("diagnostics")`](../articles/diagnostics.md) for the
model-checking workflow.
