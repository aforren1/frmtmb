# Likelihood-ratio tests between nested frmtmb fits

ML fits compare freely. REML fits compare only with each other, and only
when their fixed-effect designs span the same column space: a REML
likelihood is a likelihood for the error contrasts of that design, so
two of them are on a common scale exactly when the design is the same.
That covers the usual REML use - testing variance-component structures
with the fixed effects held fixed - and refuses the rest with the reason
(glmmTMB#776).

## Usage

``` r
# S3 method for class 'frmtmb_fit'
anova(object, ..., refit = FALSE)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Further `frmtmb_fit` objects, nested with `object`.

- refit:

  If `TRUE`, refit every REML fit in the comparison with ML and compare
  those, with a message naming what was refit. The refits reuse the
  assembled design and warm-start at the REML estimates. `FALSE` (the
  default) keeps the REML fits and refuses the comparisons a restricted
  likelihood cannot support.

## Value

An `anova` table.

## Details

`refit = TRUE` is the lme4 convenience for the refused case: every REML
fit in the comparison is refit with `REML = FALSE` and the ML fits are
compared instead. lme4 does this silently by default; here it is opt-in
and the message names the models that were refit.

When the smaller model removes a variance component, the null value sits
on the boundary of the parameter space and the usual chi-square
reference is wrong: the asymptotic null is a mixture (for one component,
half a point mass at zero and half a chi-square with one df), so the
reported p-value is conservative - up to a factor of two for a single
component. lme4 and glmmTMB report the same naive p-value; halve it for
the one-component case, or use
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
for a simulation-based reference.

## Nesting is assumed, not verified

A likelihood-ratio statistic has a chi-square null distribution only
when the smaller model is a restriction of the larger.
[`anova()`](https://rdrr.io/r/stats/anova.html) cannot verify that in
general: nesting through a nonlinear reparameterization, or through a
constraint that ties parameters across distributional parameters, is
invisible to anything the fitted objects carry. What it does check is
cheap and stated: if neither model's fixed-effect coefficient names are
a subset of the other's, and neither fixed-effect design sits inside the
other's column space, it warns. A comparison that passes that check is
not thereby verified to be nested. Two models that are genuinely not
nested are compared by AIC, not by this table.

## Examples

``` r
set.seed(1)
n <- 200
dd <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:20, 10)))
u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.5))
dd$y <- rnorm(n, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)

m0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
m1 <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd)
anova(m0, m1)
#> Likelihood-ratio tests
#> Each test assumes the smaller model is nested in the larger; see ?anova.frmtmb_fit
#> 
#>                     Df  logLik    AIC  Chisq Chi Df Pr(>Chisq)
#> y ~ x + (1 | g)      4 -334.52 677.04                         
#> y ~ x + z + (1 | g)  5 -334.05 678.10 0.9394      1     0.3324

# dropping a variance component puts the null on the boundary, so
# this p-value is conservative by up to a factor of two
m2 <- frm(bf(y ~ x) + gaussian(), data = dd)
anova(m2, m0)
#> Likelihood-ratio tests
#> Each test assumes the smaller model is nested in the larger; see ?anova.frmtmb_fit
#> 
#>                 Df  logLik    AIC  Chisq Chi Df Pr(>Chisq)    
#> y ~ x            3 -352.75 711.51                             
#> y ~ x + (1 | g)  4 -334.52 677.04 36.468      1  1.552e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# REML fits compare only when the fixed-effect designs agree, which
# is the case for a variance-component test
r0 <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd, REML = TRUE)
r1 <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd, REML = TRUE)
anova(r0, r1)
#> Likelihood-ratio tests
#> Each test assumes the smaller model is nested in the larger; see ?anova.frmtmb_fit
#> 
#>                 Df  logLik    AIC  Chisq Chi Df Pr(>Chisq)    
#> y ~ x + (1 | g)  2 -336.69 677.39                             
#> y ~ x + (x | g)  4 -325.44 658.87 22.516      2   1.29e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# differing designs are refused; refit = TRUE compares ML fits instead
rz <- frm(bf(y ~ x + z + (1 | g)) + gaussian(), data = dd, REML = TRUE)
try(anova(r0, rz))
#> Error : REML likelihoods are comparable only between fits whose fixed-effect designs span the same column space; these do not. Pass refit = TRUE (or refit with REML = FALSE) to compare fixed effects, or hold the fixed effects fixed to compare random-effect structures
anova(r0, rz, refit = TRUE)
#> anova(): refitting 2 REML models with ML: y ~ x + (1 | g); y ~ x + z + (1 | g)
#> Likelihood-ratio tests
#> Each test assumes the smaller model is nested in the larger; see ?anova.frmtmb_fit
#> 
#>                     Df  logLik    AIC  Chisq Chi Df Pr(>Chisq)
#> y ~ x + (1 | g)      4 -334.52 677.04                         
#> y ~ x + z + (1 | g)  5 -334.05 678.10 0.9394      1     0.3324
```
