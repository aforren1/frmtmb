# Single-term deletions

Drops each fixed-effect term of the primary (`mu`) formula in turn,
refits, and tabulates AIC (and likelihood-ratio tests with
`test = "Chisq"`), following
[`stats::drop1()`](https://rdrr.io/r/stats/add1.html) and lme4's
`drop1.merMod`. Random-effect, smooth, and `mo()`/`mi()` terms are not
part of the deletion scope.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
drop1(object, scope, test = c("none", "Chisq"), k = 2, ...)
```

## Arguments

- object:

  A `frmtmb_fit` from an ML fit (`REML = FALSE`) of a univariate model.

- scope:

  Terms to drop: a character vector or a right-hand-side formula.
  Defaults to all fixed-effect terms that marginality allows
  ([`stats::drop.scope()`](https://rdrr.io/r/stats/factor.scope.html)).

- test:

  `"Chisq"` adds likelihood-ratio tests.

- k:

  AIC penalty per parameter.

- ...:

  Unused.

## Value

An `anova` table with one row per dropped term.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), z = rnorm(100),
                 g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x * z + (1 | g)) + gaussian(), data = dd)

# marginality keeps the main effects out of scope while x:z is in it
drop1(fit)
#> Single term deletions
#> 
#> Model: y ~ x * z + (1 | g)
#> 
#>        Df    AIC
#> <none>    313.20
#> x:z     1 315.12
drop1(fit, test = "Chisq")
#> Single term deletions
#> 
#> Model: y ~ x * z + (1 | g)
#> 
#>        Df    AIC    LRT Pr(>Chi)  
#> <none>    313.20                  
#> x:z     1 315.12 3.9187  0.04775 *
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# name the terms to override the default scope
drop1(fit, scope = ~ x + z, test = "Chisq")
#> Single term deletions
#> 
#> Model: y ~ x * z + (1 | g)
#> 
#>        Df    AIC     LRT  Pr(>Chi)    
#> <none>    313.20                      
#> x       1 328.23 17.0327 3.674e-05 ***
#> z       1 312.12  0.9277    0.3355    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
