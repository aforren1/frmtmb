# Additional response families

Family constructors without a
[stats::family](https://rdrr.io/r/stats/family.html) equivalent,
following brms naming.
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html), and
[`Gamma()`](https://rdrr.io/r/stats/family.html) from 'stats' are
accepted directly by
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) and `+`.

## Usage

``` r
student(link = "identity")

lognormal(link = "identity")

negbinomial(link = "log")

nbinom1(link = "log")

Beta(link = "logit")

tweedie(link = "log")

compois(link = "log")

zero_inflated_poisson(link = "log")

zero_inflated_negbinomial(link = "log")

hurdle_poisson(link = "log")

multinomial(K)

cumulative(link = "logit")

beta_binomial(link = "logit")

skew_normal(link = "identity")

exgaussian(link = "identity")

bernoulli(link = "logit")

geometric(link = "log")

exponential(link = "log")

weibull(link = "log")

shifted_lognormal(link = "identity")

hurdle_gamma(link = "log")

hurdle_lognormal(link = "identity")

zero_inflated_binomial(link = "logit")

zero_inflated_beta(link = "logit")

asym_laplace(link = "identity")

zero_inflated_asym_laplace(link = "identity")

sratio(link = "logit")

cratio(link = "logit")

acat(link = "logit")
```

## Arguments

- link:

  Link for `mu`.

- K:

  For `multinomial()`: number of response categories (columns of the
  count-matrix response); category 1 is the reference.

## Value

A `frmtmb_family` object.

## Details

An ordinal family (`cumulative()`, `sratio()`, `cratio()`, `acat()`)
takes the response's level order as the category order. Supply an
ordered factor, or integer codes `1..K`: an unordered factor is
accepted, as brms accepts it, but warns and names the order it is about
to use, which is alphabetical unless the levels were set.

## Quantile regression inference

`asym_laplace()` and `zero_inflated_asym_laplace()` fit quantile
regression through a WORKING likelihood: at a fixed `quantile` the point
estimates are consistent quantile estimates (they match
[`quantreg::rq()`](https://rdrr.io/pkg/quantreg/man/rq.html)), but the
asymmetric Laplace is not the data's true density, so Wald standard
errors and [`confint()`](https://rdrr.io/r/stats/confint.html) intervals
computed from it are not calibrated. This is a property of the
asymmetric-Laplace approach, shared by brms. Use
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
for intervals you can defend. The check function's kink can also produce
a benign false-convergence warning near the optimum;
[`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)
confirms the fit when in doubt.

## Examples

``` r
set.seed(4)
n <- 120
dd <- data.frame(x = rnorm(n))

# heavier tails than gaussian(), with an estimated df
dd$y <- 1 + 0.8 * dd$x + rt(n, df = 4)
fixef(frm(bf(y ~ x) + student(), data = dd))
#> $mu
#> (Intercept)           x 
#>   0.9572315   0.8342757 
#> 
#> $sigma
#> (Intercept) 
#>  -0.2075767 
#> 
#> $nu
#> (Intercept) 
#>   0.4055075 
#> 

# counts with more spread than poisson() allows
dd$cnt <- rnbinom(n, mu = exp(0.5 + 0.4 * dd$x), size = 2)
fit <- frm(bf(cnt ~ x) + negbinomial(), data = dd)
fixef(fit)$mu
#> (Intercept)           x 
#>   0.3690288   0.3371061 

# a zero-inflated count: the zi dpar gets its own predictor
dd$zi <- ifelse(runif(n) < 0.3, 0, dd$cnt)
frm(bf(zi ~ x, zi ~ 1) + zero_inflated_poisson(), data = dd)
#> frmtmb fit: zi ~ x 
#> Family: zero_inflated_poisson   Method: ML 
#> logLik: -173.275  AIC: 352.55  nobs: 120 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      0.6070      0.2271 
#>  zi:
#> (Intercept) 
#>     -0.3201 

# an ordered response: level order is the category order
dd$grade <- cut(1 + 0.8 * dd$x + rlogis(n), 3,
                labels = c("low", "mid", "high"), ordered_result = TRUE)
frm(bf(grade ~ x) + cumulative(), data = dd)
#> frmtmb fit: grade ~ x 
#> Family: cumulative   Method: ML 
#> logLik: -81.675  AIC: 169.35  nobs: 120 
#> 
#> Fixed effects:
#>  mu:
#>     x 
#> 1.108 

# a proportion in (0, 1)
dd$p <- plogis(0.2 + 0.6 * dd$x + rnorm(n, 0, 0.3))
frm(bf(p ~ x) + Beta(), data = dd)
#> frmtmb fit: p ~ x 
#> Family: beta   Method: ML 
#> logLik: 154.805  AIC: -303.61  nobs: 120 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      0.2133      0.5771 
#>  phi:
#> (Intercept) 
#>       3.921 
```
