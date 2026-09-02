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

von_mises(link = "tan_half")

categorical(link = "logit", levels = NULL, K = NULL)

cox(link = "log", df = 5, degree = 3, intercept = TRUE)
```

## Arguments

- link:

  Link for `mu`.

- K:

  For `multinomial()`: number of response categories (columns of the
  count-matrix response); category 1 is the reference.

- levels:

  For `categorical()`: the response's category labels, in the order that
  fixes the reference category (the first) and the dpar names. Only
  needed when the family is built away from the data;
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) reads
  them off the response.

- df:

  For `cox()`: the number of M-spline basis functions in the baseline
  hazard (brms's `bhaz(df = 5)` default).

- degree:

  For `cox()`: the spline degree of that basis (cubic by default).

- intercept:

  For `cox()`: keep the basis function that is non-zero at the lower
  boundary knot.

## Value

A `frmtmb_family` object.

## Details

An ordinal family (`cumulative()`, `sratio()`, `cratio()`, `acat()`)
takes the response's level order as the category order. Supply an
ordered factor, or integer codes `1..K`: an unordered factor is
accepted, as brms accepts it, but warns and names the order it is about
to use, which is alphabetical unless the levels were set.

## Categorical (nominal) responses

`categorical()` fits a multinomial logit to an unordered factor. The
FIRST level is the reference category, its linear predictor is held at
zero, and each remaining level gets its own predictor named `mu<Level>`,
as in brms. The main model formula applies to every one of them, and any
single category is overridden by naming it: `bf(y ~ x, mub ~ z)` gives
category `"b"` its own predictor and leaves the rest on `~ x`.

The categories come from the data, so
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) reads them
off the response before it parses the formula. Constructing the family
away from a data set (to inspect it, or to reach the parser through
another entry point) needs them stated:
`categorical(levels = c("a", "b", "c"))`, or `categorical(K = 3)` for a
response already coded `1..K`, which names the dpars `mu2 ... muK`. A
character or logical response is coerced to a factor with a message
naming the level order, because that order is the model.

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
`predict(type = "response")` return the `n x K` matrix of category
probabilities, columns named by the response's own levels and rows
summing to one - the same convention the ordinal families follow.
`predict(type = "link")` and `predict(dpar =)` give the per-category
latent predictors, which is where `se.fit` lives.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws factor
levels. The same likelihood is available on a count-matrix response as
`multinomial(K)`, and a one-hot matrix gives an identical
log-likelihood.

## Circular responses

`von_mises()` models an angle in radians on `(-pi, pi]`. Its `mu` is the
mean direction and takes the `tan_half` link, which maps the whole line
onto that interval; `kappa` is the concentration and takes a log link,
with `kappa = 0` the uniform distribution on the circle. Both are brms's
choices. [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
`predict(type = "response")` report the mean direction. The normalizing
constant needs `log I0(kappa)`, which RTMB differentiates exactly
through its own `besselI` method, so nothing here is a series
approximation. Residuals are differences of angles and are NOT wrapped,
so read [`residuals()`](https://rdrr.io/r/stats/residuals.html) on a von
Mises fit with that in mind.

## Cox proportional hazards

`cox()` is the flexible-parametric proportional hazards model brms fits:
the baseline hazard is an M-spline in time, `h0(t) = sum_j s_j M_j(t)`,
and the cumulative baseline hazard is the I-spline integral of the same
basis over the same weights, which makes it monotone by construction.
The weights `s` form a simplex - that is what identifies the baseline
against the intercept - and are estimated as `sbhaz_raw`, their
`Kbhaz - 1` free softmax coordinates;
[`cox_baseline()`](https://aforren1.github.io/frmtmb/reference/cox_baseline.md)
returns the simplex itself. The default basis is brms's: `df = 5` cubic
M-splines with an intercept, internal knots on response quantiles and
boundary knots just outside the observed range.

The hazard is `h0(t) exp(eta)`, so a coefficient is a log hazard ratio,
exactly as in
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html).
Right, left, and interval censoring come through the ordinary `cens()`
addition term: an event contributes the density and a censored
observation the survivor function, which is what this family's
log-density and log-CDF are. Random effects are the point:
`time | cens(c) ~ x + (1 | g)` is a frailty model, and the Laplace
approximation integrates the frailties out.

The baseline is semiparametric only in spirit - it has `df` parameters,
not one per event time - so coefficients agree with `coxph()` closely
rather than exactly. A survival response has no mean, so
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
`predict(type = "response")` are refused; `predict(type = "link")` gives
the log hazard ratio.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) is not available.

Maximum likelihood often puts one or more baseline weights ON the
simplex boundary, at exactly zero. Their softmax coordinates then run
off to minus infinity along a flat ridge, the Hessian is singular in
those directions, and the optimizer reports singular convergence even
though the gradient is zero and the regression coefficients are at their
optimum -
[`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
names `sbhaz_raw` as the culprit. This is what an unpenalized flexible
baseline does; brms does not meet it because its Dirichlet prior keeps
the weights interior. Lower `df` until the baseline is one the data
supports, and read
[`cox_baseline()`](https://aforren1.github.io/frmtmb/reference/cox_baseline.md)
to see which weights collapsed.

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

# an unordered factor: one predictor per non-reference category,
# named after the level it belongs to
dd$pick <- factor(sample(c("ale", "stout", "lager"), n, TRUE))
cat_fit <- frm(bf(pick ~ x), family = categorical(), data = dd)
fixef(cat_fit)                     # mulager and mustout; ale is the
#> $mulager
#> (Intercept)           x 
#>   0.3644546  -0.1122166 
#> 
#> $mustout
#> (Intercept)           x 
#> -0.06738302 -0.30584099 
#> 
                                   # reference
head(fitted(cat_fit))              # n x K category probabilities
#>         ale     lager     stout
#> 1 0.3048775 0.4283931 0.2667294
#> 2 0.2752055 0.4210911 0.3037034
#> 3 0.3317258 0.4321454 0.2361288
#> 4 0.3199327 0.4308183 0.2492490
#> 5 0.3616400 0.4333564 0.2050036
#> 6 0.3236542 0.4312905 0.2450553

# one category may take its own predictor
dd$w <- rnorm(n)
frm(bf(pick ~ x, mustout ~ w), family = categorical(), data = dd)
#> frmtmb fit: pick ~ x 
#> Family: categorical   Method: ML 
#> logLik: -129.366  AIC: 266.732  nobs: 120 
#> 
#> Fixed effects:
#>  mulager:
#> (Intercept)           x 
#>     0.34532     0.03655 
#>  mustout:
#> (Intercept)           w 
#>     -0.1019     -0.1376 

# an angle: mu is the mean direction, kappa the concentration
dd$angle <- atan2(sin(0.5 + dd$x), cos(0.5 + dd$x))
vm_fit <- frm(bf(angle ~ x), family = von_mises(), data = dd)
head(fitted(vm_fit))               # the mean direction, in radians
#>           1           2           3           4           5           6 
#>  0.82530547 -0.05745388  1.41150621  1.18172554  1.83731815  1.25871540 

# proportional hazards with a spline baseline; (1 | g) is a frailty
dd$time <- rexp(n, exp(-0.5 + 0.7 * dd$x))
dd$out <- rbinom(n, 1, 0.3)        # 1 = right censored
cox_fit <- frm(bf(time | cens(out) ~ x), family = cox(), data = dd)
fixef(cox_fit)$mu                  # log hazard ratios
#> (Intercept)           x 
#>   2.4319804   0.9064254 
cox_baseline(cox_fit)              # the baseline hazard weights
#>           s1           s2           s3           s4           s5 
#> 1.292303e-02 1.489409e-01 4.591886e-01 4.732665e-09 3.789474e-01 
```
