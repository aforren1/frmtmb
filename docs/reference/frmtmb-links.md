# Link functions

The links frmtmb accepts, what each one maps, and where each one may be
used. Every family constructor takes a `link` argument for the mean and
a `link_<dpar>` argument for each of its other distributional
parameters, so `student(link_sigma = "softplus")` puts sigma on a
softplus. The roster follows brms 2.23.0: a model ported from brms is
fitted through the same inverse link.

## Value

A link is not a free-standing function. It is named to a family
constructor, which resolves it against the registry and carries it on
the fitted model, so the value a link contributes is the scale each
distributional parameter is estimated on. Read it back with `$link` and
`$link_<dpar>` on the family, and read a custom link's own fields with
the list that was supplied. This page documents the roster.

## Details

A link is defined here over plain arithmetic that RTMB overloads, not
through [`stats::make.link()`](https://rdrr.io/r/stats/make.link.html),
because [`make.link()`](https://rdrr.io/r/stats/make.link.html) clamps
at C level in ways the automatic differentiation tape cannot see.

## The roster

`g` is the link and `g^-1` its inverse. The inverse is the one that does
the work: it maps a linear predictor to the parameter, so its RANGE is
the parameter support the link fits. `eta` is the linear predictor and
`mu` the parameter.

|  |  |  |  |  |
|----|----|----|----|----|
| Link | `g(mu)` | `g^-1(eta)` | Range | Robust field |
| `identity` | `mu` | `eta` | the real line |  |
| `log` | `log(mu)` | `exp(eta)` | positive | `log_eta` |
| `softplus` | `log(expm1(mu))` | `log1p(exp(eta))` | positive | `log_eta` |
| `squareplus` | `(mu^2 - 1) / mu` | `(eta + sqrt(eta^2 + 4)) / 2` | positive | `log_eta` |
| `inverse` | `1 / mu` | `1 / eta` | positive |  |
| `1/mu^2` | `1 / mu^2` | `1 / sqrt(eta)` | positive |  |
| `sqrt` | `sqrt(mu)` | `eta^2` | non-negative | `log_eta` |
| `logit` | `log(mu / (1 - mu))` | `1 / (1 + exp(-eta))` | the unit interval | `logit_eta` |
| `probit` | `qnorm(mu)` | `pnorm(eta)` | the unit interval | `logit_eta` |
| `probit_approx` | `qnorm(mu)` | `plogis(0.07056 eta^3 + 1.5976 eta)` | the unit interval | `logit_eta` |
| `cloglog` | `log(-log(1 - mu))` | `1 - exp(-exp(eta))` | the unit interval | `logit_eta` |
| `cauchit` | `tan(pi (mu - 0.5))` | `0.5 + atan(eta) / pi` | the unit interval |  |
| `softit` | `log(expm1(mu / (1 - mu)))` | `s / (1 + s)` for `s = log1p(exp(eta))` | the unit interval | `logit_eta` |
| `logm1` | `log(mu - 1)` | `1 + exp(eta)` | above one |  |
| `log1p` | `log1p(mu)` | `expm1(eta)` | above minus one |  |
| `power12` | `log((mu - 1) / (2 - mu))` | `1 + 1 / (1 + exp(-eta))` | between one and two |  |
| `tan_half` | `tan(mu / 2)` | `2 atan(eta)` | the circle |  |

`probit_approx` is the one pair that is not an exact inverse, in frmtmb
and in brms alike. brms generates Stan code that uses `Phi_approx()`,
the logistic of a cubic, while its own R-side `inv_link()` answers
[`pnorm()`](https://rdrr.io/r/stats/Normal.html). The Stan form is the
one a ported model was fitted with, so it is the form used here.

`sqrt` is not injective on the whole line. brms allows it for a count
mean anyway, and so does frmtmb.

## Links for the mean

Each family constructor takes, for `mu`, the links brms 2.23.0 takes for
that family, and refuses the rest by name with the set it allows:
`bernoulli("sqrt")` is an error, as in brms, because a probability on
`eta^2` is not bounded by one. The sets are read out of brms by a script
rather than copied by hand. The first link listed is the default.

The link may be written unquoted, as in `student(identity)` or
`negbinomial(sqrt)`, as
[`stats::family()`](https://rdrr.io/r/stats/family.html) and brms allow.
Only the link for the mean can be; a `link_<dpar>` argument is a string.

A custom link object (see Custom links) is not held to the table: brms
has no such object, and a user who writes one has said what range the
mean has.
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
checks no set at all, because an extension family is free to mean
something else by its own `mu`.

`nbinom1`, `tweedie` and `huber` have no brms family. They take the set
of the nearest brms family: `negbinomial`, `Gamma` and `gaussian`, in
that order.

|  |  |
|----|----|
| Family | Links for `mu` |
| `gaussian`, `student`, `skew_normal` | `identity`, `log`, `inverse`, `softplus`, `squareplus`, and every unit-interval link |
| `lognormal`, `shifted_lognormal`, `hurdle_lognormal` | `identity`, `inverse` |
| `Gamma`, `weibull`, `exponential`, `hurdle_gamma` | `log`, `identity`, `inverse`, `softplus`, `squareplus` |
| `asym_laplace`, `exgaussian`, `zero_inflated_asym_laplace` | `identity`, `log`, `inverse`, `softplus`, `squareplus` |
| `poisson`, `negbinomial`, `geometric`, `compois`, and the zero-inflated and hurdle counts | `log`, `identity`, `sqrt`, `softplus`, `squareplus` |
| `binomial`, `bernoulli`, `Beta`, `zero_inflated_binomial`, `zero_inflated_beta`, `zero_one_inflated_beta` | `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit`, `softit`, `identity`, `log` |
| `beta_binomial` | the same list WITHOUT `log`: brms's own table for this one family stops at `identity` |
| `inverse.gaussian` | `1/mu^2`, `inverse`, `identity`, `log`, `softplus`, `squareplus` |
| `cox` | `log`, `identity`, `softplus`, `squareplus` |
| `von_mises` | `tan_half`, `identity` |
| `cumulative` | `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit`, `softit` |
| `sratio`, `cratio` | `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit` |
| `acat` | `logit` ONLY. brms takes the same six as `cumulative`; this is the one place frmtmb departs, and the reason is below |
| `categorical`, `multinomial` | `logit` |

An ordinal family's `link` names the cumulative distribution function
the thresholds are read through, not a link on a mean, so only a link
whose inverse maps onto the unit interval can serve.

[`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
is refused for a different reason, and it is the one place in this table
where frmtmb takes less than brms. Probit, cloglog, cauchit and softit
all map onto the unit interval and brms accepts every one of them for
`acat`, so the refusal is not about the link. It is about the density.
`brms:::inv_link_acat()` branches: on the logit a category probability
is `c(1, cumprod(exp(x)))` normalized, which is the log-linear form
frmtmb implements, and off the logit it is a product of distribution
functions times a reversed product of survivals, normalized. The second
form agrees with the first when the distribution function is logistic,
so it generalizes the same model rather than replacing it, but it is a
second expression that has to be written and taped. Substituting a
distribution function into the log-linear form does not reach it, which
is why the other three ordinal families could be routed through this
registry and
[`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
could not.
[`acat()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
says all of this when it refuses.

## Links for the other distributional parameters

`link_<dpar>` takes the values brms allows for that parameter, and the
set is fixed by what the parameter must stay inside. The first link
listed is the default.

|  |  |  |
|----|----|----|
| Parameter | Argument | Links |
| `sigma`, `shape`, `phi`, `kappa`, `beta`, `ndt`, and the `nu` of `compois` | `link_sigma` and so on | `log`, `identity`, `softplus`, `squareplus` |
| `nu` of `student` | `link_nu` | `logm1`, `identity` |
| `zi`, `hu`, `zoi`, `coi`, `quantile` | `link_zi`, `link_hu`, `link_zoi`, `link_coi`, `link_quantile` | `logit`, `identity` |
| `alpha` of `skew_normal` | `link_alpha` | `identity`, `log`, `softplus`, `squareplus` |

`identity` is in every one of those sets because brms puts it there. It
lets a parameter that must stay positive go negative, which the density
then reports as a `NaN` rather than as a bad link. It is offered so that
a brms model ports, not because it is a good choice.

The `power` of
[`tweedie()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
takes no argument. It is confined to the open interval from one to two,
`power12` is the only link in the roster that maps onto it, and brms has
no tweedie to port from.

A prior class is a distributional parameter name, so changing a
parameter's link changes where its prior sits.
`set_prior(class = "sigma")` is a density on sigma itself, and the
placement [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)
gives it is read from sigma's link, so changing `link_sigma` changes the
Jacobian that prior carries. See
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md).

## What a robust field buys

An inverse link SATURATES. `plogis(40)` is exactly one in double
precision, and so is `1 - exp(-exp(4))`. A density written over `1 - mu`
then reads `log(0)`, which is `-Inf` with an unusable gradient, where
the true log density is an ordinary `-40`. The linear predictor never
saturated, so a link that carries a robust field recovers the quantity
the density actually needs from the linear predictor directly:
`logit_eta` gives the log odds and `log_eta` the log mean.

The field is present only where the plain round trip measurably
saturates AND an exact form exists. `cauchit` has neither the problem
nor a cure: its tails are polynomial, so the plain round trip keeps
eight digits out to a linear predictor of 3.2e9. `inverse` and `1/mu^2`
cannot overflow or underflow at any representable positive linear
predictor.

Nothing has to be done to use this. A family whose link carries the
field is put on the robust form of its density automatically:
[`binomial()`](https://rdrr.io/r/stats/family.html) on
`dbinom_robust()`,
[`negbinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
on `dnbinom_robust()`.

## Custom links

Wherever a link is taken, a list is taken instead of a name. It must
carry `name`, `linkfun`, `linkinv` and `mu_eta`, the derivative of
`linkinv`, which `frm_linpred(se.fit = TRUE)` and every delta-method
interval read. `logit_eta` and `log_eta` are optional; supply one only
if it is exact.

## See also

[frmtmb-families](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
for the constructors that take these,
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
for what a link does to a prior.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(80))
d$y <- rbinom(80, 1, pnorm(0.4 + 0.8 * d$x))

# the mean on a probit rather than a logit
fixef_by_dpar(frm(bf(y ~ x), family = bernoulli(link = "probit"),
                  data = d))$mu
#> (Intercept)           x 
#>   0.2330309   1.0758635 

# a link on a parameter that is not the mean
d$z <- rnorm(80, 1 + d$x, exp(0.2 + 0.3 * d$x))
frm(bf(z ~ x, sigma ~ x), family = student(link_sigma = "log"), data = d)
#>  Family: student 
#>  Links: mu = identity; sigma = log; nu = logm1
#> 
#> Formula: z ~ x 
#>    Data: d (Number of observations: 80) 
#>  Method: ML   logLik: -129.955   AIC: 269.909   BIC: 281.82 
#> 
#> Regression Coefficients:
#>                 Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> Intercept           0.88      0.14     0.62     1.15    6.46  1.0e-10
#> sigma_Intercept     0.19      0.08     0.03     0.34    2.35    0.019
#> x                   0.83      0.16     0.51     1.15    5.08  3.9e-07
#> sigma_x             0.17      0.10    -0.02     0.36    1.74    0.082
#> 
#> Further Distributional Parameters:
#>     Estimate    Est.Error l-95% CI u-95% CI
#> nu 334785501 2.605674e+12        1      Inf
```
