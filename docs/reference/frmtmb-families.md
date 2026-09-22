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
student(link = "identity", link_sigma = "log", link_nu = "logm1")

lognormal(link = "identity", link_sigma = "log")

negbinomial(link = "log", link_shape = "log")

nbinom1(link = "log", link_phi = "log")

Beta(link = "logit", link_phi = "log")

tweedie(link = "log", link_phi = "log")

compois(link = "log", link_nu = "log")

zero_inflated_poisson(link = "log", link_zi = "logit")

zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

hurdle_poisson(link = "log", link_hu = "logit")

multinomial(K)

cumulative(link = "logit")

beta_binomial(link = "logit", link_phi = "log")

skew_normal(link = "identity", link_sigma = "log", link_alpha = "identity")

exgaussian(link = "identity", link_sigma = "log", link_beta = "log")

bernoulli(link = "logit")

geometric(link = "log")

exponential(link = "log")

weibull(link = "log", link_shape = "log")

shifted_lognormal(link = "identity", link_sigma = "log", link_ndt = "log")

hurdle_gamma(link = "log", link_shape = "log", link_hu = "logit")

hurdle_lognormal(link = "identity", link_sigma = "log", link_hu = "logit")

zero_inflated_binomial(link = "logit", link_zi = "logit")

zero_inflated_beta(link = "logit", link_phi = "log", link_zi = "logit")

asym_laplace(link = "identity", link_sigma = "log", link_quantile = "logit")

zero_inflated_asym_laplace(
  link = "identity",
  link_sigma = "log",
  link_quantile = "logit",
  link_zi = "logit"
)

huber(link = "identity", k = 1.345, link_sigma = "log")

sratio(link = "logit")

cratio(link = "logit")

acat(link = "logit")

von_mises(link = "tan_half", link_kappa = "log")

categorical(link = "logit", levels = NULL, K = NULL)

cox(link = "log", df = 5, degree = 3, intercept = TRUE)
```

## Arguments

- link:

  Link for `mu`. See
  [frmtmb-links](https://aforren1.github.io/frmtmb/reference/frmtmb-links.md).

- link_sigma, link_shape, link_phi, link_kappa, link_ndt, link_beta:

  Link for a strictly positive parameter: one of `"log"` (the default),
  `"identity"`, `"softplus"` or `"squareplus"`.

- link_nu:

  Link for `nu`. `student()`'s degrees of freedom take `"logm1"` (the
  default) or `"identity"`, which keeps them above one; `compois()`'s
  dispersion is an ordinary positive parameter and takes the positive
  set.

- link_zi, link_hu, link_quantile:

  Link for a parameter on the unit interval: `"logit"` (the default) or
  `"identity"`.

- K:

  For `multinomial()`: number of response categories (columns of the
  count-matrix response); category 1 is the reference.

- link_alpha:

  Link for `skew_normal()`'s skewness, which is signed: `"identity"`
  (the default), `"log"`, `"softplus"` or `"squareplus"`.

- k:

  For `huber()`: Huber's tuning constant, the residual size in units of
  `sigma` where the density stops being gaussian and becomes Laplace. It
  is FIXED, not estimated - the default 1.345 is
  [`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html)'s, which gives
  95% efficiency against a gaussian.

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
ordered factor, or integer codes `1..K`. An unordered factor is refused,
as brms refuses it, because its level order is alphabetical unless
someone set it, and that order is the model.

A response with only two outcomes gets brms's message suggesting
`bernoulli()`: an ordinal or categorical response with two categories,
and a [`binomial()`](https://rdrr.io/r/stats/family.html),
`beta_binomial()` or `zero_inflated_binomial()` response whose trials
are all one.

Every constructor has brms's fields: `$link` is the name of the link for
the mean, and `$link_<dpar>` the link of each other parameter, as in
`beta_binomial()$link_phi`. See
[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md).

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
`frm_linpred(type = "response")` return the `n x K` matrix of category
probabilities, columns named by the response's own levels and rows
summing to one - the same convention the ordinal families follow.
`frm_linpred(type = "link")` and `frm_linpred(dpar =)` give the
per-category latent predictors, which is where `se.fit` lives.
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
`frm_linpred(type = "response")` report the mean direction. The
normalizing constant needs `log I0(kappa)`, which RTMB differentiates
exactly through its own `besselI` method, so nothing here is a series
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
`frm_linpred(type = "response")` are refused;
`frm_linpred(type = "link")` gives the log hazard ratio.
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

## Degrees of freedom that run off

`student()` estimates `nu` on the `logm1` link, which holds it above one
but puts no ceiling on it. A student-t reaches
[`gaussian()`](https://rdrr.io/r/stats/family.html) only in the limit
`nu -> Inf`, so on data with no heavy tails the likelihood keeps rising
as `nu` grows and the maximum is never attained. The fit then reports
whatever `nu` the optimizer last reached, with a standard error to
match: values of `1e9` and above are ordinary, and two runs of the same
model on the same data in a different row order can differ by orders of
magnitude in `nu` while agreeing on every coefficient to the last bit.

That is the model saying the data shows no heavy tails, not a failure to
converge.
[`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
names it under "Distributional parameter at the end of its link". Read
the fit as the gaussian one it has become. If you want a number you can
report, refit with [`gaussian()`](https://rdrr.io/r/stats/family.html),
or hold `nu` somewhere finite with a prior (see
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md));
brms does the same with its default `gamma(2, 0.1)`.

`nu` runs off the OTHER end too, and there the reading is opposite. On
tails heavier than any identified `nu` can hold, Cauchy data for
instance, `nu` goes down to one instead: `log(nu - 1)` around -20, a
standard error in the thousands, and a natural-scale value that prints
as `1`.
[`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
names that under the same heading. The remedy is not
[`gaussian()`](https://rdrr.io/r/stats/family.html), which is the worst
fit available for such data. There the data is the message; a prior is
still the way to hold `nu` finite.

The likelihood itself stays accurate the whole way. The log density is
formed so that `log Gamma((nu + 1) / 2) - log Gamma(nu / 2)` never
cancels, which holds it to 7e-15 of a 300-bit reference for every `nu`
up to `1e50`.

## Robust regression

`huber()` fits Huber's least-favorable distribution: gaussian within `k`
residual standard deviations of `mu` and Laplace outside, so a far-out
point pulls on the fit with a bounded influence instead of its squared
distance. It is a proper normalized density, not a penalty, so this is
ordinary maximum likelihood and
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and the likelihood-ratio
machinery all mean what they say.

`k` is a fixed constant of the family, `huber(k = 1.345)`, not a
distributional parameter. It states where the analyst draws the line
between a residual and an outlier, which is a modelling choice rather
than something the data identifies;
[`MASS::rlm()`](https://rdrr.io/pkg/MASS/man/rlm.html) treats it the
same way. Estimating it would let the likelihood buy fit by widening the
gaussian core, which is the opposite of the point. As `k` grows the
family collapses to [`gaussian()`](https://rdrr.io/r/stats/family.html).

Point estimates track `MASS::rlm(psi = psi.huber)` closely but not
exactly, and the difference is the scale: `rlm()` fixes the scale at a
MAD-type estimate and iterates the location, while `huber()` estimates
`sigma` by maximum likelihood jointly with `mu`. The coefficients agree
to about `1e-2` on well-behaved data; the `sigma` estimates need not.

The working-likelihood caveat that applies to `asym_laplace()` applies
here for the same reason. If the data are not actually
Huber-distributed - and the family is chosen precisely because the error
distribution is unknown - then the model is misspecified, the
information matrix is not the variance of the score, and Wald standard
errors and [`confint()`](https://rdrr.io/r/stats/confint.html) intervals
from it are not calibrated. The point estimates stay consistent for the
location. Use
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
for intervals you can defend. `cens()` and
[`trunc()`](https://rdrr.io/r/base/Round.html) are unavailable: the CDF
is a three-piece function of a parameter-dependent residual and has no
branch-free form for the tape.

`rho` has a kink at `|u| = k`, so the objective is only piecewise smooth
and the optimizer often stops with a maximum absolute gradient around
`1e-4` and the accompanying false-convergence warning. That is the kink,
not a bad fit: the same thing happens to `asym_laplace()`. The estimates
satisfy Huber's own estimating equations, `X' psi(u) = 0` with
`psi(u) = min(max(u, -k), k)`, to the same order.
[`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)
confirms the fit when in doubt.

## Links

Every constructor takes `link` for the mean and a `link_<dpar>` for each
of its other distributional parameters, following brms:
`student(link_sigma = "softplus")`,
`zero_inflated_poisson(link_zi = "identity")`.
[frmtmb-links](https://aforren1.github.io/frmtmb/reference/frmtmb-links.md)
lists the whole roster, what each link maps, which families take it for
the mean, and which set each parameter admits.

The four families 'stats' owns,
[`gaussian()`](https://rdrr.io/r/stats/family.html),
[`poisson()`](https://rdrr.io/r/stats/family.html),
[`binomial()`](https://rdrr.io/r/stats/family.html) and
[`Gamma()`](https://rdrr.io/r/stats/family.html), have no frmtmb
constructor to carry these. Reach their links through
[`brmsfamily()`](https://aforren1.github.io/frmtmb/reference/brmsfamily.md):
`brmsfamily("gaussian", link_sigma = "softplus")`.

The link for the mean may be unquoted, `student(identity)`, and is
refused by name when brms does not allow it for the family:
`bernoulli("sqrt")` is an error.

An ordinal family's `link` is not a link on a mean. It names the
distribution function the thresholds are read through, so
`cumulative()`, `sratio()` and `cratio()` take `logit`, `probit`,
`probit_approx`, `cloglog` and `cauchit` (and `cumulative()` also takes
`softit`), and refuse anything else. `acat()` takes `logit` alone,
because brms defines its other links by a different density rather than
by substituting a distribution function.

## Examples

``` r
set.seed(4)
n <- 120
dd <- data.frame(x = rnorm(n))

# heavier tails than gaussian(), with an estimated df
dd$y <- 1 + 0.8 * dd$x + rt(n, df = 4)
fixef(frm(bf(y ~ x) + student(), data = dd))
#>            Estimate  Est.Error      Q2.5    Q97.5
#> Intercept 0.9572315 0.09371241 0.7735586 1.140904
#> x         0.8342757 0.10463054 0.6292036 1.039348

# counts with more spread than poisson() allows
dd$cnt <- rnbinom(n, mu = exp(0.5 + 0.4 * dd$x), size = 2)
fit <- frm(bf(cnt ~ x) + negbinomial(), data = dd)
fixef_by_dpar(fit)$mu
#> (Intercept)           x 
#>   0.3690288   0.3371061 

# a zero-inflated count: the zi dpar gets its own predictor
dd$zi <- ifelse(runif(n) < 0.3, 0, dd$cnt)
frm(bf(zi ~ x, zi ~ 1) + zero_inflated_poisson(), data = dd)
#>  Family: zero_inflated_poisson 
#>  Links: mu = log; zi = logit
#> 
#> Formula: zi ~ x 
#>    Data: dd (Number of observations: 120) 
#>  Method: ML   logLik: -173.275   AIC: 352.55   BIC: 360.913 
#> 
#> Regression Coefficients:
#>              Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> Intercept        0.61      0.12     0.38     0.84    5.20    2e-07
#> zi_Intercept    -0.32      0.24    -0.80     0.16   -1.31    0.190
#> x                0.23      0.11     0.02     0.44    2.13    0.033

# an ordered response: level order is the category order
dd$grade <- cut(1 + 0.8 * dd$x + rlogis(n), 3,
                labels = c("low", "mid", "high"), ordered_result = TRUE)
frm(bf(grade ~ x) + cumulative(), data = dd)
#>  Family: cumulative 
#>  Links: cdf = logit
#> 
#> Formula: grade ~ x 
#>    Data: dd (Number of observations: 120) 
#>  Method: ML   logLik: -81.675   AIC: 169.35   BIC: 177.713 
#> 
#> Regression Coefficients:
#>              Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> Intercept[1]    -0.88      0.22    -1.32    -0.45   -4.01  6.0e-05
#> Intercept[2]     3.68      0.51     2.69     4.68    7.26  3.9e-13
#> x                1.11      0.26     0.61     1.61    4.33  1.5e-05

# a proportion in (0, 1)
dd$p <- plogis(0.2 + 0.6 * dd$x + rnorm(n, 0, 0.3))
frm(bf(p ~ x) + Beta(), data = dd)
#>  Family: beta 
#>  Links: mu = logit; phi = log
#> 
#> Formula: p ~ x 
#>    Data: dd (Number of observations: 120) 
#>  Method: ML   logLik: 154.805   AIC: -303.61   BIC: -295.247 
#> 
#> Regression Coefficients:
#>           Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> Intercept     0.21      0.03     0.16     0.27    8.07  7.1e-16
#> x             0.58      0.03     0.52     0.64   19.06  < 2e-16
#> 
#> Further Distributional Parameters:
#>     Estimate Est.Error l-95% CI u-95% CI
#> phi    50.46      6.45    39.27    64.83

# bounded influence: a few wild points barely move the slope
dd$rob <- 1 + 0.8 * dd$x + rnorm(n)
dd$rob[1:5] <- dd$rob[1:5] + 30
fixef_by_dpar(frm(bf(rob ~ x), family = huber(), data = dd))$mu
#> (Intercept)           x 
#>   1.1663906   0.8551163 
fixef_by_dpar(frm(bf(rob ~ x), family = gaussian(), data = dd))$mu
#> (Intercept)           x 
#>    2.231839    1.490395 

# an unordered factor: one predictor per non-reference category,
# named after the level it belongs to
dd$pick <- factor(sample(c("ale", "stout", "lager"), n, TRUE))
cat_fit <- frm(bf(pick ~ x), family = categorical(), data = dd)
fixef(cat_fit)                     # mulager and mustout; ale is the
#>                      Estimate Est.Error       Q2.5      Q97.5
#> mulager_Intercept  0.28202930 0.2503092 -0.2085678 0.77262640
#> mustout_Intercept  0.56454273 0.2363012  0.1014009 1.02768452
#> mulager_x         -0.08958449 0.2669129 -0.6127242 0.43355526
#> mustout_x         -0.43738957 0.2601125 -0.9472007 0.07242157
                                   # reference
head(fitted(cat_fit))              # n x K category probabilities
#> , , P(Y = ale)
#> 
#>       Estimate  Est.Error      Q2.5     Q97.5
#> [1,] 0.2564172 0.04058006 0.1768818 0.3359527
#> [2,] 0.2163825 0.04641332 0.1254141 0.3073510
#> [3,] 0.2928210 0.05531228 0.1844109 0.4012311
#> [4,] 0.2768567 0.04647044 0.1857763 0.3679371
#> [5,] 0.3327697 0.08839368 0.1595213 0.5060181
#> [6,] 0.2819019 0.04890836 0.1860433 0.3777606
#> 
#> , , P(Y = lager)
#> 
#>       Estimate  Est.Error      Q2.5     Q97.5
#> [1,] 0.3334248 0.04389351 0.2473951 0.4194545
#> [2,] 0.3011704 0.05134162 0.2005426 0.4017981
#> [3,] 0.3584389 0.05855308 0.2436769 0.4732008
#> [4,] 0.3479778 0.04992207 0.2501324 0.4458233
#> [5,] 0.3810589 0.08940322 0.2058318 0.5562860
#> [6,] 0.3513701 0.05233584 0.2487938 0.4539465
#> 
#> , , P(Y = stout)
#> 
#>       Estimate  Est.Error      Q2.5     Q97.5
#> [1,] 0.4101580 0.04632187 0.3193688 0.5009472
#> [2,] 0.4824471 0.05570475 0.3732678 0.5916264
#> [3,] 0.3487401 0.05923470 0.2326422 0.4648380
#> [4,] 0.3751655 0.05205583 0.2731379 0.4771930
#> [5,] 0.2861714 0.07905345 0.1312295 0.4411134
#> [6,] 0.3667279 0.05416057 0.2605751 0.4728807
#> 

# one category may take its own predictor
dd$w <- rnorm(n)
frm(bf(pick ~ x, mustout ~ w), family = categorical(), data = dd)
#>  Family: categorical 
#>  Links: mulager = identity; mustout = identity
#> 
#> Formula: pick ~ x 
#>    Data: dd (Number of observations: 120) 
#>  Method: ML   logLik: -128.572   AIC: 265.145   BIC: 276.295 
#> 
#> Regression Coefficients:
#>                   Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> mulager_Intercept     0.24      0.24    -0.24     0.72    0.99    0.321
#> mustout_Intercept     0.51      0.23     0.06     0.97    2.22    0.027
#> mulager_x             0.19      0.21    -0.24     0.61    0.86    0.388
#> mustout_w             0.10      0.19    -0.28     0.48    0.52    0.600

# an angle: mu is the mean direction, kappa the concentration
dd$angle <- atan2(sin(0.5 + dd$x), cos(0.5 + dd$x))
vm_fit <- frm(bf(angle ~ x), family = von_mises(), data = dd)
head(fitted(vm_fit))               # the mean direction, in radians
#>         Estimate  Est.Error        Q2.5       Q97.5
#> [1,]  0.82530547 0.01603449  0.79387845  0.85673250
#> [2,] -0.05745388 0.01824795 -0.09321921 -0.02168855
#> [3,]  1.41150621 0.01768376  1.37684668  1.44616573
#> [4,]  1.18172554 0.01719206  1.14802971  1.21542136
#> [5,]  1.83731815 0.01706177  1.80387769  1.87075861
#> [6,]  1.25871540 0.01740168  1.22460872  1.29282207

# proportional hazards with a spline baseline; (1 | g) is a frailty
dd$time <- rexp(n, exp(-0.5 + 0.7 * dd$x))
dd$out <- rbinom(n, 1, 0.3)        # 1 = right censored
cox_fit <- frm(bf(time | cens(out) ~ x), family = cox(), data = dd)
fixef_by_dpar(cox_fit)$mu                  # log hazard ratios
#> (Intercept)           x 
#>   2.1265524   0.7045597 
cox_baseline(cox_fit)              # the baseline hazard weights
#>           s1           s2           s3           s4           s5 
#> 3.668968e-03 1.910337e-01 2.186349e-01 4.645067e-09 5.866624e-01 
```
