# Extract random-effect standard deviations and correlations

brms's `VarCorr()`, with frequentist content. The return value is a list
with one entry per GROUPING FACTOR, named by the factor as brms and lme4
name it (`"patient"`, not the term label `"1 | patient"`), followed by
`residual__` when `sigma` is a parameter nobody wrote a formula for (one
row per response on a multivariate model, with the residual correlations
under `rescor = TRUE`). Each entry has:

## Usage

``` r
# S3 method for class 'frmtmb_fit'
VarCorr(
  x,
  sigma = 1,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)
```

## Arguments

- x:

  A `frmtmb_fit`.

- sigma:

  Ignored, as in brms. It is carried by nlme's generic, which frmtmb
  shares.

- summary:

  Must be `TRUE`; see above.

- robust:

  Must be `FALSE`; see above.

- probs:

  The quantiles to report.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

## Value

A named list, as above.

## Details

- `sd`: a matrix with one row per coefficient and the columns
  `Estimate`, `Est.Error` and one quantile column per `probs` (`Q2.5`,
  `Q97.5`);

- `cor` and `cov`, when the group has correlations: arrays of
  `coefficient x statistic x coefficient`, the same statistics.

The names are brms's: the entry is the group as brms spells it (`g:h`
for an interaction), and the coefficients are `Intercept`, `x`, and for
a distributional or nonlinear parameter `sigma_Intercept`, as in
`sd_<group>__sigma_Intercept`.

## What the columns mean on a maximum-likelihood fit

brms summarizes posterior draws. A fit has one estimate and its sampling
distribution, so the same columns carry:

- `Estimate`: the maximum-likelihood (or REML) estimate.

- `Est.Error`: its delta-method standard error, from the joint
  covariance of the covariance parameters
  ([`vcov()`](https://rdrr.io/r/stats/vcov.html) with `full = TRUE`).

- `Q<p>`: the Wald quantile `Estimate + qnorm(p) * Est.Error`, on the
  natural scale. That is the interval
  [`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
  reports for the same quantity. It can cross zero for a standard
  deviation near its boundary;
  [`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
  with `method = "profile"` gives an interval that respects the
  boundary.

A diagonal entry of `cor` is 1 with error 0.

`summary = FALSE` and `robust = TRUE` are refused by name: brms returns
the draws, or their median and MAD, and a fit has no draws.
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
gives an object whose `VarCorr()` does both.

Several random-effect blocks on one grouping factor merge into its
entry, as brms merges every term on a group, and a pair across two
blocks has correlation 0. Two terms that give one group the same
coefficient, as in `(1 | gr(id, cov = A)) + (1 | id)`, are refused when
the model is built, as brms refuses them: give the second term a copy of
the grouping column under another name. Smooth, Gaussian-process, CAR
and SPDE blocks are not in brms's `VarCorr()` and are not here:
[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
reports them. On a `gr(dist = "student")` block the `sd` row is the
SCALE, as brms's `sd_` is.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:20, 10)))
u <- cbind(rnorm(20, 0, 0.8), rnorm(20, 0, 0.4))
dd$y <- rnorm(200, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x, 1)
fit <- frm(bf(y ~ x + (x | g)) + gaussian(), data = dd)

vc <- VarCorr(fit)
names(vc)                       # "g" and "residual__", as in brms
#> [1] "g"          "residual__"
vc$g$sd                         # standard deviations
#>            Estimate Est.Error      Q2.5     Q97.5
#> Intercept 0.8861648 0.1581149 0.5762653 1.1960643
#> x         0.4512135 0.1239615 0.2082534 0.6941735
vc$g$cor["Intercept", "Estimate", "x"]
#> [1] 0.05273605
vc$residual__$sd
#>   Estimate  Est.Error      Q2.5    Q97.5
#>  0.9641289 0.05388581 0.8585147 1.069743
```
