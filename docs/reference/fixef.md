# Extract fixed effects

Two shapes, and the difference is the NAMES. The default is brms's
summary matrix: one row per population-level coefficient under brms's
own name (`Intercept`, `sigma_Intercept`, `x`), and the columns
`Estimate`, `Est.Error`, `Q2.5` and `Q97.5`. A maximum-likelihood fit
has no draws, so `Est.Error` is the standard error and the `Q` columns
are the Wald interval at those probabilities, which is the interval
[`confint()`](https://rdrr.io/r/stats/confint.html) reports.
`flatten = TRUE` is one vector of estimates in the INTERNAL parameter
name, which is what [`confint()`](https://rdrr.io/r/stats/confint.html),
[`par_template()`](https://aforren1.github.io/frmtmb/reference/par_template.md),
`start` and `newparams` use: `dpar_column`, with the location
parameter's own coefficients left unprefixed (`x`, not `mu_x`), and the
response prefixed ahead of that in a multivariate fit.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
fixef(
  object,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  pars = NULL,
  ...,
  flatten = FALSE
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- summary, robust:

  brms's arguments, in brms's positions so that a positional brms call
  asks the same question. brms answers `summary = FALSE` with the
  posterior draws and `robust = TRUE` with their median and MAD, and a
  maximum-likelihood fit has no draws, so both are refused by name with
  the reason. The default of each is accepted and changes nothing.

- probs:

  Probabilities of the two quantile columns. brms takes the posterior
  quantiles there; here they are the ends of the Wald interval at those
  probabilities.

- pars:

  Row names to keep, in the order given, as in brms. A name the fit does
  not have is an error listing the ones it does.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

- flatten:

  If `TRUE`, one named vector of ESTIMATES in the
  [`confint()`](https://rdrr.io/r/stats/confint.html) spelling instead
  of the summary matrix. Coefficients come in linear-predictor order,
  which need not be the matrix's row order; index by name.
  Distributional parameters held at a constant are included here and are
  absent from [`vcov()`](https://rdrr.io/r/stats/vcov.html) and from the
  matrix, which cover the ESTIMATED population-level coefficients only.

## Value

A coefficients-by-four matrix, or with `flatten = TRUE` a single named
vector of estimates.

## Details

A distributional parameter nobody wrote a formula for is not a
population-level coefficient in brms and is not a row here: `sigma` of a
plain gaussian fit is reported by
[`summary()`](https://rdrr.io/r/base/summary.html) under
`Further Distributional Parameters`, on its own natural scale.

[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
and
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
use a THIRD vocabulary, brms's parameter names: `b_Intercept`,
`b_sigma_Intercept` for a `sigma` formula, and `sigma` on its natural
scale when no formula was written for it.
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
puts `class = "b"`'s `b_` in front of a bare name, so its default
spelling of those is `Intercept` and `sigma_Intercept`.

`unlist(fixef(fit, flatten = TRUE))` is none of them. Use
`flatten = TRUE` to line coefficients up with a prior or with
[`confint()`](https://rdrr.io/r/stats/confint.html), and the default
matrix to read an estimate with its standard error.

## See also

[`vcov.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/vcov.frmtmb_fit.md),
which names its rows as this matrix does, and
[`confint.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/confint.frmtmb_fit.md),
which names them as `flatten = TRUE` names its entries.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# brms's summary matrix, one row per population-level coefficient
fit <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
fixef(fit)
#>                   Estimate  Est.Error       Q2.5      Q97.5
#> Intercept        1.2368754 0.25319216  0.7406279 1.73312292
#> sigma_Intercept -0.0959919 0.07493873 -0.2428691 0.05088532
#> x                0.5943860 0.10979868  0.3791845 0.80958743
#> sigma_x          0.2082866 0.07165806  0.0678394 0.34873384
exp(fixef(fit)["sigma_Intercept", "Estimate"])  # sigma is on the log
#> [1] 0.9084714

# flatten to the vector confint() names its rows by
fixef(fit, flatten = TRUE)
#>       (Intercept)                 x sigma_(Intercept)           sigma_x 
#>         1.2368754         0.5943860        -0.0959919         0.2082866 
all(names(fixef(fit, flatten = TRUE)) %in% rownames(confint(fit)))
#> [1] TRUE

# the standard error column is vcov()'s diagonal
all.equal(fixef(fit)[, "Est.Error"], sqrt(diag(vcov(fit))))
#> [1] TRUE
```
