# Extract fixed effects

Two shapes, and the difference is the NAMES. The default is a list keyed
by distributional parameter, whose entries are named by DESIGN COLUMN:
`fixef(fit)$sigma[["(Intercept)"]]`. `flatten = TRUE` is one vector in
the INTERNAL parameter name, which is what
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
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

- summary, robust, probs, pars:

  brms's arguments, in brms's positions so that a positional brms call
  asks the same question. brms answers `summary = FALSE` with the
  posterior draws and `robust = TRUE` with their median and MAD, and a
  maximum-likelihood fit has no draws, so both are refused by name with
  the reason. The default of each is accepted and changes nothing.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

- flatten:

  If `TRUE`, one named vector in the
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) /
  [`confint()`](https://rdrr.io/r/stats/confint.html) spelling instead
  of the per-dpar list. Coefficients come in linear-predictor order,
  which need not be [`vcov()`](https://rdrr.io/r/stats/vcov.html)'s row
  order; index by name. Distributional parameters held at a constant are
  included here and are absent from
  [`vcov()`](https://rdrr.io/r/stats/vcov.html), which covers the
  ESTIMATED coefficients only - so dividing by a
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) diagonal is `NA` for
  those entries, and `intersect(names(cf), rownames(vcov(fit)))` selects
  the ones a standard error exists for.

## Value

A named list of coefficient vectors, one per dpar, or with
`flatten = TRUE` a single named vector.

## Details

[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
and
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
use a THIRD vocabulary, brms's parameter names: `b_Intercept`,
`b_sigma_Intercept` for a `sigma` formula, and `sigma` on its natural
scale when no formula was written for it.
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
puts `class = "b"`'s `b_` in front of a bare name, so its default
spelling of those is `Intercept` and `sigma_Intercept`.

`unlist(fixef(fit))` is none of them. It is base R's composite of a list
KEY and an element name, `mu.x`, and it names no parameter of the model:
the separator differs and the location parameter is named where the
model does not name it. Use `flatten = TRUE` to line coefficients up
with a covariance matrix or a prior.

## See also

[`vcov.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/vcov.frmtmb_fit.md)
and
[`confint.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/confint.frmtmb_fit.md),
which name their rows the way `flatten = TRUE` names its entries.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)

# one entry per distributional parameter, each on its link scale
fit <- frm(bf(y ~ x + (1 | g), sigma ~ x) + gaussian(), data = dd)
fixef(fit)
#> $mu
#> (Intercept)           x 
#>    1.236875    0.594386 
#> 
#> $sigma
#> (Intercept)           x 
#>  -0.0959919   0.2082866 
#> 
exp(fixef(fit)$sigma[["(Intercept)"]])   # sigma is modeled on the log
#> [1] 0.9084714

# flatten to the vector vcov() and confint() name their rows by
fixef(fit, flatten = TRUE)
#>       (Intercept)                 x sigma_(Intercept)           sigma_x 
#>         1.2368754         0.5943860        -0.0959919         0.2082866 
all(names(fixef(fit, flatten = TRUE)) %in% rownames(confint(fit)))
#> [1] TRUE

# so a standard error goes with its coefficient by name
cf <- fixef(fit, flatten = TRUE)
cf / sqrt(diag(vcov(fit))[names(cf)])
#>       (Intercept)                 x sigma_(Intercept)           sigma_x 
#>          4.885125          5.413417         -1.280938          2.906674 
```
