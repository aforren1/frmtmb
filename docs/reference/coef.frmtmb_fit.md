# Per-group coefficients (fixed effects plus conditional modes)

Follows the lme4/glmmTMB/brms convention: for each random-effect
grouping factor, the fixed effects of its linear predictor broadcast
over the group levels, with the conditional modes added to the matching
columns. Use
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) for
the fixed effects alone.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
coef(object, summary = TRUE, robust = FALSE, probs = c(0.025, 0.975), ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- summary, robust, probs:

  brms's arguments, in brms's positions so that a positional brms call
  asks the same question. brms answers `summary = FALSE` with the
  posterior draws and `robust = TRUE` with their median and MAD, and a
  maximum-likelihood fit has no draws, so both are refused by name with
  the reason. The default of each is accepted and changes nothing.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

## Value

A named list of data frames, one per grouping factor, each with one row
per group level and one column per coefficient. When random effects
appear in more than one linear predictor, the list is nested one level
deeper, keyed as in
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md). A fit
without random effects returns the per-predictor coefficient vectors
instead.

## Details

The result is a list of data frames keyed by grouping factor. When
random effects appear in more than one dpar (or response), an outer
layer keyed like
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) is
added. Smooth terms are excluded. A fit without random effects returns
the coefficient vector of the location predictor (when there is one
linear predictor), or one vector per predictor.

The coefficients are named as brms names them, which is how
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) name their rows:
`Intercept`, not `(Intercept)`. Anything that pairs
[`coef()`](https://rdrr.io/r/stats/coef.html) with
[`vcov()`](https://rdrr.io/r/stats/vcov.html) by name,
[`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html)
among them, needs the two to agree. Use
[`fixef_by_dpar()`](https://aforren1.github.io/frmtmb/reference/fixef_by_dpar.md)
for the design-column spelling.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# one row per group: the fixed effects with the modes added in
head(coef(fit)$g)
#>     Intercept         x
#> 1  0.44665008 0.6706869
#> 2  1.43337233 0.6706869
#> 3  0.54268306 0.6706869
#> 4  0.87303571 0.6706869
#> 5 -0.02536969 0.6706869
#> 6  2.28544238 0.6706869
# which is fixef() plus ranef(), the lme4 identity, under the same
# names fixef() and vcov() use
all.equal(coef(fit)$g[["Intercept"]],
          fixef(fit)["Intercept", "Estimate"] + ranef(fit)$g[, 1],
          check.attributes = FALSE)
#> [1] TRUE

# without random effects there are no groups, so coef() is the
# coefficient vector of the location predictor
coef(frm(bf(y ~ x) + gaussian(), data = dd))
#> Intercept         x 
#> 1.2457711 0.6136201 
```
