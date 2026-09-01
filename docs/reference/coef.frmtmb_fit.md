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
coef(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

## Value

A named list of data frames, one per grouping factor, each with one row
per group level and one column per coefficient. When random effects
appear in more than one linear predictor, the list is nested one level
deeper, keyed as in
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md). A fit
without random effects returns the
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) value
instead.

## Details

The result is a list of data frames keyed by grouping factor. When
random effects appear in more than one dpar (or response), an outer
layer keyed like
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) is
added. Smooth terms are excluded. A fit without random effects returns
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) (the
single coefficient vector when there is one linear predictor).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# one row per group: the fixed effects with the modes added in
head(coef(fit)$g)
#>   (Intercept)         x
#> 1  0.44665008 0.6706869
#> 2  1.43337233 0.6706869
#> 3  0.54268306 0.6706869
#> 4  0.87303571 0.6706869
#> 5 -0.02536969 0.6706869
#> 6  2.28544238 0.6706869
# which is fixef() plus ranef(), the lme4 identity
all.equal(coef(fit)$g[["(Intercept)"]],
          fixef(fit)$mu[["(Intercept)"]] + ranef(fit)$g[, 1],
          check.attributes = FALSE)
#> [1] "Numeric: lengths (10, 0) differ"

# without random effects there are no groups, so coef() is fixef()
coef(frm(bf(y ~ x) + gaussian(), data = dd))
#> (Intercept)           x 
#>   1.2457711   0.6136201 
```
