# Number of levels per random-effect grouping factor

Number of levels per random-effect grouping factor

## Usage

``` r
ngrps(object, ...)

# S3 method for class 'frmtmb_fit'
ngrps(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

## Value

brms's named list, one integer per grouping factor, or `NULL` for a fit
with no grouping factor. Smooth, Gaussian-process, CAR and SPDE blocks
are random-effect blocks here and are not grouping factors in brms, so
they are excluded; their parameters are in
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
and
[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100),
                 g = factor(rep(1:10, 10)),
                 h = factor(rep(1:4, each = 25)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g) + (1 | h)) + gaussian(), data = dd)

# one count per distinct grouping factor
ngrps(fit)
#> $g
#> [1] 10
#> 
#> $h
#> [1] 4
#> 
# the count that decides whether a variance component is trustworthy,
# and the unit influence() deletes when given `groups`
ngrps(fit)[["h"]]
#> [1] 4
```
