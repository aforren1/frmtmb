# Number of levels per random-effect grouping factor

Number of levels per random-effect grouping factor

## Usage

``` r
# S3 method for class 'frmtmb_draws'
ngrps(object, ...)

ngrps(object, ...)

# S3 method for class 'frmtmb_fit'
ngrps(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

## Value

A named integer vector (smooth terms are excluded).

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
#>  g  h 
#> 10  4 
# the count that decides whether a variance component is trustworthy,
# and the unit influence() deletes when given `groups`
ngrps(fit)[["h"]]
#> [1] 4
```
