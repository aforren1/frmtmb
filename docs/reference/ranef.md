# Extract random-effect modes

Extract random-effect modes

## Usage

``` r
ranef(object, ...)

# S3 method for class 'frmtmb_fit'
ranef(object, condVar = FALSE, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

- condVar:

  If `TRUE`, attach the conditional SDs of the modes (from the Laplace
  posterior) as a `"condSD"` attribute on each matrix, in matching
  layout.

## Value

A named list of levels-by-coefficients matrices, one per random-effect
term, KEYED BY THE GROUPING FACTOR as brms and lme4 key it (so
`ranef(fit)$g` and `coef(fit)$g` name the same group). Each matrix
carries its block label in a `"term"` attribute, which is also the key
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
uses and the `grp` column of
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html). That
long form (with a `condsd` column when `condVar = TRUE` was used) is
what broom.mixed-style code reads.

Several terms on ONE grouping factor give several entries under one
name, and `$` and `[[` then take the block label instead:
`ranef(fit)[["chi: 1 | id"]]`. The bare factor name is refused there,
and names the labels, rather than returning the first block. Positional
indexing reaches every block either way.

Three ordinary spellings put more than one block on one factor, so this
is not a nonlinear-model corner: an uncorrelated slope `(1 + x || g)`,
which desugars to `1 | g` and `0 + x | g`; two bars written out on the
same factor, `(1 | g) + (0 + x | g)`, the same two blocks; and a random
effect on more than one distributional or nonlinear parameter, such as
`bf(y ~ x + (1 | g), sigma ~ (1 | g))` (`1 | g` and `sigma: 1 | g`) or
`(1 | id)` on three nonlinear parameters. A CORRELATED slope
`(1 + x | g)` is one block and is addressed by the factor name as
before.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100, 1 + 0.5 * dd$x + rnorm(10, 0, 0.8)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)

# one matrix per random-effect term, levels by coefficients
ranef(fit)
#> $g   (1 | g)
#>    (Intercept)
#> 1   -0.7929115
#> 2    0.1938107
#> 3   -0.6968785
#> 4   -0.3665259
#> 5   -1.2649313
#> 6    1.0458808
#> 7    0.8996844
#> 8    0.2424397
#> 9   -0.2011950
#> 10   0.9405886
#> 
ranef(fit)$g[1:3, ]                 # keyed by the grouping factor
#>          1          2          3 
#> -0.7929115  0.1938107 -0.6968785 
attr(ranef(fit)$g, "term")          # the block it came from
#> [1] "1 | g"
ranef(fit)[["1 | g"]][1:3, ]        # or by that label
#>          1          2          3 
#> -0.7929115  0.1938107 -0.6968785 

# condVar adds the conditional SDs a caterpillar plot needs
re <- as.data.frame(ranef(fit, condVar = TRUE))
head(re)
#>     grp        term level    condval    condsd
#> 1 1 | g (Intercept)     1 -0.7929115 0.3795333
#> 2 1 | g (Intercept)     2  0.1938107 0.3736784
#> 3 1 | g (Intercept)     3 -0.6968785 0.3773146
#> 4 1 | g (Intercept)     4 -0.3665259 0.3788799
#> 5 1 | g (Intercept)     5 -1.2649313 0.3829424
#> 6 1 | g (Intercept)     6  1.0458808 0.3798946
with(re[order(re$condval), ],
     plot(condval, seq_along(condval), pch = 16,
          xlim = range(condval - 2 * condsd, condval + 2 * condsd),
          xlab = "conditional mode", ylab = "group"))
```
