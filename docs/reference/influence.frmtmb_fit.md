# Influence measures by case deletion

Refits the model with one group (or observation) left out at a time,
warm-started at the full-data estimates, and collects the fixed-effect
and covariance-parameter changes.
[`cooks.distance()`](https://rdrr.io/r/stats/influence.measures.html) on
the result gives the scaled fixed-effect displacement (calling it on the
fit itself runs
[`influence()`](https://rdrr.io/r/stats/lm.influence.html) first);
[`dfbeta()`](https://rdrr.io/r/stats/influence.measures.html) and
[`dfbetas()`](https://rdrr.io/r/stats/influence.measures.html) give the
per-unit coefficient changes, raw and scaled by the coefficient standard
errors (the lme4 influence surface).
[`plot.frmtmb_influence()`](https://aforren1.github.io/frmtmb/reference/plot.frmtmb_influence.md)
draws all three.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
influence(model, groups = NULL, data = NULL, force = FALSE, ...)

# S3 method for class 'frmtmb_fit'
cooks.distance(model, ...)

# S3 method for class 'frmtmb_influence'
dfbeta(model, ...)

# S3 method for class 'frmtmb_influence'
dfbetas(model, ...)
```

## Arguments

- model:

  A `frmtmb_fit`.

- groups:

  Name of a random-effect grouping factor (see
  [`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.md)) to
  delete level-wise; `NULL` deletes single observations (refuses for
  large n unless `force = TRUE`).

- data:

  The original model data; defaults to the stored model frame, which
  works unless the formula uses variables that are not stored raw (e.g.
  inside [`poly()`](https://rdrr.io/r/stats/poly.html)).

- force:

  Allow observation-wise deletion for n \> 500.

- ...:

  Unused.

## Value

A `frmtmb_influence` object: `fixed` and `theta` matrices (one row per
deleted unit) plus the full-data reference.

## Examples

``` r
# \donttest{
set.seed(7)
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
infl <- influence(fit, groups = "g")
cooks.distance(infl)
#>          1          2          3          4          5          6 
#> 0.04211219 0.49008216 0.04168487 0.12687101 0.02619522 0.18339080 
# }
```
