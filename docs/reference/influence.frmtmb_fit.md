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

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

## Value

A `frmtmb_influence` object: `fixed` and `theta` matrices (one row per
deleted unit) plus the full-data reference.

## Details

Each deletion refits the same model, so an ordinal response keeps the
threshold count of the full-data fit, and the per-level counts of
`thres(gr = )`, whether the response is coded as integers, as a
character vector or as an ordered factor. Deleting the last observation
in ANY category, bottom, interior or top, leaves that category's
threshold in the model and unidentified, which shows as a large
displacement in its column, never as a shorter row or a coefficient
under the next column's name.

Two deletions cannot refit the fitted model at all and are refused by
name. A `groups = ` deletion that removes a whole level of
`thres(gr = )` takes away every observation behind that level's
thresholds, so they could only be invented. `data = ` holding a response
category outside the fitted model's threshold layout, or a
`thres(gr = )` level the fit never saw, describes a different model. A
category the fit never OBSERVED is not outside that layout when
`thres(K)` declared it, so `data = ` may reach one of those.

A refit that fails for any reason, these two included, leaves its row
`NA`, and the failures are counted. When some units failed the table
comes back with a warning carrying the count and the first reason; when
every unit failed there is no table to hand back, so that first reason
is raised as an error instead of a silent matrix of `NA`.

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
