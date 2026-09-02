# Set up a model formula

Specifies a model with brms-compatible syntax. Distributional parameters
(dpars) can get their own formulas with the full predictor grammar, or
be fixed to constants: `bf(y ~ x + (1 | g), sigma ~ z + (1 | g))`,
`bf(y ~ x, sigma = 1)`. Nonlinear formulas (`nl = TRUE`) and
multivariate models (see
[`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md)) use the
same grammar.

## Usage

``` r
bf(formula, ..., family = NULL, nl = FALSE)
```

## Arguments

- formula:

  The model formula for `mu`.

- ...:

  Two-sided formulas for other dpars (the left-hand side names the dpar,
  e.g. `sigma ~ z`, or several sharing one right-hand side, e.g.
  `b1 + b2 ~ 1`), or named scalars fixing a dpar to a constant on the
  response scale (e.g. `sigma = 1`).

- family:

  Optional family; can also be attached with `+` or passed to
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md), which
  uses [`gaussian()`](https://rdrr.io/r/stats/family.html) when nothing
  names one.

- nl:

  Nonlinear-formula flag: the main formula becomes a nonlinear
  expression of named parameters, each given its own `...` formula with
  the full predictor grammar.

## Value

An object of class `frmtmb_formula`.

## Details

The left-hand side accepts addition terms after `|`:
`y | weights(w) ~ ...` and `y | trials(n) ~ ...`. Every linear predictor
accepts lme4-style random effects `(1 | g)`, `(1 + x | g)`, `(x || g)`,
and explicit covariance-structure wrappers `us(x | g)` and
`diag(x | g)`.

Attach a family with `+`, for example `bf(y ~ x) + gaussian()`, or pass
one to [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md). A
model that names no family is gaussian.

## Examples

``` r
# brms-style model formulas: attach a family with `+`
bf(y ~ x + (1 | g)) + gaussian()
#> y ~ x + (1 | g) 
#> Family: gaussian 
# distributional parameters get their own formulas or constants
bf(y ~ x, sigma ~ x)
#> y ~ x 
#> sigma ~ x 
bf(y ~ x, shape = 2) + Gamma()
#> y ~ x 
#> shape = 2 
#> Family: Gamma 
# nonlinear models declare parameter formulas and nl = TRUE
bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1 + (1 | g), nl = TRUE)
#> y ~ a * exp(-b * x) 
#> a ~ 1 
#> b ~ 1 + (1 | g) 
```
