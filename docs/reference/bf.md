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
bf(formula, ..., family = NULL, nl = NULL, center = NULL)
```

## Arguments

- formula:

  The model formula for `mu`, or a formula `bf()` already built. Given
  one of those and nothing else, `bf()` returns it unchanged, as brms's
  `bf()` does; further arguments add to it.

- ...:

  Two-sided formulas for other dpars (the left-hand side names the dpar,
  e.g. `sigma ~ z`, or several sharing one right-hand side, e.g.
  `b1 + b2 ~ 1`), one-sided formulas named by their dpar (`sigma = ~ z`,
  the same formula), or named scalars fixing a dpar to a constant on the
  response scale (e.g. `sigma = 1`).

- family:

  Optional family; can also be attached with `+` or passed to
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md), which
  uses [`gaussian()`](https://rdrr.io/r/stats/family.html) when nothing
  names one.

- nl:

  Nonlinear-formula flag: the main formula becomes a nonlinear
  expression of named parameters, each given its own `...` formula with
  the full predictor grammar. `NULL`, the default as in brms, means
  `FALSE` for a new formula and keeps the setting of a formula `bf()`
  already built.

- center:

  Whether the location formula's intercept is brms's class
  `"Intercept"`, a density on the intercept at the means of the
  predictors. `NULL`, the default, means `TRUE` for a new formula and
  keeps the setting of a formula `bf()` already built. `FALSE` makes the
  intercept an ordinary coefficient, class `"b"` with coef
  `"Intercept"`, as `0 + Intercept` in the formula does; see the section
  on brms's reserved `Intercept` below. It changes priors only: the
  likelihood and the maximum likelihood fit are the same. It applies to
  the location formula alone, as in brms; give a parameter formula its
  own with [`lf()`](https://aforren1.github.io/frmtmb/reference/lf.md).

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

## brms's reserved `Intercept`

In a formula without an intercept, `Intercept` is a reserved name, as in
brms: `y ~ 0 + Intercept + x` is the model `y ~ 1 + x`, with the
intercept as an ordinary population-level coefficient. Factors get
treatment contrasts, as they do beside an intercept. The likelihood is
the same, so a maximum likelihood fit is the same; what changes is the
prior. brms places a class `"Intercept"` prior at the means of the
predictors, and a class `"b"` prior on this coefficient at zero, and so
does frmtmb: `set_prior(..., class = "b")` reaches it,
`class = "Intercept"` does not. `bf(y ~ x, center = FALSE)` gives the
same model. The spelling works in every linear formula: the location, a
distributional parameter (`sigma ~ 0 + Intercept + z`) and a nonlinear
parameter (`a ~ 0 + Intercept + x`, where it changes nothing, because
brms never centers a nonlinear parameter).

`Intercept` must be a term of its own; `Intercept:x` is refused. An
ordinal family refuses `0 + Intercept`, as brms does, because its
thresholds take the intercept's place. A data column named `Intercept`
must hold only ones in such a model, as brms requires; in a formula with
an intercept, `Intercept` is an ordinary variable. Prediction needs no
`Intercept` column in `newdata`.

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
#> y ~ a * exp(-b * x) (nonlinear)
#> a ~ 1 
#> b ~ 1 + (1 | g) 
# an intercept that a class "b" prior reaches; the two are one model
bf(y ~ 0 + Intercept + x)
#> y ~ 0 + Intercept + x
bf(y ~ x, center = FALSE)
#> Error: Cannot interpret bf() argument 'center': expected a dpar formula or a named numeric constant
```
