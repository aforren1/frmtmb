# Set up a model formula

Specifies a model with brms-compatible syntax. Distributional parameters
(dpars) can get their own formulas with the full predictor grammar, or
be fixed to constants: `bf(y ~ x + (1 | g), sigma ~ z + (1 | g))`,
`bf(y ~ x, sigma = 1)`. Nonlinear (`nl = TRUE`) and multivariate
formulas arrive in later versions and signal an error for now.

## Usage

``` r
bf(formula, ..., family = NULL, nl = FALSE)
```

## Arguments

- formula:

  The model formula for `mu`.

- ...:

  Two-sided formulas for other dpars (the left-hand side names the dpar,
  e.g. `sigma ~ z`), or named scalars fixing a dpar to a constant on the
  response scale (e.g. `sigma = 1`).

- family:

  Optional family; can also be attached with `+` or passed to
  [`frm()`](frm.md).

- nl:

  Nonlinear-formula flag (not yet supported).

## Value

An object of class `frmtmb_formula`.

## Details

The left-hand side accepts addition terms after `|`:
`y | weights(w) ~ ...` and `y | trials(n) ~ ...`. Every linear predictor
accepts lme4-style random effects `(1 | g)`, `(1 + x | g)`, `(x || g)`,
and explicit covariance-structure wrappers `us(x | g)` and
`diag(x | g)`.

Attach a family with `+`, for example `bf(y ~ x) + gaussian()`.
