# Add parameter formulas to a model formula

brms's `lf()`: one or more two-sided formulas for distributional or
nonlinear parameters, added to a
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) with `+`. It
is sugar for passing the same formulas to
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) directly,
and is useful when the parameter formulas are built somewhere else than
the response formula.

## Usage

``` r
lf(..., resp = NULL, center = NULL)
```

## Arguments

- ...:

  Two-sided formulas naming the parameter on the left, e.g. `sigma ~ x`
  or (with `nl = TRUE` on the
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md)) a
  nonlinear parameter's formula `a ~ 1 + (1 | g)`, or one-sided formulas
  named by their parameter, `sigma = ~ x`.

- resp:

  The response the formulas belong to, when the `lf()` is added to a
  multivariate formula. `NULL` (the default) adds them to the
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) on the
  left of the `+`, which must then be a single formula.

- center:

  `FALSE` makes the intercept of each of these formulas an ordinary
  coefficient, class `"b"` with coef `"Intercept"`, instead of brms's
  class `"Intercept"`: the same as `0 + Intercept` in the formula. See
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md). `NULL`,
  the default, leaves the intercept as class `"Intercept"`.

## Value

An object of class `frmtmb_lf`, to be added to a
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md).

## Details

`bf(y ~ x) + lf(sigma ~ z)` and `bf(y ~ x, sigma ~ z)` give the same
model. In a multivariate model, add an `lf()` to the
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) of the
response it belongs to, or name that response with `resp =`:
`bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + lf(sigma ~ z, resp = "y3")`.

## Examples

``` r
# the two spellings are the same model
bf(y ~ x) + lf(sigma ~ z)
#> y ~ x
#> sigma ~ z 
bf(y ~ x, sigma ~ z)
#> y ~ x
#> sigma ~ z 

# nonlinear parameter formulas can arrive the same way
bf(y ~ a * exp(-b * x), a ~ 1, nl = TRUE) + lf(b ~ 1 + (1 | g))
#> y ~ a * exp(-b * x) (nonlinear)
#> a ~ 1 
#> b ~ 1 + (1 | g) 

# in a multivariate formula, resp = says which response it modifies
bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + lf(sigma ~ z, resp = "y3")
#> y1 ~ x
#> y2 ~ x
#> y3 ~ x
#> sigma ~ z 
#> rescor: FALSE 
```
