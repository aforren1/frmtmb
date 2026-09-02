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
lf(...)
```

## Arguments

- ...:

  Two-sided formulas naming the parameter on the left, e.g. `sigma ~ x`
  or (with `nl = TRUE` on the
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md)) a
  nonlinear parameter's formula `a ~ 1 + (1 | g)`.

## Value

An object of class `frmtmb_lf`, to be added to a
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md).

## Details

`bf(y ~ x) + lf(sigma ~ z)` and `bf(y ~ x, sigma ~ z)` give the same
model. In a multivariate model an `lf()` must be added to the
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) of the
response it belongs to, before the responses are combined.

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
```
