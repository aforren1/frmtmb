# Add to a model formula

One `+` method serves
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) and
[`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md) objects
alike, so a sum of any number of formulas reads left to right. The sum
of `bf(y1 ~ x)`, `bf(y2 ~ x)` and `bf(y3 ~ x)` is
`mvbf(bf(y1 ~ x), bf(y2 ~ x), bf(y3 ~ x))`.

## Usage

``` r
# S3 method for class 'frmtmb_bform'
e1 + e2
```

## Arguments

- e1:

  A [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) or
  [`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md)
  object.

- e2:

  The object to add.

## Value

A `frmtmb_formula` or a `frmtmb_mvformula`.

## Details

The right-hand side can be another formula, a family,
[`lf()`](https://aforren1.github.io/frmtmb/reference/lf.md),
[`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) or
[`set_rescor()`](https://aforren1.github.io/frmtmb/reference/mvbf.md). A
family added to a single formula sets its family. A family added to a
multivariate formula goes to every response that has no family yet, so
in the sum of `bf(o ~ x)`,
[`cumulative()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
`bf(y ~ x)` and [`gaussian()`](https://rdrr.io/r/stats/family.html) the
response `o` stays ordinal. brms instead gives the last family to every
response. An [`lf()`](https://aforren1.github.io/frmtmb/reference/lf.md)
or [`nlf()`](https://aforren1.github.io/frmtmb/reference/nlf.md) added
to a multivariate formula names its response with `resp =`.

## Examples

``` r
# three responses, one family for all of them
bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + gaussian()
#> Warning: Incompatible methods ("+.frmtmb_mvformula", "+.frmtmb_formula") for "+"
#> Error in bf(y1 ~ x) + bf(y2 ~ x) + bf(y3 ~ x): non-numeric argument to binary operator

# a family per response, and a dpar formula for the third response
bf(o ~ x) + cumulative() + bf(y1 ~ x) + gaussian() + bf(y2 ~ x) +
  lf(sigma ~ x, resp = "y2")
#> Warning: Incompatible methods ("+.frmtmb_mvformula", "+.frmtmb_formula") for "+"
#> Error in bf(o ~ x) + cumulative() + bf(y1 ~ x) + gaussian() + bf(y2 ~     x): non-numeric argument to binary operator
```
