# An empty prior specification

brms's `empty_prior()`: a prior object with no specifications, to add
specifications to with `+` or [`c()`](https://rdrr.io/r/base/c.html).
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) given it
applies no prior.

## Usage

``` r
empty_prior()
```

## Value

A `frmtmb_priorlist` of length zero.

## Examples

``` r
pr <- empty_prior()
pr <- pr + set_prior("normal(0, 1)", class = "b")
pr
#> b ~ normal(0, 1)
```
