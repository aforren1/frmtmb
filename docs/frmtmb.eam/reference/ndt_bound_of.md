# The non-decision-time bound a family already carries

`NULL` when the family has none, and `NULL` when the one it has is still
PENDING, which is the state a family constructed outside
[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
is in before the response has been seen.

## Usage

``` r
ndt_bound_of(fam)
```

## Arguments

- fam:

  A `frmtmb_family` object, fitted or not.

## Value

An object of class `"frmtmb_eam_ndt_bound"`, or `NULL`.

## Details

A settled bound is KEPT rather than re-derived, and that is what this is
for.
[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
runs `family_finalize` again on new rows for
[`influence()`](https://rdrr.io/r/stats/lm.influence.html),
`frm_simulate()` and the prior-predictive path, so a bound re-derived
from fewer rows would make a leave-one-out refit a refit of a DIFFERENT
model. A `family_finalize` written against this seam therefore asks
first and returns the family untouched when the answer is not `NULL`.

## See also

[`ndt_bound()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound.md),
[`ndt_bound_attach()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_attach.md).

## Examples

``` r
# a family built outside frm() has no settled bound yet
ndt_bound_of(wiener())
#> NULL
ndt_bound_of(wiener(max_ndt = 0.4))[["ub"]]
#> [1] 0.4
```
