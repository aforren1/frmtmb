# Check a custom family's log-density for AD safety

Tapes the family's `lpdf` on test values and compares the AD gradient
against central finite differences. A mismatch usually means the lpdf
uses operations the tape cannot see (base
[`matrix()`](https://rdrr.io/r/base/matrix.html)/[`c()`](https://rdrr.io/r/base/c.html)
on advectors, branching on parameter values, `min`/`max`, clamping).

## Usage

``` r
check_custom_family(family, y, dpars, aterms = list(), tol = 1e-04)
```

## Arguments

- family:

  A `frmtmb_family` (from
  [`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
  /
  [`custom_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)).

- y:

  A response vector of test data.

- dpars:

  Named list of numeric test values, one entry per dpar (each of length
  1 or `length(y)`).

- aterms:

  Named list of addition-term values (e.g. `trials`).

- tol:

  Maximum relative gradient error.

## Value

Invisibly `TRUE`; signals an error on failure.
