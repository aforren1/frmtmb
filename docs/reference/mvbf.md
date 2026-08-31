# Combine formulas into a multivariate model

Each response keeps its own formula, family, dpar formulas, and addition
terms. Residual correlation between gaussian responses is requested with
`rescor = TRUE` or `set_rescor()`. Random-effect correlation across
responses uses the brms `|ID|` syntax, e.g. `(1 | p | g)` in several
formulas correlates their `g` effects.

## Usage

``` r
mvbf(..., rescor = FALSE)

set_rescor(rescor_value = TRUE)
```

## Arguments

- ...:

  [`bf()`](bf.md) formulas, each with a family attached (or supply one
  `family` to [`frm()`](frm.md) for all of them).

- rescor:

  Model residual correlation between the responses (gaussian only).

- rescor_value:

  For `set_rescor()`: turn residual correlation on or off.

## Value

An object of class `frmtmb_mvformula`.
