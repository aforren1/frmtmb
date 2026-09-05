# Build a drift term

The extension point behind
[`gddm_drift_constant()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md)
and its siblings.

## Usage

``` r
gddm_drift_term(label, dpars, fn, aterms = character(0), base = FALSE)
```

## Arguments

- label:

  One word naming the term, used when the family prints.

- dpars:

  Named list, one entry per free parameter, each a list with `link` (a
  link name or a link object) and `init`, a `function(y, aterms)`
  returning a starting value on the natural scale.

- fn:

  `function(x, t, p, cov)` returning the drift at each node.

- aterms:

  Character vector of addition-term values `fn` reads, named as they
  reach the density: `"vreal1"`.

- base:

  `TRUE` if the term supplies `mu` and can stand alone.

## Value

A `gddm_component`.

## See also

[gddm-drift](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md)
