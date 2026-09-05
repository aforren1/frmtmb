# Build a starting-point term

The extension point behind
[`gddm_start_point()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-start.md)
and its sibling.

## Usage

``` r
gddm_start_term(label, dpars, fn)
```

## Arguments

- label:

  One word naming the term, used when the family prints.

- dpars:

  Named list, one entry per free parameter, each a list with `link` and
  `init`; see
  [`gddm_drift_term()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_drift_term.md).

- fn:

  `function(y, h, p)` returning the starting density at the rescaled
  grid nodes.

## Value

A `gddm_component`.

## See also

[gddm-start](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-start.md)
