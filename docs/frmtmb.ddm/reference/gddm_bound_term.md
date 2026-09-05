# Build a boundary term

The extension point behind
[`gddm_bound_constant()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md)
and its siblings.

## Usage

``` r
gddm_bound_term(label, dpars, fn)
```

## Arguments

- label:

  One word naming the term, used when the family prints.

- dpars:

  Named list, one entry per free parameter, each a list with `link` and
  `init`; see
  [`gddm_drift_term()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_drift_term.md).

- fn:

  `function(t, p, ctl)` returning a list with `B` and `dlogB`.

## Value

A `gddm_component`.

## See also

[gddm-bound](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md)
