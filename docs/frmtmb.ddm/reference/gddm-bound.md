# Boundary specifications for [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)

The decision boundaries sit at plus and minus `B(t)`, so the separation
between them is `2 B(t)`. All three components below name the separation
`bs`, as
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
does, and set `B(0) = bs / 2`, so a boundary separation estimated here
is directly comparable with one from the analytic family.

## Usage

``` r
gddm_bound_constant()

gddm_bound_exponential()

gddm_bound_linear()
```

## Value

A `gddm_component`.

## Details

- `gddm_bound_constant()`:

  `B(t) = bs / 2`. Free parameter `bs`, log link.

- `gddm_bound_exponential()`:

  `B(t) = (bs / 2) exp(-t / tau)`. Free parameters `bs` and the time
  constant `tau`, both log link. A large `tau` is a bound that barely
  collapses, which is why it is where the fit starts. Note PyDDM
  parameterizes the same bound by a rate, the reciprocal of `tau`.

- `gddm_bound_linear()`:

  `B(t) = (bs / 2) (1 - kappa t / t_max)`, where `kappa` is the fraction
  of the bound lost by the end of the modeled window. Free parameters
  `bs` (log link) and `kappa` (logit link). The usual spelling of a
  linear collapse clips the bound at zero, which is a comparison against
  a parameter and cannot be taped; bounding the fraction below one
  instead keeps the boundary strictly positive across the window by
  construction, with no clipping.

## Writing your own

A bound term is the value of
[`gddm_bound_term()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_bound_term.md).
Its `fn` is called on the tape as `fn(t, p, ctl)` for a plain number
`t`, and returns a list with `B`, the boundary at that time, and
`dlogB`, its logarithmic derivative `B'(t) / B(t)`. The solver needs the
second because the change of variable that pins the walls contributes a
`-y B'(t)/B(t)` term to the drift; supplying it rather than differencing
`B` keeps the rescaled drift exact.

## See also

[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md),
[gddm-drift](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md),
[`gddm_start_point()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-start.md)
