# Starting-point specifications for [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)

Where the accumulator starts, as a distribution over the rescaled state
at time zero. `bias` is the relative start point in `(0, 1)`, the
fraction of the boundary separation above the lower boundary, exactly as
in
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md);
0.5 is unbiased.

## Usage

``` r
gddm_start_point()

gddm_start_uniform()
```

## Value

A `gddm_component`.

## Details

- `gddm_start_point()`:

  All mass at `bias`. Free parameter `bias`, logit link.

- `gddm_start_uniform()`:

  Uniform on an interval of half-width `sz` around `bias`, with `sz`
  measured as a fraction of the half separation. Free parameters `bias`
  and `sz`, both logit link. This is the tape-safe counterpart of the
  Ratcliff start-point variability, which the analytic family can only
  reach by quadrature.

A point start is a delta function, which no fixed grid can hold. Both
components spread the initial mass with the same cubic B-spline used for
the non-decision-time shift, so the starting distribution is smooth in
`bias` and no grid index depends on a parameter. Both are then
normalized to unit mass on the grid, which also keeps the family honest
when `bias` and `sz` together push mass past a boundary.

## Writing your own

A start term is the value of
[`gddm_start_term()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_start_term.md).
Its `fn` is called as `fn(y, h, p)`, where `y` are the rescaled grid
nodes on `(-1, 1)` and `h` is the node spacing, and it returns a density
at those nodes.

## See also

[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md),
[gddm-drift](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md),
[gddm-bound](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md)
