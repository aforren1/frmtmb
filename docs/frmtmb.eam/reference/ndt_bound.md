# The bound a non-decision time is measured against

Derives the upper bound
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
and its siblings put on the non-decision time, so that a family in
another package gets the same bound, the same `ndt_group()` behavior and
the same refusals without writing them again. Pair it with
[`ndt_bound_attach()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_attach.md),
which puts the result on a family object.

## Usage

``` r
ndt_bound(y, aterms = list(), max_ndt = NULL, what = "ndt_bound()")
```

## Arguments

- y:

  The response, a numeric vector of times.

- aterms:

  The addition-term values for this response, as a density receives
  them. Only `ndt_group` is read.

- max_ndt:

  An absolute bound in the response's own units, or `NULL` to take the
  fastest response. Above the fastest response it is refused, and it
  cannot be combined with `ndt_group()`.

- what:

  The family's name, used in the refusals so that the sentence a user
  reads names the family they wrote.

## Value

An object of class `"frmtmb_eam_ndt_bound"`: a list with `ub`, the
scalar bound; `floors`, the per-group bounds or `NULL`; `sizes`, the
trial count behind each one; `pending`, always `FALSE` here; and `what`.

## Details

The bound is one of two things.

- With no `ndt_group()` in `aterms` it is ONE number: `max_ndt` when
  that is given, and the fastest response otherwise. `ndt` is then a
  TIME, the link carries the bound, and a consumer's density reads
  `dpars$ndt` in the response's own units.

- With `ndt_group()` it is one number PER GROUP, that group's own
  fastest response. `ndt` is then a FRACTION of the row's own bound, on
  a plain logit, and the density multiplies it back out by the per-row
  vector that arrives in the addition-term values under the reserved
  name `ndt_floor`.

`max_ndt` and `ndt_group()` together are refused: they set the same
bound to different things.

## Why the group's own floor and not the data set's

The information about one subject's non-decision time is that subject's
own fastest response. At 30 subjects by 400 trials with a
between-subject spread of 26 ms on a mean of 250 ms, 20 of the 30
subjects have a true non-decision time above the GLOBAL minimum while
none is above its own, so a single bound runs the fit to the wall.
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
and `dev/ndt-findings.md` carry the measurements.

## See also

[`ndt_bound_attach()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_attach.md)
to put it on a family,
[`ndt_bound_of()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_of.md)
to read one back,
[`ndt_time()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_time.md)
to report a fitted non-decision time in the response's own units.

## Examples

``` r
rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
g <- c("a", "a", "a", "b", "b")
# `[[` and not `$`: the record has five names and `bd$s` would
# partial-match `sizes`
ndt_bound(rt)[["ub"]]
#> [1] 0.31
ndt_bound(rt, list(ndt_group = ndt_bound_key(g)))[["floors"]]
#> 1627390186 1644167403 
#>       0.31       0.61 
```
