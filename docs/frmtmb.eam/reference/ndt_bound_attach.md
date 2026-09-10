# Put a non-decision-time bound on a family

Records `bound` on `fam`, sets the `ndt` link that makes the bound
structural, and, for a per-group bound, arranges for the per-row floor
to reach the family's likelihood as data. It does NOT touch the
likelihood: what a consumer's density has to do is three lines, written
below, and doing it there rather than through a wrapper is what keeps
this seam free of every assumption about where that density lives.

## Usage

``` r
ndt_bound_attach(fam, bound)
```

## Arguments

- fam:

  A `frmtmb_family` with a distributional parameter named `ndt`.

- bound:

  The result of
  [`ndt_bound()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound.md).

## Value

`fam`, modified.

## What the consumer's density must do

Call
[`ndt_apply()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_apply.md)
wherever it would have read `dpars[["ndt"]]`, in the density and in the
mean and in the simulator alike. A family that attaches a bound here and
then does NOT call it fits, converges and reports a non-decision time
measured at 37.5 percent too small, with nothing refusing;
[`?ndt_apply`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_apply.md)
has that measurement. Nothing on this side can check it. That is the
whole of it, and it is one line rather than three, because the three
lines of arithmetic come with a refusal that is longer than they are and
that is the part that matters. An earlier version of this seam
documented the arithmetic and left the refusal in prose; a consumer who
copied the first and not the second got a fraction silently read as
seconds, which is the defect the per-group bound exists to remove.

## What it changes on the family

`links$ndt`, `ndt_bound`, and for a per-group bound `aterm_data` and the
`ndt` starting value, which becomes the fraction one half. A scalar
bound leaves the family's own starting value alone. It also records the
family's pre-attach `aterm_data` and `ndt` starting value in a reserved
slot named `ndt_seam`, which is what makes a second call replace the
first bound instead of composing with it; a family with a slot of that
name of its own would lose it. Nothing else moves, and `st` is not
handled: a family with an across-trial range on the non-decision time
sets that bound itself.

Calling it twice replaces the first bound rather than stacking on it,
which matters because
[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
runs `family_finalize` again on the
[`influence()`](https://rdrr.io/r/stats/lm.influence.html),
`frm_simulate()` and prior-predictive paths.

## See also

[`ndt_bound()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound.md),
[`ndt_bound_of()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_of.md).

## Examples

``` r
# what a family_finalize() written against this seam looks like
finalize <- function(fam, y, aterms) {
  if (!is.null(ndt_bound_of(fam))) return(fam)
  ndt_bound_attach(fam, ndt_bound(y, aterms, what = "my_family"))
}
```
