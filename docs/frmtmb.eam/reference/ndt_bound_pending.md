# The bound a family carries before it has seen a response

A family object built outside
[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
has no data, so it has no fastest response and no bound. That state is
not neutral, and this is the seam entry that settles it, so that a
family which reaches
[`ndt_bound_attach()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_attach.md)
at CONSTRUCTION is in one of two honest states rather than carrying
whatever link its constructor happened to declare.

## Usage

``` r
ndt_bound_pending(max_ndt = NULL, what = "ndt_bound()")
```

## Arguments

- max_ndt:

  An absolute bound in the response's own units, or `NULL` for the
  pending state.

- what:

  The family's name, for the refusals.

## Value

An object of class `"frmtmb_eam_ndt_bound"`.

## Details

- `max_ndt = NULL` gives a PENDING bound. The `ndt` link then refuses by
  name until
  [`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
  settles it. That refusal is load-bearing rather than decorative:
  [`frmtmb::mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.html)
  never finalizes its components, so a component with an unsettled bound
  would read `ndt` on a scale nothing had set; and `bf(ndt = )` runs
  `linkfun()` at PARSE time, where the range check is against whichever
  link the family is carrying, so a constant outside the SETTLED bound
  passes a check against a placeholder and reaches the objective as
  `NaN`.

- `max_ndt = <a number>` gives a SETTLED bound, and the family scores
  outside `frm()`: `ndt` is a time measured against that bound,
  `predict(dpar = "ndt", type = "response")` reports seconds, and a
  pinned constant means what it says.

The second is what a package reaching into another family's density slot
needs, and it is what
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md),
[`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
and
[`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
do for themselves through their own `max_ndt` argument.

## See also

[`ndt_bound()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound.md)
for the bound derived once the response is in hand,
[`ndt_bound_attach()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_attach.md)
to put either on a family.

## Examples

``` r
ndt_bound_pending()[["pending"]]
#> [1] TRUE
ndt_bound_pending(0.4)[["ub"]]
#> [1] 0.4
```
