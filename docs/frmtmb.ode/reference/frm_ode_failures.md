# Groups whose ODE solve failed in the last `frm_ode()` call

Reports the groups that
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
could not solve on its most recent call, and whose rows therefore hold
the `penalty` value rather than a solution.
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
also warns when that happens; this function is for reading the record
afterwards, for example after a
[`predict()`](https://rdrr.io/r/stats/predict.html) whose warnings were
suppressed.

## Usage

``` r
frm_ode_failures()
```

## Value

`NULL` when the last call solved every group. Otherwise a list with
`groups` (the labels that failed), `n_groups` (how many groups the call
had), `penalty` (the value written into their rows), and `when`.

## Details

The record covers the last call only, from anywhere: a fit, a
[`predict()`](https://rdrr.io/r/stats/predict.html), or a direct call.
It is reset by every
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
call, so a clean call clears it. See the "Failed solves" section of
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
for why a failure during fitting usually cannot be recorded at all.

## See also

[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)

## Examples

``` r
# NULL until a solve fails
frm_ode_failures()
#> NULL
```
