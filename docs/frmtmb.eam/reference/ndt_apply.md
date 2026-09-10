# The non-decision time on the response's own scale

The whole of what a per-group bound costs a consumer's density, and it
is one multiplication and one refusal. Call it wherever the density, the
mean or the simulator would have read `dpars$ndt`.

## Usage

``` r
ndt_apply(dpars, aterms = dpars, what = "this family")
```

## Arguments

- dpars:

  The distributional parameters, as the density receives them. Only
  `ndt` is read.

- aterms:

  The addition-term values. Only `ndt_floor` and `ndt_group` are read.
  Defaults to `dpars`, for a family whose engine merges the two into one
  list, which is what
  [`frmtmb.learn::rlddm()`](https://aforren1.github.io/frmtmb/frmtmb.learn/reference/rlddm.html)
  does.

- what:

  The family's name, so that the refusal names the family the user
  wrote.

## Value

The non-decision time in the response's own units, recycled to the
length the arguments imply. An 'RTMB' advector passes through: there is
no comparison and no branch on a parameter here, only two
[`is.null()`](https://rdrr.io/r/base/NULL.html) tests on data.

## Details

With a per-group bound `dpars$ndt` is a FRACTION and the row's own bound
arrives in the addition-term values as `ndt_floor`, so this returns
`ndt * ndt_floor`. With a scalar bound `ndt_floor` is absent,
`dpars$ndt` is already a time, and this returns it untouched, so a model
that did not opt in is the arithmetic it always was.

## Why this is exported rather than documented

The arithmetic is three lines. The REFUSAL is not, and it is the part
that matters: a grouped model whose density is reached without
`ndt_floor` would read the fraction as seconds and report a converged
fit at a non-decision time several times too small. That is the silent
wrong answer the per-group bound exists to remove, and a consumer who
copies the arithmetic and not the refusal reintroduces it. The grouping
itself is in the addition-term values on every path that can reach a
likelihood, so `ndt_group` present with `ndt_floor` absent is exactly
that condition, and this refuses on it.

## What happens if a consumer does not call this

Nothing refuses, and the fit is silently wrong. This package cannot
check another package's arithmetic: a family that attaches a per-group
bound through
[`ndt_bound_attach()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_attach.md)
and then reads `dpars[["ndt"]]` directly is reading a fraction as a
time, and everything downstream still looks right.

Measured, on two families identical in every line but this one, a
shifted gamma over four groups of 300 rows with true shifts 0.60 to
0.90: both fit, both report convergence code 0, a positive definite
Hessian and no non-finite standard errors, both give the IDENTICAL
log-likelihood 1345.58425, and the one that skipped this call has an
[`ndt_time()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_time.md)
up to **37.5 percent too small**. Every fitted time is still below its
own group's floor on both, which is why that check cannot catch it: it
tests the floor lookup, not the arithmetic.

The fraction one half that
[`ndt_bound_attach()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_attach.md)
sets as the `ndt` starting value is a partial guard and only by luck. On
a faster design the same broken family died at once with
`NA/NaN gradient evaluation`, because half a second is past most of the
data; on a slower one it converged silently. A guard that works only
when the responses are fast is not a guard.

## What it does NOT do

It returns a value rather than a modified `dpars`, because where the
value goes is the consumer's own layout. Two consequences a consumer
owns. A family with an across-trial RANGE on the non-decision time
scales that separately; this package's own families use twice the bound.
And a family that reads a linear predictor back through frmtmb's robust
dpar accessors must drop `dpars$.eta_ndt` itself after calling this,
because the linear predictor beside `ndt` no longer maps to the value
that was used.

## See also

[`ndt_bound_attach()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/ndt_bound_attach.md),
which puts `ndt_floor` where this reads it.

## Examples

``` r
# no grouping: `ndt` is already a time
ndt_apply(list(ndt = 0.2))
#> [1] 0.2
# grouped: a fraction of each row's own bound
ndt_apply(list(ndt = 0.5), list(ndt_floor = c(0.4, 0.6)))
#> [1] 0.2 0.3
```
