# Drift specifications for [`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md)

The drift is the mean rate of evidence accumulation, `a(x, t)`. A GDDM
lets it depend on the current state `x`, on time, and on a covariate.
Terms are summed: pass one, or a list of several.

## Usage

``` r
gddm_drift_constant()

gddm_drift_coherence(cmax = 1, cov = 1L)

gddm_drift_leak()
```

## Arguments

- cmax:

  Scale the coherence covariate is divided by, so that the nonlinearity
  is anchored at a coherence a subject sees. The paper uses the largest
  coherence in the design.

- cov:

  Which `vreal()` value carries the coherence, counting from 1.

## Value

A `gddm_component`.

## Details

- `gddm_drift_constant()`:

  `a = mu`. The textbook drift-diffusion drift. Free parameter `mu`,
  identity link, so it is signed.

- `gddm_drift_coherence()`:

  `a = sign(C) mu (|C| / cmax)^alpha`, the coherence nonlinearity of
  Shinn et al. (2020). The covariate `C` is supplied through `vreal()`
  and `cmax` scales it. Free parameters `mu` (identity link) and `alpha`
  (log link).

- `gddm_drift_leak()`:

  `a = -leak x`, added to whichever base term is used. Positive `leak`
  is leaky integration, which pulls the accumulator back toward zero;
  negative `leak` is unstable integration, which pushes it away.
  Identity link. Note the sign is the paper's `l` and the negative of
  PyDDM's `leak`.

Exactly one base term, `gddm_drift_constant()` or
`gddm_drift_coherence()`, must be present and must come first, because
it supplies `mu`, the parameter that receives the model formula.

## Writing your own

A drift term is the value of
[`gddm_drift_term()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_drift_term.md).
Its `fn` is called on the tape as `fn(x, t, p, cov)` where `x` is the
accumulator state at the grid nodes in the original coordinate, `t` is
the time as a plain number, `p` is a named list holding one value of
each free parameter for the condition being solved, and `cov` is that
condition's `vreal()` values. It returns the drift at each node, and
must contain no comparison against a parameter.

## See also

[`gddm()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm.md),
[`gddm_bound_constant()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md),
[`gddm_start_point()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-start.md)
