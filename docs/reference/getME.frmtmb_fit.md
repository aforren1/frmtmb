# Extract components of a fit, lme4 style

A small [`lme4::getME()`](https://rdrr.io/pkg/lme4/man/getME.html)
vocabulary, for downstream code written against merMod objects.
Registered on lme4's generic, so call it as `lme4::getME(fit, "X")` or
load lme4 first.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
getME(object, name, resp = NULL, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- name:

  One or more names from the vocabulary above. A vector returns a named
  list.

- resp:

  Response name, for the design extractors on a multivariate fit.

- ...:

  Unused.

## Value

The requested component, or a named list when `name` names several.

## Details

Supported names:

- `"X"`:

  Fixed-effect design matrix of the `mu` predictor.

- `"Z"`, `"Zt"`:

  The sparse random-effect design of the `mu` predictor and its
  transpose. Columns (rows of `Zt`) span the whole random-effect
  coefficient vector, so a block belonging to another distributional
  parameter contributes zero columns here.

- `"beta"`, `"fixef"`:

  The primary (`mu`-family) fixed-effect coefficients, named.
  Coefficients of auxiliary distributional parameters are a separate
  vector; use
  [`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) for
  all of them.

- `"b"`:

  Conditional modes in coefficient space, aligned with the columns of
  `Z`. Reduced-rank (`rr()`) blocks are expanded through their loadings,
  so this is not the internal parameter vector.

- `"theta"`:

  Covariance parameters on the internal (unconstrained) scale, as in
  [`confint()`](https://rdrr.io/r/stats/confint.html). These are not
  lme4's relative-covariance-factor entries.

- `"lower"`:

  Lower bounds on `theta`. The internal parameterization is unbounded,
  so this is a vector of `-Inf`, not lme4's mixture of `0` and `-Inf`.
  Code that tests `theta == lower` to detect a singular fit will never
  fire; use
  [`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
  or
  [`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
  instead.

- `"sigma"`:

  Residual standard deviation
  ([`sigma.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/sigma.frmtmb_fit.md)).

- `"flist"`:

  The grouping factors, one per distinct grouping variable. Smooth and
  Gaussian-process blocks are excluded: their levels are basis
  functions, not groups. There is no `"assign"` attribute.

- `"n_rtrms"`, `"n_rfacs"`:

  Number of random-effect terms and of distinct grouping factors.

Multivariate fits have one design per response, so `"X"`, `"Z"` and
`"Zt"` need `resp`; the other names answer without it.
