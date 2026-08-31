# Define a model family

Constructs a family object for use with [`frm()`](frm.md). The
log-density function must be vectorized and AD-compatible: it is
evaluated on RTMB 'advector' objects during taping, so it must use
RTMB-overloaded operations (RTMB and RTMBdist `d*` functions, plain
arithmetic) and must not branch on parameter values.

## Usage

``` r
frmtmb_family(
  family,
  dpars,
  links,
  lpdf,
  valid_y = NULL,
  init_dpars = list(),
  type = "continuous",
  post = list(),
  sim = NULL,
  primary_dpars = "mu",
  lcdf = NULL,
  extra_pars = NULL,
  drop_intercept = FALSE
)

custom_family(
  family,
  dpars,
  links,
  lpdf,
  valid_y = NULL,
  init_dpars = list(),
  type = "continuous",
  post = list(),
  sim = NULL,
  primary_dpars = "mu",
  lcdf = NULL,
  extra_pars = NULL,
  drop_intercept = FALSE
)
```

## Arguments

- family:

  Character name of the family.

- dpars:

  Character vector of distributional parameter names. The first entry
  must be `"mu"`.

- links:

  Named list mapping each dpar to a link name (see
  `frmtmb:::frmtmb_links`) or a link object.

- lpdf:

  Function `(y, dpars, aterms)` returning the vectorized log-density.
  `dpars` is a named list of advector vectors; `aterms` is a named list
  of numeric addition-term values (for example `trials`).

- valid_y:

  Optional function `(y, aterms)` that signals an error for invalid
  responses. Called once at assembly time.

- init_dpars:

  Optional named list of functions `(y, aterms)` giving a response-scale
  starting value per dpar (applied to the intercept through the link).

- type:

  One of `"continuous"`, `"discrete"`, `"ordinal"`.

- post:

  Named list of numeric helper functions `mean_fn(dpars, aterms)` and
  `var_fn(dpars, aterms)`, used by
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html).

- sim:

  Optional numeric simulator `(dpars, aterms, n)` returning `n` response
  draws; used by [`simulate()`](https://rdrr.io/r/stats/simulate.html).

- primary_dpars:

  Which dpars receive the main model formula (default `"mu"`). Families
  with several location predictors (for example multinomial's
  per-category `mu2`, `mu3`, ...) list them all; these live in the
  `beta` parameter vector and are integrated out under REML.

- lcdf:

  Optional vectorized AD log-safe CDF `(q, dpars, aterms)` returning
  probabilities; enables `cens()` and
  [`trunc()`](https://rdrr.io/r/base/Round.html) addition terms.

- extra_pars:

  Optional function `(y, aterms)` returning a named list of numeric
  starting vectors for family-level parameters outside the dpar system
  (for example ordinal thresholds). They join the parameter template
  under their own names and reach `lpdf` as its fourth argument.

- drop_intercept:

  If `TRUE`, the intercept column is removed from the main formula's
  design matrix (ordinal families: thresholds take its place).

## Value

An object of class `frmtmb_family`.
