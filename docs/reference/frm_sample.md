# Sample the fitted model with NUTS

Runs
[`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html) on
the fitted objective, initialized at the ML estimates (which shortens
warmup considerably), and returns the draws with frmtmb coefficient
names. Without priors this samples the likelihood with flat improper
priors on the outer parameters - the random effects get their proper
hierarchical Gaussian terms - so treat the result as an ML diagnostic
(see [`check_laplace()`](check_laplace.md)) rather than a full Bayesian
analysis; posteriors can be improper for variance components with few
groups.

## Usage

``` r
frm_sample(
  fit,
  ...,
  priors = NULL,
  lower = NULL,
  upper = NULL,
  init = "last.par.best"
)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- ...:

  Passed to
  [`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html)
  (`chains`, `iter`, `laplace`, `cores`, ...).

- priors:

  Optional named list of priors (see
  [`prior_normal()`](frmtmb-priors.md)); names are parameter names as in
  the draws (or whole components: `"beta"`, `"theta"`, ...). Parameters
  without a prior keep the flat improper default. The objective is
  re-taped with the prior terms added; the ML fit itself is unchanged.

- lower, upper:

  Optional named numeric vectors of hard bounds on outer parameters
  (brms `lb`/`ub`), applied on the internal scale through Stan's
  constrained transforms.

- init:

  Initialization; defaults to the ML mode.

## Value

An object of class `frmtmb_draws`: list with the `stanfit`, a draws
matrix with named columns
([`as.matrix()`](https://rdrr.io/r/base/matrix.html) method), and the
originating fit.
