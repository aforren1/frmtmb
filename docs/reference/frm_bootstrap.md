# Parametric bootstrap

Simulates `nsim` response vectors from the fitted model (by default with
new random effects each draw), refits the model to each through
[`refit()`](refit.md) (warm-started, no re-parsing), and collects `FUN`
of every refit. Draws whose refit fails are kept as `NA` rows; draws
whose optimizer does not report convergence are kept but flagged.

## Usage

``` r
frm_bootstrap(
  fit,
  FUN = function(f) unlist(fixef(f)),
  nsim = 500,
  seed = NULL,
  re.form = NA
)
```

## Arguments

- fit:

  A `frmtmb_fit` for a univariate model.

- FUN:

  Function of a `frmtmb_fit` returning a numeric vector. Default: the
  flattened fixed effects.

- nsim:

  Number of bootstrap draws.

- seed:

  Optional seed.

- re.form:

  Passed to [`simulate()`](https://rdrr.io/r/stats/simulate.html); the
  default `NA` simulates marginally (new random effects), which is the
  standard parametric bootstrap for mixed models.

## Value

A `frmtmb_boot` object: `t0` (FUN at the original fit), `t` (`nsim` x
`length(t0)` matrix), and `converged`.
[`confint()`](https://rdrr.io/r/stats/confint.html) gives percentile
intervals.

## Details

There is no standard `bootstrap` generic to implement
([`boot::boot`](https://rdrr.io/pkg/boot/man/boot.html) and
[`lme4::bootMer`](https://rdrr.io/pkg/lme4/man/bootMer.html) are plain
functions), hence the `frm_` prefix.
