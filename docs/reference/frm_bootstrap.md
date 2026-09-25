# Parametric bootstrap

Simulates `nsim` response vectors from the fitted model (by default with
new group effects and new smooth coefficients each draw; see
`re_formula`), refits the model to each through
[`refit()`](https://aforren1.github.io/frmtmb/reference/refit.md)
(warm-started, no re-parsing), and collects `FUN` of every refit. Draws
whose refit fails are kept as `NA` rows; draws whose optimizer does not
report convergence are kept but flagged.

## Usage

``` r
frm_bootstrap(
  fit,
  FUN = function(f) fixef(f, flatten = TRUE),
  nsim = 500,
  seed = NULL,
  re_formula = NA
)
```

## Arguments

- fit:

  A `frmtmb_fit` for a univariate model.

- FUN:

  Function of a `frmtmb_fit` returning a numeric vector. Default:
  `fixef(f, flatten = TRUE)`, the fixed effects named as
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) names its rows, so a
  bootstrap standard error lines up with the Wald one by name.

- nsim:

  Number of bootstrap draws.

- seed:

  Optional seed.

- re_formula:

  Which random terms stay at their fitted values; the rest are redrawn
  in every replicate. The default `NA` (or `~0`, `~1`) redraws every
  group-level effect AND the penalized coefficients of every smooth,
  `gp()` and `hsgp()` term from their fitted laws, a whole-model
  parametric bootstrap in the manner of lme4's `bootMer(use.u = FALSE)`
  with the smooths treated as random effects. `NULL` conditions on the
  fitted random effects and smooths and redraws the observation noise
  alone. A one-sided formula keeps the group-level terms it names at
  their fitted values and the smooths with them, and redraws the other
  group-level terms, as
  [`simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.md)
  reads it.

  This is not `simulate(re_formula = NA)`, which holds a population
  smooth at its fitted curve, because a posterior-predictive check needs
  the curve the model estimated.

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

## Examples

``` r
set.seed(3)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
bs <- frm_bootstrap(fit, nsim = 20, seed = 1)
bs
#> Parametric bootstrap: 20 refits, 0 failed or not converged
#> 
#>                   estimate      bias       se       lwr     upr
#> (Intercept)       1.235000 -0.017843 0.137530  0.943890 1.40490
#> x                 0.570230  0.028611 0.126600  0.423170 0.88196
#> sigma_(Intercept) 0.055309  0.017652 0.086742 -0.071832 0.24671
confint(bs)
#>                           lwr       upr       est
#> (Intercept)        0.94389441 1.4049435 1.2349782
#> x                  0.42317128 0.8819600 0.5702318
#> sigma_(Intercept) -0.07183161 0.2467123 0.0553095
```
