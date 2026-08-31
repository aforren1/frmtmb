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
#> mu.(Intercept)    1.235000 -0.017843 0.137530  0.943890 1.40490
#> mu.x              0.570230  0.028611 0.126600  0.423170 0.88196
#> sigma.(Intercept) 0.055309  0.017652 0.086742 -0.071832 0.24671
confint(bs)
#>                           lwr       upr       est
#> mu.(Intercept)     0.94389441 1.4049435 1.2349782
#> mu.x               0.42317128 0.8819600 0.5702318
#> sigma.(Intercept) -0.07183161 0.2467123 0.0553095
```
