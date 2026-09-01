# Sample the fitted model with NUTS

Runs
[`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html) on
the fitted objective, initialized at the ML estimates (which shortens
warmup considerably), and returns the draws with frmtmb coefficient
names. Without priors this samples the likelihood with flat improper
priors on the outer parameters - the random effects get their proper
hierarchical Gaussian terms - so treat the result as an ML diagnostic
(see
[`check_laplace()`](https://aforren1.github.io/frmtmb/reference/check_laplace.md))
rather than a full Bayesian analysis; posteriors can be improper for
variance components with few groups.

## Usage

``` r
frm_sample(
  fit,
  ...,
  priors = NULL,
  lower = NULL,
  upper = NULL,
  init = "last.par.best",
  init_jitter = 0.25
)
```

## Arguments

- fit:

  A `frmtmb_fit`.

- ...:

  Passed to
  [`tmbstan::tmbstan()`](https://rdrr.io/pkg/tmbstan/man/tmbstan.html)
  (`chains`, `iter`, `laplace`, `cores`, ...). On Windows `cores > 1`
  falls back to sequential chains with a warning: parallel chains run on
  socket workers, which can evaluate neither the RTMB tape nor the
  objective closure (the known RTMB limitation of tmbstan, tmbstan#27).
  Fork clusters on unix can, so `cores` works there.

- priors:

  Optional named list of priors (see
  [`prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md));
  names are parameter names as in the draws (or whole components:
  `"beta"`, `"theta"`, ...). Parameters without a prior keep the flat
  improper default. The objective is re-taped with the prior terms
  added; the ML fit itself is unchanged.

- lower, upper:

  Optional named numeric vectors of hard bounds on outer parameters
  (brms `lb`/`ub`), applied on the internal scale through Stan's
  constrained transforms. Chain starting values are clamped strictly
  inside the bounds; a bound that excludes the ML mode itself warns,
  because the chains then no longer start there.

- init:

  Initialization; the default starts chain 1 exactly at the ML mode and
  every further chain at the mode plus a normal perturbation of sd
  `init_jitter` on the unconstrained scale. The mode anchor keeps warmup
  short; the jitter keeps the chains overdispersed enough for Rhat to
  retain power against multimodality (the standard objection to
  identical mode starts). `"random"` requests Stan's own overdispersed
  initialization.

- init_jitter:

  Per-chain perturbation sd for the default init; `0` starts every chain
  exactly at the mode. Draws from the R session's RNG, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) makes the inits
  reproducible.

## Value

An object of class `frmtmb_draws`: list with the `stanfit`, a draws
matrix with named columns
([`as.matrix()`](https://rdrr.io/r/base/matrix.html) method), and the
originating fit.

## Multimodal posteriors

For
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
fits the posterior is multimodal by construction (label switching at
minimum). Mode-centered inits, jittered or not, leave every chain in one
symmetry branch, so Rhat cannot flag the others; use `init = "random"`
there and inspect chains individually.

## Examples

``` r
# \donttest{
set.seed(9)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
ds <- frm_sample(fit, chains = 1, iter = 500, refresh = 0)
#> Warning: There were 3 divergent transitions after warmup. See
#> https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
#> to find out why this is a problem and how to eliminate them.
#> Warning: Examine the pairs() plot to diagnose sampling problems
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
summary(ds)
#>                         mean        sd       2.5%      97.5%     n_eff
#> (Intercept)        0.8452223 0.1318395  0.5597141  1.0977453 132.09215
#> x                  0.5997377 0.1124987  0.3760346  0.7664300  90.55402
#> sigma_(Intercept) -0.0330158 0.0800034 -0.1722470  0.1183057 122.47510
#> theta_1           -2.3309557 0.8768774 -3.8887436 -0.8218252  16.81438
#>                        Rhat
#> (Intercept)       1.0168741
#> x                 0.9990810
#> sigma_(Intercept) 0.9960028
#> theta_1           0.9980670
fixef(ds)
#>                     Estimate Est.Error       Q2.5     Q97.5
#> (Intercept)        0.8452223 0.1318395  0.5597141 1.0977453
#> x                  0.5997377 0.1124987  0.3760346 0.7664300
#> sigma_(Intercept) -0.0330158 0.0800034 -0.1722470 0.1183057
hypothesis(ds, "sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)")
#> Hypothesis tests (method = posterior)
#>                                         hypothesis estimate      se       lwr
#>  sd_g__Intercept^2 / (sd_g__Intercept^2 + sigma^2)  0.03404 0.05642 0.0004609
#>     upr      z        p
#>  0.2145 0.6032 0.007968
# }
```
