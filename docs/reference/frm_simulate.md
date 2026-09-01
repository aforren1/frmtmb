# Simulate responses from a formula and parameters

Builds the design from `formula` and `data` exactly as
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) would,
sets the parameters from `newparams` (and/or draws them from `priors`),
and simulates responses. Random effects are redrawn from their
covariance for every simulation unless `newparams$b` supplies them.

## Usage

``` r
frm_simulate(
  formula,
  data,
  family = NULL,
  newparams = NULL,
  priors = NULL,
  nsim = 1,
  seed = NULL,
  data2 = list()
)
```

## Arguments

- formula:

  A [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formula
  (with a family attached) or a plain formula plus `family`.

- data:

  Model data, including a dummy response column.

- family:

  Family, when `formula` does not carry one.

- newparams:

  Named list of parameters, in either spelling (see Details). Optional
  when `priors` pins everything.

- priors:

  A
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  specification to draw parameters from, once per simulation.

- nsim, seed:

  As in [`simulate()`](https://rdrr.io/r/stats/simulate.html).

- data2:

  Structural objects, as in
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md).

## Value

A data frame with `nsim` columns of simulated responses, carrying the
drawn parameters in `attr(., "pars")` when `priors` is used.

## Details

`data` must contain a response column with values that are valid for the
family (any dummy values do; they only anchor the design).

## Two spellings for `newparams`

*Natural scale* (recommended): the names
[`variables()`](https://aforren1.github.io/frmtmb/reference/variables.md)
and
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
use, one number each.

- fixed coefficients under their
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) names with parentheses
  stripped (`Intercept`, `x`, `sigma_Intercept`, ...), on the LINK
  scale;

- `sigma` (and any other intercept-only dispersion dpar: `shape`, `phi`,
  `zi`, ...) on the RESPONSE scale, so `sigma = 0.7` is a residual SD of
  0.7;

- `sd_<group>__<term>` for random-effect standard deviations,
  `cor_<group>__<t1>__<t2>` for their correlations, and `sds_<label>`
  for a smooth's smoothing SD. Unset correlations are 0.

*Internal scale*: named after the `par_template` components - `beta`,
`betad`, `theta`, and optionally `b` - each a full-length vector, on the
internal parameterization (`theta` holds log SDs and Cholesky
correlation parameters, `betad` holds dispersion dpars on their link
scale). Inspect the layout with
`frm(formula, data, dry_run = "frame")$par_template`.

The two spellings cannot be mixed: `newparams` is read as internal when
every name is a `par_template` component (or `b`), and as natural
otherwise. The natural spelling covers `us`, `diag`, `homdiag`,
`smooth`, `gr_cov` and `gr_prec` random-effect structures; other
structures (`ar1`, `cs`, `toep`, GPs, ...) have no inverse map from SDs
and correlations and need the internal spelling.

## Prior-predictive simulation

With `priors`, each of the `nsim` simulations draws a fresh parameter
vector from the
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
specification and simulates from it - the analog of brms's
`sample_prior = "only"` followed by
[`posterior_predict()`](https://aforren1.github.io/frmtmb/reference/posterior_epred.md).
Draws are taken on the scale each class defines: class
`"b"`/`"Intercept"` on the link scale, class `"sd"` on the natural SD
scale (mapped to `theta` afterwards), class `"theta"` on the internal
scale. `lb`/`ub` truncate by rejection. The drawn values come back as
`attr(result, "pars")`, one row per simulation, so a prior-predictive
check can relate parameters to outcomes.

Parameters without a prior keep their `newparams` value. Whenever
`priors` are used, or `newparams` uses the natural spelling, every fixed
coefficient and every random-effect SD must be pinned by one or the
other; an unpinned parameter is an error rather than a silent zero
effect or unit SD.

## Examples

``` r
# power analysis: simulate from a design with chosen parameters
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)), y = 0)
sims <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                     newparams = list(Intercept = 1, x = 0.5,
                                      sigma = 0.7,
                                      sd_g__Intercept = 0.5),
                     nsim = 3, seed = 1)
head(sims)
#>      sim_1      sim_2      sim_3
#> 1 1.007314 0.50436944  0.7886221
#> 2 1.731551 1.20181120  0.3238773
#> 3 0.817993 0.03184368 -0.3540123
#> 4 1.327337 2.03378663  0.7988518
#> 5 1.936121 0.64057031 -0.6028489
#> 6 1.028016 0.81116086  0.5767491
# the same thing on the internal scale
frm(bf(y ~ x + (1 | g)) + gaussian(), dd,
    dry_run = "frame")$par_template   # the required layout
#> $beta
#> (Intercept)           x 
#>           0           0 
#> 
#> $betad
#> sigma_(Intercept) 
#>                 0 
#> 
#> $b
#> [1] 0 0 0 0 0 0
#> 
#> $theta
#> [1] 0
#> 
sims2 <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                      newparams = list(beta = c(1, 0.5),
                                       betad = log(0.7),
                                       theta = log(0.5)),
                      nsim = 3, seed = 1)
# prior-predictive draws
pp <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                   priors = set_prior("normal(0, 1)", class = "b") +
                     set_prior("normal(0, 2)", class = "Intercept") +
                     set_prior("exponential(1)", class = "sd") +
                     set_prior("exponential(1)", class = "Intercept",
                               dpar = "sigma"),
                   nsim = 4, seed = 1)
head(attr(pp, "pars"))
#>            x  Intercept sd_g__Intercept sigma_Intercept
#> 1 -0.6264538  0.3672866       1.9997498       0.4580301
#> 2  0.4755095 -1.4198929       0.4586192       0.8145358
#> 3 -1.9143594  2.3531666       3.3072809       1.4953203
#> 4  0.5101084 -0.3287517       0.3260219       1.5316271
```
