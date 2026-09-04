# Simulate responses from a formula and parameters

Builds the design from `formula` and `data` exactly as
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) would,
sets the parameters from `newparams` (and/or draws them from `prior`),
and simulates responses. Random effects are redrawn from their
covariance for every simulation unless `newparams$b` supplies them.

## Usage

``` r
frm_simulate(
  formula,
  data,
  family = NULL,
  newparams = NULL,
  prior = NULL,
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

  A data frame of model data, including a dummy response column.

- family:

  Family, when `formula` does not carry one.

- newparams:

  Named list of parameters, in either spelling (see Details).
  [`par_template()`](https://aforren1.github.io/frmtmb/reference/par_template.md)
  discovers the internal spelling for a formula and data. Optional when
  `prior` pins everything.

- prior:

  A
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
  specification to draw parameters from, once per simulation, or a
  `brmsprior` object brms built, which is translated row by row. The
  argument takes brms's spelling, `prior`; the `priors` of releases
  before 0.43 is gone rather than aliased.

- nsim, seed:

  As in [`simulate()`](https://rdrr.io/r/stats/simulate.html).

- data2:

  Structural objects, as in
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md).

## Value

A data frame with `nsim` columns of simulated responses, carrying the
drawn parameters in `attr(., "pars")` when `prior` is used.

## Details

`data` must contain a response column with values that are valid for the
family (any dummy values do; they only anchor the design).

Draws come back in the response's own type, exactly as
[`simulate()`](https://rdrr.io/r/stats/simulate.html)'s do: an ordered
factor for an ordinal family, an unordered one for a categorical family,
and a matrix column for a matrix response
([`multinomial()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md),
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md),
[`lca()`](https://aforren1.github.io/frmtmb/reference/lca.md)).

The structured families draw here through the same implementation
[`simulate()`](https://rdrr.io/r/stats/simulate.html) uses (see its
Structured draws section):
[`hmm()`](https://aforren1.github.io/frmtmb/reference/hmm.md) walks its
chain per sequence, `mixture(groups = ~g)` takes one class per group,
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
uses its class covariances, and a residual correlation term
([`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, ...) contributes one
correlated residual per group. The de novo frame carries those
structures, so nothing is lost relative to a fit;
[`ar()`](https://rdrr.io/r/stats/ar.html) and friends need their
`thetaac` entry in the internal `newparams` spelling, since a
correlation parameter has no natural-scale name here.

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

*Internal scale*: named after the parameter components - `beta`,
`betad`, `theta`, and optionally `b` - each a full-length vector, on the
internal parameterization (`theta` holds log SDs and Cholesky
correlation parameters, `betad` holds dispersion dpars on their link
scale).
[`par_template()`](https://aforren1.github.io/frmtmb/reference/par_template.md)
returns that layout for a formula and data without fitting anything,
filled with the default values and carrying the name of every entry;
edit it and pass it back as `newparams`.

The two spellings cannot be mixed: `newparams` is read as internal when
every name is a `par_template` component (or `b`), and as natural
otherwise. The natural spelling covers `us`, `diag`, `homdiag`,
`smooth`, `gr_cov` and `gr_prec` random-effect structures; other
structures (`ar1`, `cs`, `toep`, GPs, ...) have no inverse map from SDs
and correlations and need the internal spelling.

## Prior-predictive simulation

With `prior`, each of the `nsim` simulations draws a fresh parameter
vector from the
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
specification and simulates from it - the analog of brms's
`sample_prior = "only"` followed by `posterior_predict()`. Draws are
taken on the scale each class defines: class `"b"`/`"Intercept"` on the
link scale, class `"sd"` on the natural SD scale (mapped to `theta`
afterwards), class `"theta"` on the internal scale. `lb`/`ub` truncate
by rejection. The drawn values come back as `attr(result, "pars")`, one
row per simulation, so a prior-predictive check can relate parameters to
outcomes.

Parameters without a prior keep their `newparams` value. Whenever
`prior` are used, or `newparams` uses the natural spelling, every fixed
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
#>       sim_1        sim_2      sim_3
#> 1 0.8476542  0.344709864  0.6289625
#> 2 1.6862256  1.156485736  0.2785519
#> 3 0.2897796 -0.496369784 -0.8822257
#> 4 2.2704655  2.976915206  1.7419804
#> 5 1.7904259  0.494875611 -0.7485436
#> 6 0.2114016 -0.005453893 -0.2398657
# the same thing on the internal scale
par_template(bf(y ~ x + (1 | g)) + gaussian(), dd)   # the layout
#> <frmtmb parameter template> starting values
#> $beta  (2)
#> (Intercept)           x 
#>           0           0 
#> $betad  (1)
#> sigma_(Intercept) 
#>                 0 
#> $b  (6)
#> b_1 b_2 b_3 b_4 b_5 b_6 
#>   0   0   0   0   0   0 
#> $theta  (1)
#> theta_1 
#>       0 
#> Edit and pass back as frm(start = ) or frm_simulate(newparams = ). Parentheses are optional in names you supply.
sims2 <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                      newparams = list(beta = c(1, 0.5),
                                       betad = log(0.7),
                                       theta = log(0.5)),
                      nsim = 3, seed = 1)
# prior-predictive draws
pp <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                   prior = set_prior("normal(0, 1)", class = "b") +
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
