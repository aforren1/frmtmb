# Within-group residual correlation (R-side autocorrelation)

[`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`, `cosy()`
and `unstr()` are written as terms of the model formula, next to the
fixed and random effects, and make the residuals of one group a single
correlated draw instead of independent ones:

## Value

[`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`, `cosy()`
and `unstr()` are formula terms, not free-standing functions:
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) reads them
at parse time, and the value they contribute is the fitted
autocorrelation block of the model, reachable through
[`autocor_matrix()`](https://aforren1.github.io/frmtmb/reference/autocor_matrix.md),
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md)
and
[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md).
This page itself documents the term grammar and returns nothing.

## Details

\$\$y_g \sim N(\mu_g,\\ D_g R D_g), \qquad D_g = \mathrm{diag}(\sigma_i,
i \in g).\$\$

`R` is a unit-diagonal correlation matrix over the time points, so
`sigma` keeps its usual meaning - the marginal residual standard
deviation - and a `sigma ~ ...` distributional model enters through the
diagonal. Nothing is added to the linear predictor, so
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`predict.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/predict.frmtmb_fit.md)
and `se.fit` are unchanged; what changes is the likelihood. This is the
model `nlme::gls(correlation = corAR1())` fits, and the one brms fits
under `cov = TRUE`.

## Structures

- `ar(time, gr, p = 1, cov = TRUE)`:

  Stationary AR(`p`) correlation. `p = 1` gives \\R\_{ij} =
  \rho^{\|i-j\|}\\.

- `ma(time, gr, q = 1, cov = TRUE)`:

  Invertible MA(`q`): correlation dies after lag `q`.

- `arma(time, gr, p = 1, q = 1, cov = TRUE)`:

  Stationary, invertible ARMA(`p`, `q`).

- `cosy(time, gr)`:

  Compound symmetry - one correlation shared by every pair. Equivalent
  in fit to a random intercept per group, but the correlation may also
  be negative (down to `-1 / (d - 1)`), which a variance component
  cannot express.

- `unstr(time, gr)`:

  One free correlation per pair of time levels, through the same
  Cholesky parameterization the `us()` random-effect structure uses.

The AR and MA coefficients are estimated through partial
autocorrelations (the Monahan/Jones transform
[`nlme::corARMA`](https://rdrr.io/pkg/nlme/man/corARMA.html) also uses),
so every parameter value is a stationary and invertible process and the
optimizer cannot leave the parameter space. The ARMA autocorrelation
function is exact, not a truncated MA(\\\infty\\) expansion.

## Arguments

The argument order is brms's, so the FIRST positional argument is `time`
and the second is `gr`: write `cosy(gr = subj)`, not `cosy(subj)`.

- `time`:

  The variable whose levels index the correlation. Omit it (brms's
  `time = NA`) to use each row's position within its group.

- `gr`:

  The grouping variable the residual factorizes over. Omit it to treat
  the whole data set as one series.

- `p`, `q`:

  Autoregressive and moving-average orders.

- `cov`:

  Must be `TRUE`. brms's default `cov = FALSE` is a different likelihood
  (a residual regression that conditions on each group's first rows),
  which is not implemented; the call is refused rather than silently
  reinterpreted.

## Families

[`gaussian()`](https://rdrr.io/r/stats/family.html) and
[`student()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
only - the two families with a real residual, and exactly the two brms
treats this way
([`student()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.md)
gets the multivariate-t analog, with one `nu` per group, so a predicted
`nu ~ ...` is refused). brms accepts the same spelling for other
families but fits a different model there: a latent Gaussian AR process
added to the linear predictor. That model has a spelling of its own
here - a random effect over the time factor,
`+ ar1(factor(week) + 0 | subj)`, or `toep()` / `us()` for a freer lag
structure - and the refusal names it.

## Time points, gaps and ragged groups

The lag between two rows is the distance between their positions in the
GLOBAL set of time levels, so a group missing week 3 gets
`cor(week2, week4) = rho^2`. That is `nlme`'s reading
(`corAR1(form = ~ week | subj)`) and what the agreement tests pin down;
brms instead indexes
[`ar()`](https://rdrr.io/r/stats/ar.html)/`ma()`/`arma()`/`cosy()` by
position within the group and treats a missing row as no gap (only its
`unstr()` carries a time index). On complete balanced groups the two
coincide. Time levels that are whole numbers but not consecutive warn,
because the lag is then not the one the labels suggest.

Groups of different sizes are handled by construction: the density is
evaluated one PATTERN at a time (a pattern being a distinct set of
present time levels), so a balanced design costs one Cholesky per
gradient evaluation and ragged data costs one per distinct pattern.
Duplicate `(gr, time)` pairs are refused, as they are in brms.

## What cannot be combined with it

The likelihood is a joint density over each group, so it no longer
factorizes into per-row contributions.
[`weights()`](https://rdrr.io/r/stats/weights.html), `cens()`,
[`trunc()`](https://rdrr.io/r/base/Round.html), `se()` and `mi()` on the
same response are therefore refused, as are `rescor = TRUE`, mixture
families, `quadrature = TRUE`,
[`frm_simulate()`](https://aforren1.github.io/frmtmb/reference/frm_simulate.md)
and `residuals(type = "osa")`. brms refuses the same core set. Random
effects ARE allowed and are the point of the feature: the marginal
likelihood is a Laplace approximation over the modes with the correlated
residual density inside, which reproduces
`nlme::lme(random = ~ 1 | subj, correlation = corAR1())`.

## Where the parameters appear

On the internal (unconstrained) scale they are the `thetaac_*` rows of
[`confint.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/confint.frmtmb_fit.md)
and of `vcov(full = TRUE)`. On the natural scale,
[`summary()`](https://rdrr.io/r/base/summary.html) prints them under
"Within-group residual correlation" and
[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
reports one row per parameter under brms's names - `ar[1]`, `ma[1]`,
`cosy`, `cortime__<t1>__<t2>` - with a delta-method interval (Fisher-z
for the bounded ones, the identity for the coefficients of a
higher-order AR/MA process).
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
sees them as `ar1`, `ma1`, `cosy` and `cortime__<t1>__<t2>`.
[`autocor_matrix()`](https://aforren1.github.io/frmtmb/reference/autocor_matrix.md)
returns the fitted `R`. They are NOT part of
[`VarCorr()`](https://aforren1.github.io/frmtmb/reference/VarCorr.md),
which reports random-effect blocks.
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
cannot target them yet.

## Divergence from brms

brms parameterizes [`ar()`](https://rdrr.io/r/stats/ar.html), `ma()` and
`arma()` by the INNOVATION standard deviation (its `cholesky_cor_ar1()`
divides by `1 - ar^2`), while `cosy()` and `unstr()` use the marginal
one. Here every structure uses the marginal `sigma`, so that
[`sigma()`](https://rdrr.io/r/stats/sigma.html), pearson residuals and a
`sigma ~ x` model mean one thing throughout. The correlation parameters
agree with brms exactly; the scales relate by
`sigma_marginal = sigma_innovation / sqrt(1 - phi^2)` for AR(1) and
`sigma_marginal = sigma_innovation * sqrt(1 + theta^2)` for MA(1). brms
also limits `cov = TRUE` to order one ("Covariance formulation of ARMA
structures is only possible for effects of maximal order one"); higher
`p` and `q` are supported here.

## See also

[`autocor_matrix()`](https://aforren1.github.io/frmtmb/reference/autocor_matrix.md)
for the fitted correlation matrix,
[`confint_varcorr()`](https://aforren1.github.io/frmtmb/reference/confint_varcorr.md)
for natural-scale intervals, and
[`vignette("brms-migration")`](https://aforren1.github.io/frmtmb/articles/brms-migration.md)
for the porting notes.

## Examples

``` r
set.seed(1)
d <- expand.grid(week = 1:5, subj = factor(1:30))
d$x <- rnorm(150)
e <- as.vector(vapply(1:30, function(i) {
  as.vector(stats::filter(rnorm(5), 0.6, "recursive"))
}, numeric(5)))
d$y <- 1 + 0.5 * d$x + e

fit <- frm(bf(y ~ x + ar(week, subj, cov = TRUE)) + gaussian(),
           data = d)
summary(fit)
#> Family: gaussian 
#> Formula: y ~ x + ar(week, subj, cov = TRUE) 
#> Method: ML   nobs: 150 
#> logLik: -214.257  AIC: 436.514  BIC: 448.556 
#> 
#> Within-group residual correlation: ar(week, subj, cov = TRUE)
#>       Estimate  2.5 % 97.5 %
#> ar[1]   0.5226 0.3675 0.6494
#> 
#> Coefficients (mu):
#>             Estimate Std. Error z value  Pr(>|z|)
#> (Intercept)  1.10730    0.13956  7.9343 2.116e-15
#> x            0.45411    0.08061  5.6334 1.767e-08
#> 
#> Coefficients (sigma):
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept) 0.137029   0.071079  1.9278  0.05387
autocor_matrix(fit)
#>            1         2         3         4          5
#> 1 1.00000000 0.5225995 0.2731103 0.1427273 0.07458923
#> 2 0.52259955 1.0000000 0.5225995 0.2731103 0.14272731
#> 3 0.27311029 0.5225995 1.0000000 0.5225995 0.27311029
#> 4 0.14272731 0.2731103 0.5225995 1.0000000 0.52259955
#> 5 0.07458923 0.1427273 0.2731103 0.5225995 1.00000000

# compound symmetry, and the unstructured correlation over the five
# weeks
frm(bf(y ~ x + cosy(week, subj)) + gaussian(), data = d)
#> frmtmb fit: y ~ x + cosy(week, subj) 
#> Family: gaussian   Method: ML 
#> logLik: -220.426  AIC: 448.852  nobs: 150 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      1.0901      0.4085 
#>  sigma:
#> (Intercept) 
#>      0.1398 
frm(bf(y ~ x + unstr(week, subj)) + gaussian(), data = d)
#> frmtmb fit: y ~ x + unstr(week, subj) 
#> Family: gaussian   Method: ML 
#> logLik: -213.317  AIC: 452.634  nobs: 150 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      1.1178      0.4386 
#>  sigma:
#> (Intercept) 
#>      0.1371 

# a random intercept alongside the correlated residual is allowed
frm(bf(y ~ x + (1 | subj) + ar(week, subj, cov = TRUE)) + gaussian(),
    data = d)
#> frmtmb fit: y ~ x + (1 | subj) + ar(week, subj, cov = TRUE) 
#> Family: gaussian   Method: ML 
#> logLik: -214.232  AIC: 438.465  nobs: 150 
#> 
#> Fixed effects:
#>  mu:
#> (Intercept)           x 
#>      1.1060      0.4518 
#>  sigma:
#> (Intercept) 
#>       0.109 
#> 
#> Random effects:
#>   1 | subj 
#>         Name Std.Dev.
#>  (Intercept)  0.26669
```
