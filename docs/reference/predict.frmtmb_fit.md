# Predictions from a frmtmb fit

brms's [`predict()`](https://rdrr.io/r/stats/predict.html): a summary of
the PREDICTIVE distribution of the response, with observation noise
included, in the columns `Estimate`, `Est.Error`, `Q2.5` and `Q97.5`. It
is a much wider interval than
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html)'s, which
summarizes the expected response.

A maximum-likelihood fit has no posterior, so the draws are simulated.
Each replicate draws the outer parameter vector from its asymptotic
normal law, `N(theta_hat, vcov(object, full = TRUE))`, evaluates every
distributional parameter there, and draws a response from the family's
own simulator, the one
[`simulate()`](https://rdrr.io/r/stats/simulate.html) uses. The interval
is the empirical quantile of those draws, so it carries the family's
skew and its discreteness: a poisson predictive interval is on the
counts.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
predict(
  object,
  newdata = NULL,
  re_formula = NULL,
  transform = NULL,
  resp = NULL,
  negative_rt = FALSE,
  ndraws = NULL,
  draw_ids = NULL,
  sort = FALSE,
  ntrys = 100L,
  cores = NULL,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...,
  allow_new_levels = FALSE,
  sample_new_levels = NULL,
  param_uncertainty = TRUE
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- newdata:

  Optional data frame to predict on. Defaults to the training data.

- re_formula:

  `NULL` (default) keeps the random effects, so the draw is conditional
  on the modes; `NA` or `~0` draws at the population level.

- transform:

  A function applied to the draws before they are summarized, as in
  brms.

- resp:

  For multivariate fits: which response, or `NULL` (default) for all of
  them, which is brms's `nrow x 4 x nresp` array with the responses
  named on the third dimension.

- negative_rt:

  brms's sign convention for its `wiener` family. Refused, with the
  reason: this package's evidence-accumulation families return the
  response and the time as they declare them.

- ndraws:

  Number of simulated replicates. Defaults to 1000. brms thins its
  posterior with this; here it sets how many draws are taken, so a
  larger number narrows the Monte Carlo error of the quantile columns
  and costs proportionally more.

- draw_ids, sort, cores:

  Refused by name. There are no stored draws to index, rows always come
  back in the order of the data, and the simulation is not parallelized.

- ntrys:

  Rejection-sampling attempts per row for a truncated response, brms's
  spelling of the family simulator's `max_iter`.

- summary:

  If `FALSE`, the `ndraws x nrow` matrix of simulated draws instead of
  the summary.

- robust:

  If `TRUE`, the median and the median absolute deviation of the draws
  instead of the mean and the standard deviation. brms's argument, and
  it is answerable here because these draws exist.

- probs:

  Probabilities of the quantile columns.

- ...:

  Refused. An argument this method does not have is an error naming it,
  and `type`, `se.fit`, `dpar` and `scale` are refused with the name of
  the function that took over each one.

- allow_new_levels:

  Predict rows whose grouping-factor level the fit never saw, instead of
  erroring. Each replicate draws that level's effect from the block's
  estimated covariance.

- sample_new_levels:

  brms's argument. `"gaussian"`, which is what happens, or `NULL`.
  `"uncertainty"` and `"old_levels"` resample the posterior draws of the
  levels that were seen, and a maximum-likelihood fit has none, so both
  are refused by name.

- param_uncertainty:

  If `FALSE`, simulate at the estimates alone (the plug-in predictive
  distribution) instead of drawing the parameters first.

## Value

With `summary = TRUE` (the default) an `nrow x 4` matrix with the
columns `Estimate`, `Est.Error` and one per entry of `probs`; for an
ordinal or categorical response an `nrow x K` matrix of simulated
category proportions; for a multivariate fit an `nrow x 4 x nresp` array
with the responses named. With `summary = FALSE` the `ndraws x nrow`
matrix of draws, or `ndraws x nrow x nresp` for a multivariate fit.

## Which uncertainty is in the interval

Two of the three sources.

- Observation noise, from the family's simulator. This is almost always
  the largest term.

- Uncertainty in the estimates, from the asymptotic covariance of the
  outer parameter vector. `param_uncertainty = FALSE` drops it and
  simulates at the estimates alone.

- Uncertainty in the random-effect modes of a level the fit SAW is NOT
  included. The draw is conditional on them, the convention every other
  method here follows. On a design with few observations per grouping
  level that term is not small, and the interval is then narrower than a
  fully marginal one. `dev/shapes-findings.md` reports the measured
  coverage for a design with and without a grouping factor.

A level the fit did NOT see is different: there is no mode to condition
on, so `allow_new_levels = TRUE` draws its effect from the block's own
estimated covariance, once per replicate, at that replicate's
parameters. This is brms's `sample_new_levels = "gaussian"` and it is
the whole between-group variance, so the interval at an unseen level is
wider than at a known one.

A fit whose covariance could not be recovered from the Hessian has no
second term to draw: the simulation falls back to the estimates and
warns, rather than drawing from a matrix of `NaN`.

Under `REML = TRUE` (or `frmtmb_control(profile = TRUE)`) the fixed
effects are integrated out of the outer problem, so they are not in
`vcov(object, full = TRUE)` and are not drawn either. The interval there
carries the covariance parameters' uncertainty and the observation
noise, and holds the coefficients fixed.

## Ordinal and categorical responses

There is no mean to summarize, so
[`predict()`](https://rdrr.io/r/stats/predict.html) returns the
simulated proportion of each category, one column per category named
`P(Y = k)`, which is brms's shape.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
MODELLED category probabilities with their standard errors instead,
which is the smoother quantity.

## See also

[`fitted.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/fitted.frmtmb_fit.md)
for the expected response,
[`frm_linpred()`](https://aforren1.github.io/frmtmb/reference/frm_linpred.md)
for the linear predictor and the glmmTMB scale vocabulary,
[`simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.md)
for the draws themselves, and
[frmtmb-scales](https://aforren1.github.io/frmtmb/reference/frmtmb-scales.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(100))
dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x))
fit <- frm(bf(y ~ x) + poisson(), data = dd)

# the predictive summary: observation noise is in it
head(predict(fit, ndraws = 200))
#>      Estimate Est.Error Q2.5 Q97.5
#> [1,]    0.765 0.8738381    0 3.000
#> [2,]    1.270 1.0876650    0 4.000
#> [3,]    0.765 0.8445956    0 3.000
#> [4,]    2.430 1.6087387    0 6.000
#> [5,]    1.235 1.1072057    0 4.000
#> [6,]    0.665 0.8099724    0 2.025
# against the expected response, which is much tighter
head(fitted(fit))
#>       Estimate Est.Error      Q2.5     Q97.5
#> [1,] 0.7751024 0.1130250 0.5535775 0.9966273
#> [2,] 1.1849807 0.1126043 0.9642804 1.4056811
#> [3,] 0.6946365 0.1132732 0.4726251 0.9166479
#> [4,] 2.4828747 0.3477984 1.8012024 3.1645469
#> [5,] 1.2791025 0.1161520 1.0514488 1.5067563
#> [6,] 0.7001766 0.1132798 0.4781522 0.9222009

# the draws themselves
dim(predict(fit, summary = FALSE, ndraws = 50))
#> [1]  50 100

# the linear predictor is frm_linpred() now
head(frm_linpred(fit, type = "link"))
#>          1          2          3          4          5          6 
#> -0.2547601  0.1697265 -0.3643666  0.9094170  0.2461587 -0.3564227 
```
