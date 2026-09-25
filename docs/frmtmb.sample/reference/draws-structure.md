# Model structure behind a set of draws

[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`formula()`](https://rdrr.io/r/stats/formula.html),
[`family()`](https://rdrr.io/r/stats/family.html),
[`getCall()`](https://rdrr.io/r/stats/update.html) and
[`ngrps()`](https://aforren1.github.io/frmtmb/reference/ngrps.html) on a
`frmtmb_draws` report the model the sampler ran, by delegating to the
fit stored inside it. They read structure only, so they work on draws
from the formula route, which has no maximum-likelihood estimate.

## Usage

``` r
# S3 method for class 'frmtmb_draws'
fixef(
  object,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  pars = NULL,
  ...
)

# S3 method for class 'frmtmb_draws'
VarCorr(
  x,
  sigma = 1,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  ...
)

# S3 method for class 'frmtmb_draws'
ranef(
  object,
  summary = TRUE,
  robust = FALSE,
  probs = c(0.025, 0.975),
  pars = NULL,
  groups = NULL,
  ...
)

# S3 method for class 'frmtmb_draws'
nobs(object, ...)

# S3 method for class 'frmtmb_draws'
formula(x, ...)

# S3 method for class 'frmtmb_draws'
family(object, ...)

# S3 method for class 'frmtmb_draws'
getCall(x, ...)

# S3 method for class 'frmtmb_draws'
coef(object, summary = TRUE, robust = FALSE, probs = c(0.025, 0.975), ...)
```

## Arguments

- object, x:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md).

- summary:

  If `TRUE` (brms's default), summaries; otherwise the draws, as above.

- robust:

  If `TRUE`, median and MAD instead of mean and SD.

- probs:

  The quantiles to report.

- pars:

  For [`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) and
  [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html), the
  coefficients to keep, by name without the `b_` prefix, as in brms.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

- sigma:

  Ignored, as in brms.

- groups:

  For [`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html), the
  grouping factors to keep.

## Value

As for the corresponding `frmtmb_fit` method.

## Details

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html),
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) and
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) are posterior
quantities, not structural ones, and they are brms's methods on these
draws, compared [`identical()`](https://rdrr.io/r/base/identical.html)
against brms's own installed methods in `dev/brmsnames-findings.md`.
[`fixef()`](https://rdrr.io/pkg/nlme/man/fixed.effects.html) is a
coefficients x statistics matrix;
[`ranef()`](https://rdrr.io/pkg/nlme/man/random.effects.html) and
[`coef()`](https://rdrr.io/r/stats/coef.html) are a list keyed by
grouping factor of `levels x statistics x coefficients` arrays, and
[`coef()`](https://rdrr.io/r/stats/coef.html) broadcasts every
population-level coefficient over the levels, as brms's does;
[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html) is a list keyed
by grouping factor with `sd`, and `cor` and `cov` when the group has
correlations, then `residual__`. With `summary = FALSE` each returns
brms's raw draws instead: a draws x coefficients matrix, a
`draws x levels x coefficients` array, or a
`draws x coefficients x coefficients` array.

[`VarCorr()`](https://rdrr.io/pkg/nlme/man/VarCorr.html)'s standard
deviations and correlations are computed per draw from the sampled
covariance parameters, because the sampler stores `theta` and not brms's
`sd_` and `cor_` draws.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)
  ds <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                   data = dd, chains = 1, iter = 500, refresh = 0)
  nobs(ds)
  ngrps(ds)
  coef(ds)$g[1:3, , "Intercept"]
  dim(ranef(ds, summary = FALSE)$g)
  VarCorr(ds)$g$sd
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   sigma              student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: The largest R-hat is 1.08, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#>            Estimate Est.Error       Q2.5   Q97.5
#> Intercept 0.3791911 0.2658016 0.01011916 1.05593
# }
```
