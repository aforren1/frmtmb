# Refusals for the refit-based and marginal-likelihood brmsfit methods

These `brmsfit` methods have frmtmb spellings that do not exist yet.
They are defined so that a ported script fails with the reason and the
alternative rather than with "could not find function", and they are
documented so the reason is findable.

## Usage

``` r
loo_moment_match(x, ...)

# S3 method for class 'frmtmb_draws'
loo_moment_match(x, ...)

loo_subsample(x, ...)

# S3 method for class 'frmtmb_draws'
loo_subsample(x, ...)

reloo(x, ...)

# S3 method for class 'frmtmb_draws'
reloo(x, ...)

kfold(x, ...)

# S3 method for class 'frmtmb_draws'
kfold(x, ...)

bridge_sampler(x, ...)

# S3 method for class 'frmtmb_draws'
bridge_sampler(x, ...)

bayes_factor(x, ...)

# S3 method for class 'frmtmb_draws'
bayes_factor(x, ...)

post_prob(x, ...)

# S3 method for class 'frmtmb_draws'
post_prob(x, ...)
```

## Arguments

- x, ...:

  Ignored; these methods always stop.

## Value

These functions never return; they signal an error.

## Details

- `loo_moment_match()`, `loo_subsample()`, `reloo()` and `kfold()` all
  need to refit the model on modified data.
  [`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
  is the resampling machinery frmtmb does have, and
  [`AIC()`](https://rdrr.io/r/stats/AIC.html) on the maximum-likelihood
  fits answers the comparison question directly.

- `bridge_sampler()`, `bayes_factor()` and `post_prob()` are
  marginal-likelihood quantities. A marginal likelihood is an integral
  against the PRIOR, so it is undefined under `prior = "flat"`, and even
  under the default priors the bridge-sampling estimator needs a
  normalized log-posterior evaluator that the RTMB tape does not expose.

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE)) {
  set.seed(1)
  dd <- data.frame(x = rnorm(40))
  dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
  fit <- frm(bf(y ~ x), family = gaussian(), data = dd)
  ds <- frm_sample(fit, chains = 1, iter = 400, refresh = 0)
  # each refusal names its reason and the replacement
  try(reloo(ds))
  try(bayes_factor(ds, ds))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 1.1, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Error : reloo() is not implemented for frmtmb draws: it re-runs the sampler once per observation with a high Pareto k, and frm_sample() has no stored program to re-run on modified data. Read loo()'s Pareto k table and treat a bad k as the diagnostic it is (usually many group-level parameters left to the data alone; see the prior section of ?loo), or compare the maximum-likelihood fits with AIC()
#> Error : bayes_factor() is not available for frmtmb draws: it is a ratio of the marginal likelihoods bridge_sampler() would have to estimate, and those are undefined under prior = "flat" and unavailable from the tape. hypothesis() gives the posterior probability of a directional claim, and loo() the predictive comparison
# }
```
