# Approximate leave-one-out cross-validation

[`loo()`](https://aforren1.github.io/frmtmb/reference/loo.html) runs
Pareto-smoothed importance-sampling LOO and
[`waic()`](https://aforren1.github.io/frmtmb/reference/loo.html) the
widely applicable information criterion, both on the
[`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md)
matrix, by handing it to
[`loo::loo.matrix()`](https://mc-stan.org/loo/reference/loo.html) and
[`loo::waic.matrix()`](https://mc-stan.org/loo/reference/waic.html)
unchanged. The returned objects are the loo package's own, so
[`print()`](https://rdrr.io/r/base/print.html) and
[`loo::pareto_k_table()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
work on them directly.
[`loo_compare()`](https://aforren1.github.io/frmtmb/reference/loo.html)
computes the criterion for each draws object it is given and ranks them;
handed criteria instead of draws, it is
[`loo::loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
itself. `psis()` returns the smoothed importance weights alone.

## Usage

``` r
# S3 method for class 'frmtmb_draws'
loo(x, ndraws = NULL, resp = NULL, ...)

# S3 method for class 'frmtmb_draws'
waic(x, ndraws = NULL, resp = NULL, ...)

# S3 method for class 'frmtmb_draws'
loo_compare(x, ..., criterion = c("loo", "waic"), model_names = NULL)

psis(log_ratios, ...)

# S3 method for class 'frmtmb_draws'
psis(log_ratios, ndraws = NULL, resp = NULL, ...)

# S3 method for class 'frmtmb_draws'
LOO(x, ...)

# S3 method for class 'frmtmb_draws'
WAIC(x, ...)
```

## Arguments

- x:

  A `frmtmb_draws` from
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md),
  or (for
  [`loo_compare()`](https://aforren1.github.io/frmtmb/reference/loo.html))
  already-computed criteria.

- ndraws, resp:

  Passed to
  [`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md).

- ...:

  Further models for
  [`loo_compare()`](https://aforren1.github.io/frmtmb/reference/loo.html);
  otherwise passed to the loo package function.

- criterion:

  Which criterion
  [`loo_compare()`](https://aforren1.github.io/frmtmb/reference/loo.html)
  computes for each draws object.

- model_names:

  Row names for the comparison; the default deparses the arguments, as
  loo does.

- log_ratios:

  For `psis()`, the draws object whose negative pointwise log-likelihood
  supplies the importance ratios.

## Value

A `loo`, `waic`, `compare.loo` or `psis` object from the loo package.

## Details

[`LOO()`](https://aforren1.github.io/frmtmb/reference/loo.html) and
[`WAIC()`](https://aforren1.github.io/frmtmb/reference/loo.html) are
brms's deprecated capitalized spellings and are defined only to name
their replacements.

## Priors, and what these numbers mean

These are posterior quantities, and they inherit the standing of the
draws they are computed from. Both of
[`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)'s
routes apply brms's default priors, so these numbers are regularized the
way brms regularizes them unless the call opted out with
`prior = "flat"`. Under that opt-out the elpd is likelihood-shaped and
unregularized: expect Pareto k warnings for models with many group-level
parameters, because a flat prior leaves those to be identified by the
data alone, and an influential observation then moves them a long way.
The maximum-likelihood answer to the same question is
[`AIC()`](https://rdrr.io/r/stats/AIC.html) on the fits, or
[`frmtmb::frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.html).

## Relative efficiency

`r_eff` defaults to
[`loo::relative_eff()`](https://mc-stan.org/loo/reference/relative_eff.html)
on the chain structure of the draws, which is what brms does. Thinning
with `ndraws` breaks that structure, so `r_eff` is then dropped and the
estimate is the one loo computes without an autocorrelation correction.

## See also

[`log_lik()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/log_lik.md),
[`frmtmb::bayes_R2()`](https://aforren1.github.io/frmtmb/reference/bayes_R2.html)

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    requireNamespace("loo", quietly = TRUE)) {
  set.seed(9)
  dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)))
  dd$y <- rnorm(60, 1 + 0.5 * dd$x + rnorm(6, 0, 0.5)[dd$g], 1)

  # sample with priors: an elpd is a posterior quantity
  d1 <- frm_sample(bf(y ~ x + (1 | g)), family = gaussian(),
                   data = dd, chains = 1, iter = 500, refresh = 0)
  d2 <- frm_sample(bf(y ~ 1 + (1 | g)), family = gaussian(),
                   data = dd, chains = 1, iter = 500, refresh = 0)
  loo_compare(d1, d2)
  # the same thing, one step at a time
  loo_compare(loo(d1), loo(d2))
}
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#> Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.
#>   model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>  model1       0.0     0.0      NA           3 k_psis > 0.58
#>  model2      -5.2     2.9    0.96   N < 100                
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
