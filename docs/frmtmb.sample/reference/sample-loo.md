# Approximate leave-one-out cross-validation

[`loo()`](https://mc-stan.org/loo/reference/loo.html) runs
Pareto-smoothed importance-sampling LOO and
[`waic()`](https://mc-stan.org/loo/reference/waic.html) the widely
applicable information criterion, both on the
[`frmtmb::log_lik()`](https://aforren1.github.io/frmtmb/reference/log_lik.html)
matrix, by handing it to
[`loo::loo.matrix()`](https://mc-stan.org/loo/reference/loo.html) and
[`loo::waic.matrix()`](https://mc-stan.org/loo/reference/waic.html)
unchanged. The returned objects are the loo package's own, so
[`print()`](https://rdrr.io/r/base/print.html) and
[`loo::pareto_k_table()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)
work on them directly.
[`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
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
psis(
  log_ratios,
  newdata = NULL,
  resp = NULL,
  model_name = NULL,
  ndraws = NULL,
  ...
)

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
  [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html))
  already-computed criteria.

- ndraws, resp:

  Passed to
  [`frmtmb::log_lik()`](https://aforren1.github.io/frmtmb/reference/log_lik.html).

- ...:

  Further models for
  [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html);
  otherwise passed to the loo package function.

- criterion:

  Which criterion
  [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
  computes for each draws object.

- model_names:

  Row names for the comparison; the default deparses the arguments, as
  loo does.

- log_ratios:

  For `psis()`, the draws object whose negative pointwise log-likelihood
  supplies the importance ratios.

- newdata:

  For `psis()`, accepted in brms's own second position and refused,
  because
  [`frmtmb::log_lik()`](https://aforren1.github.io/frmtmb/reference/log_lik.html)
  does not take it.

- model_name:

  For `psis()`, brms's label for the model. Accepted and unused: a
  `psis` object has nothing to label.

## Value

A `loo`, `waic`, `compare.loo` or `psis` object from the loo package.

## Details

[`LOO()`](https://paulbuerkner.com/brms/reference/loo.brmsfit.html) and
[`WAIC()`](https://paulbuerkner.com/brms/reference/waic.brmsfit.html)
are brms's deprecated capitalized spellings and are defined only to name
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

[`frmtmb::log_lik()`](https://aforren1.github.io/frmtmb/reference/log_lik.html),
[`frmtmb::bayes_R2()`](https://aforren1.github.io/frmtmb/reference/bayes_R2.html)

## Examples

``` r
# \donttest{
if (requireNamespace("tmbstan", quietly = TRUE) &&
    requireNamespace("rstan", quietly = TRUE) &&
    requireNamespace("loo", quietly = TRUE) &&
    !frmtmb.sample:::tmbstan_build_broken()) {
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
#>   sigma              student_t(3, 0, 2.5)  [natural scale]
#>   sd                 student_t(3, 0, 2.5)  [natural sd scale]
#>   b                  (flat), as brms leaves slopes
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
#>   Intercept          student_t(3, 0.8, 2.5)
#>   sigma              student_t(3, 0, 2.5)  [natural scale]
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
#>  model1       0.0     0.0      NA                          
#>  model2      -5.3     2.8    0.97   N < 100 2 k_psis > 0.58
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
# }
```
